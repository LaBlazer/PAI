
using MPI;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.Threading;
using System.Threading.Tasks;

namespace LUDecomp
{
    class Program
    {

        static readonly object llll = new object();

        static readonly Random rnd = new Random();

        static ConcurrentQueue<double[]>[] m_broadcastQueue;
        static Thread[] m_threads;
        static int m_threadDone;

        Intracommunicator world = Communicator.world;

        static void DoolittleSerial(double[,] mtx)
        {
            // Source: https://ieeexplore.ieee.org/document/7154772

            if(mtx.GetLength(0) != mtx.GetLength(1))
            {
                Console.WriteLine("Not a block matrix!");
                return;
            }

            int n = mtx.GetLength(0);
            double sum = 0;
            double[,] lower = new double[n, n];
            double[,] upper = new double[n, n];

            for (int i = 0; i < n; i++)
            {
                lower[i, i] = 1; // Diagonal 1

                // Calculate upper
                for (int k = i; k < n; k++)
                {
                    // Sum L(i, j) * U(j, k)
                    sum = 0;
                    for (int j = 0; j < i; j++)
                        sum += (lower[i, j] * upper[j, k]);

                    upper[i, k] = mtx[i, k] - sum;
                }

                // Calculate lower
                for (int k = i; k < n; k++)
                {
                    // Sum L(k, j) * U(j, i)
                    sum = 0;
                    for (int j = 0; j < i; j++)
                        sum += (lower[k, j] * upper[j, i]);

                    lower[k, i] = (mtx[k, i] - sum) / upper[i, i];
                }
            }

            Console.WriteLine("Upper mtx:");
            PrintMatrix(upper);

            Console.WriteLine("Lower mtx:");
            PrintMatrix(lower);
        }

        static void GaussianSerial(double[,] mtx)
        {
            if (mtx.GetLength(0) != mtx.GetLength(1))
            {
                Console.WriteLine("Not a block matrix!");
                return;
            }

            int n = mtx.GetLength(0);

            for (int k = 0; k < n; k++)
            {
                //upper and lower decompositions are saved in one matrix
                // 
                //         u u u u u
                //         l u u u u
                //   mtx = l l u u u
                //         l l l u u
                //         l l l l u

                for (int i = k + 1; i < n; i++)
                {
                    double lik = mtx[i, k] / mtx[k, k];

                    for (int j = k; j < n; j++)
                        mtx[i, j] -= lik * mtx[k, j];

                    mtx[i, k] = lik;
                }
            }

            PrintMatrix(mtx);

            double[,] lower, upper;
            SplitLU(mtx, out lower, out upper);
            //Console.WriteLine("============ Upper mtx ===========");
            //PrintMatrix(upper);
            //Console.WriteLine("============ Lower mtx ===========");
            //PrintMatrix(lower);
        }

        static void Broadcast(ref double[] data, int source, int root)
        {
            if (source == root)
            {
                // broadcast
                for (int i = 0; i < m_broadcastQueue.Length; i++)
                {
                    if (i != source)
                        m_broadcastQueue[source].Enqueue(data);
                }
            } else
            {
                if (m_broadcastQueue[source].Count > 0)
                {
                    // else receive
                    m_broadcastQueue[source].TryDequeue(out data);
                }
            }
        }

        static IEnumerable<double[]> Receive(int myId)
        {
            double[] data;
            while (m_broadcastQueue[myId].TryDequeue(out data))
            {
                yield return data;
            }
        }

        static bool AreThreadsRunning(int tid)
        {
            for (int i = 0; i < m_broadcastQueue.Length; i++)
            {
                if (i != tid && m_broadcastQueue[i].Count > 0)
                    return true;
            }

            return false;
        }

        static void GaussianParallelJob(int threadId, int threadCount, double[,] mtx)
        {
            Console.WriteLine($"running job {threadId + 1}/{threadCount}");

            int n = mtx.GetLength(0);

            int counter = 0;
            double[] buffer = new double[n];

            // 1-d column agglomeration
            for (int k = 0; k < n - 1; k++)
            {             
                for (int i = k + 1; i < n; i++)
                {
                    if ((i % threadCount) == threadId)
                    {
                        // gaussian forward elimination
                        counter++;

                        double lik = mtx[i, k] / mtx[k, k];

                        for (int j = k + 1; j < n; j++)
                            mtx[i, j] -= lik * mtx[k, j];

                        mtx[i, k] = lik;
                    }
                }

                for (int i = k + 1; i < n; i++)
                {
                    // copy buffer
                    for (int j = k; j < n; j++)
                        buffer[j] = mtx[i, j];

                    Console.WriteLine($"thread {threadId} before broadcast {string.Join(", ", buffer)}");

                    Broadcast(ref buffer, threadId, i % threadCount);

                    Console.WriteLine($"thread {threadId} after broadcast {string.Join(", ", buffer)}");

                    // apply buffer
                    for (int j = k; j < n; j++)
                        mtx[i, j] = buffer[j];
                }
            }

            m_threadDone++;

            while (m_threadDone < m_threads.Length) { }

            Console.WriteLine($"thread {threadId} did {counter} ops");

            //if(threadId == 0)
            //{
            Thread.Sleep(500 * threadId);
            PrintMatrix(mtx);
            Console.WriteLine();
            //}

            return;
        }

        static void GaussianParallel(double[,] mtx, int threadCount)
        {
            // Source: https://relate.cs.illinois.edu/course/cs554-f21/f/slides/slides_06.pdf

            if (mtx.GetLength(0) != mtx.GetLength(1))
            {
                Console.WriteLine("Not a block matrix!");
                return;
            }

            int n = mtx.GetLength(0);
            threadCount = Math.Min(threadCount, n);


            m_broadcastQueue = new ConcurrentQueue<double[]>[threadCount];
            m_threads = new Thread[threadCount];

            for(int i = 0; i < threadCount; i++)
            {
                m_broadcastQueue[i] = new ConcurrentQueue<double[]>();

                int idx = i;
                m_threads[i] = new Thread((i) => GaussianParallelJob(idx, threadCount, (double[,])mtx.Clone()));
                
            }

            for (int i = 0; i < threadCount; i++)
            {
                m_threads[i].Start();
            }

            for (int i = 0; i < threadCount; i++)
            {
                m_threads[i].Join();
            }

            //PrintMatrix(mtx);
        }

        static void SplitLU(double[,] mtx, out double[,] lower, out double[,] upper)
        {
            int n = mtx.GetLength(0);
            lower = new double[n, n];
            upper = new double[n, n];


            for (int i = 0; i < n; i++)
            {
                lower[i, i] = 1;

                for (int j = 0; j < i; j++)
                {
                    lower[i, j] = mtx[i, j];
                }

                for (int j = i; j < n; j++)
                {
                    upper[i, j] = mtx[i, j];
                }
            }
        }

        static void PrintMatrix(double[,] mtx)
        {
            if (mtx.Length < 100)
            {
                for (int i = 0; i < mtx.GetLength(0); i++)
                {
                    for (int j = 0; j < mtx.GetLength(1); j++)
                        Console.Write(string.Format("{0,18:F14}", mtx[i, j]) + (((j + 1) < mtx.GetLength(1)) ? " " : "\n"));
                }
            }
        }

        static double[,] GenerateRandomMatrix(int n, double min = 0, double max = 1)
        {
            double[,] mtx = new double[n, n];

            for(int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    mtx[i, j] = rnd.NextDouble() * (max - min) + min;
                }
            }

            return mtx;
        }

        public static void Main(String[] arr)
        {
            Stopwatch sw = new Stopwatch();

            double[,] mtx = { { 1, 2, 3 },
                        { 4, 5, 3 },
                        { 1, 2, 6 } };

            double[,] mtx2 = GenerateRandomMatrix(5, -10, 10);

            Console.WriteLine("Input mtx:");
            //PrintMatrix(mtx2);

            //Console.WriteLine("============\nDoolittle serial result: ");
            //sw.Restart();
            //DoolittleSerial(mtx2);
            //sw.Stop();
            //Console.WriteLine($"Time taken: {sw.ElapsedMilliseconds}ms");

            Console.WriteLine("============\nGaussian serial result: ");
            sw.Restart();
            GaussianSerial((double[,])mtx2.Clone());
            sw.Stop();
            Console.WriteLine($"Time taken: {sw.ElapsedMilliseconds}ms");


            Console.WriteLine("============\nParallel algorithm result: ");
            sw.Restart();
            GaussianParallel((double[,])mtx2.Clone(), 4);
            sw.Stop();
            Console.WriteLine($"Time taken: {sw.ElapsedMilliseconds}ms");

            //double[,] lower, upper;
            //SplitLU(mtx, out lower, out upper);
            //Console.WriteLine("============ Upper");
            //PrintMatrix(upper);
            //Console.WriteLine("============ Lower");
            //PrintMatrix(lower);
            
        }
    }
}
