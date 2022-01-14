
using MPI;
using System;
using System.Diagnostics;
using System.IO;

namespace LUDecomp
{
    static class ExtensionMethods
    {
        public static T[][] ToJaggedArray<T>(this T[,] twoDimensionalArray)
        {
            int rowsFirstIndex = twoDimensionalArray.GetLowerBound(0);
            int rowsLastIndex = twoDimensionalArray.GetUpperBound(0);
            int numberOfRows = rowsLastIndex - rowsFirstIndex + 1;

            int columnsFirstIndex = twoDimensionalArray.GetLowerBound(1);
            int columnsLastIndex = twoDimensionalArray.GetUpperBound(1);
            int numberOfColumns = columnsLastIndex - columnsFirstIndex + 1;

            T[][] jaggedArray = new T[numberOfRows][];
            for (int i = 0; i < numberOfRows; i++)
            {
                jaggedArray[i] = new T[numberOfColumns];

                for (int j = 0; j < numberOfColumns; j++)
                {
                    jaggedArray[i][j] = twoDimensionalArray[i + rowsFirstIndex, j + columnsFirstIndex];
                }
            }
            return jaggedArray;
        }

        public static T[][] CloneJaggedArray<T>(T[][] source)
        {
            var len = source.Length;
            var dest = new T[len][];

            for (var x = 0; x < len; x++)
            {
                var inner = source[x];
                var ilen = inner.Length;
                var newer = new T[ilen];
                Array.Copy(inner, newer, ilen);
                dest[x] = newer;
            }

            return dest;
        }
    }

    class Program
    {

        static readonly Random rnd = new Random();

        static void DoolittleSerial(double[,] mtx)
        {
            // Source: https://ieeexplore.ieee.org/document/7154772

            if (mtx.GetLength(0) != mtx.GetLength(1))
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

        static void GaussianSerial(double[][] mtx)
        {
            if (mtx.Length != mtx[0].Length)
            {
                Console.WriteLine("Not a block matrix!");
                return;
            }

            int n = mtx.Length;

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
                    double lik = mtx[i][k] / mtx[k][k];

                    for (int j = k + 1; j < n; j++)
                        mtx[i][j] -= lik * mtx[k][j];

                    mtx[i][k] = lik;
                }
            }

            PrintMatrix(mtx);

            //double[,] lower, upper;
            //SplitLU(mtx, out lower, out upper);
            //Console.WriteLine("============ Upper mtx ===========");
            //PrintMatrix(upper);
            //Console.WriteLine("============ Lower mtx ===========");
            //PrintMatrix(lower);
        }

        static void GaussianParallelJob(double[][] mtx)
        {
            // Source https://relate.cs.illinois.edu/course/cs554-f21/f/slides/slides_06.pdf

            Intracommunicator world = Communicator.world;
            int threadCount = world.Size;
            int threadId = world.Rank;

            //Console.WriteLine($"running job {threadId + 1}/{threadCount}");

            int n = mtx.Length;

            world.Barrier(); // wait for others

            // 1-d row agglomeration
            for (int k = 0; k < n - 1; k++)
            {
                for (int i = k + 1; i < n; i++)
                {
                    if ((i % threadCount) == threadId)
                    {
                        // gaussian forward elimination
                        double lik = mtx[i][k] / mtx[k][k];

                        for (int j = k + 1; j < n; j++)
                            mtx[i][j] -= lik * mtx[k][j];

                        mtx[i][k] = lik;
                    }
                }

                for (int i = k + 1; i < n; i++)
                {
                    world.Broadcast(ref mtx[i], i % threadCount);
                }
            }

            return;
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

        static void PrintMatrix(double[][] mtx)
        {
            if (mtx.Length < 100)
            {
                for (int i = 0; i < mtx.Length; i++)
                {
                    for (int j = 0; j < mtx.Length; j++)
                        Console.Write(string.Format("{0,18:F14}", mtx[i][j]) + (((j + 1) < mtx.Length) ? " " : "\n"));
                }
            }
        }

        static double[,] GenerateRandomMatrix(int n, double min = 0, double max = 1)
        {
            double[,] mtx = new double[n, n];

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    mtx[i, j] = rnd.NextDouble() * (max - min) + min;
                }
            }

            return mtx;
        }

        static void SaveMatrix(double[,] mtx, string filename)
        {
            int n = mtx.GetLength(0);
            System.Text.StringBuilder sb = new System.Text.StringBuilder(n * n * 15);

            for (int i = 0; i < n; i++)
            {
                int j;
                for (j = 0; j < n - 1; j++)
                {
                    sb.Append(mtx[i, j]);
                    sb.Append(' ');
                }

                sb.Append(mtx[i, j]);
                sb.Append('\n');
            }

            File.WriteAllText(filename, sb.ToString());
        }

        static double[,] LoadMatrix(string filename)
        {
            string[] lines = File.ReadAllLines(filename);
            int n = lines[0].Split(' ').Length;
            double[,] mtx = new double[n, n];

            for (int i = 0; i < lines.Length; i++)
            {
                string[] nums = lines[i].Split(' ');
                for (int j = 0; j < nums.Length; j++)
                {
                    if (!double.TryParse(nums[j], out mtx[i, j]))
                    {
                        Console.WriteLine("Error while parsing!");
                        return mtx;
                    }

                }
            }

            return mtx;
        }

        static bool CompareMatrix(double[][] a, double[][] b)
        {
            if (a.Length != b.Length || a[0].Length != b[0].Length)
                return false;

            int n = a.Length;

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (a[i][j] != b[i][j])
                    {
                        return false;
                    }
                }
            }

            return true;
        }

        public static void Main(String[] args)
        {
            //Console.WriteLine(string.Join(" + ", args));

            if (args.Length >= 3)
            {
                if (args[0] == "-g")
                {
                    Console.WriteLine($"Generating matrix of size {args[1]} to '{args[2]}'");

                    double[,] mtx = GenerateRandomMatrix(int.Parse(args[1]), -50, 50);

                    Console.WriteLine("Saving...");
                    SaveMatrix(mtx, args[2]);

                    return;
                }
            }

            using (new MPI.Environment(ref args))
            {
                Stopwatch sw = new Stopwatch();

                if (args.Length < 1)
                {
                    Console.WriteLine("Bad args! Usage: input.txt");
                    return;
                }

                double[][] mtx = ExtensionMethods.ToJaggedArray(LoadMatrix(args[0]));

                if (Intercommunicator.world.Rank == 0) // if we are root process
                {
                    Console.WriteLine($"Loading from {args[0]}");
                    Console.WriteLine("Input mtx");
                    PrintMatrix(mtx);

                    Console.WriteLine("============\nGaussian serial algorithm");
                    sw.Restart();
                    double[][] serialMtx = ExtensionMethods.CloneJaggedArray(mtx);
                    GaussianSerial(serialMtx);
                    sw.Stop();
                    Console.WriteLine($"Time taken: {sw.ElapsedMilliseconds}ms");


                    Console.WriteLine("============\nGaussian parallel algorithm");

                    Console.WriteLine("MPI world size: " + Communicator.world.Size);


                    sw.Restart();
                    GaussianParallelJob(mtx);
                    sw.Stop();

                    Console.WriteLine($"Time taken: {sw.ElapsedMilliseconds}ms");
                    Console.WriteLine();

                    PrintMatrix(mtx);

                    if (CompareMatrix(serialMtx, mtx))
                    {
                        Console.WriteLine("Result is correct");
                    }
                    else
                    {
                        Console.WriteLine("Result is not correct");
                    }
                }
                else
                {
                    GaussianParallelJob(mtx);
                }
            }
        }
    }
}
