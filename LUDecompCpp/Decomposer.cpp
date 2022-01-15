#include "Decomposer.h"

//Program::m_barrier;
//Barrier* Program::barrier;

void Decomposer::GaussianSerial(DataMatrix& mtx)
{
	const int n = mtx.size();

	for (int k = 0; k < n; k++)
	{
		//upper and lower decompositions are saved in one matrix
		// 
		//         u u u u u
		//         l u u u u
		//   mtx = l l u u u
		//         l l l u u
		//         l l l l u

		// Data\[(\w{1})\]\[(\w{1})\]

		for (unsigned int i = k + 1; i < n; i++)
		{
			double lik = mtx.at(i, k) / mtx.at(k, k);

			for (unsigned int j = k + 1; j < n; j++)
			{
				mtx.at(i, j) -= lik * mtx.at(k, j);
			}

			mtx.at(i, k) = lik;
		}
	}

	PrintMatrix(mtx);
}

void Decomposer::GaussianParallel(DataMatrix& mtx, int threadCount)
{
	m_barrier = std::make_unique<std::barrier<>>(threadCount);

	std::vector<std::thread> threads;
	threads.reserve(threadCount);

	for (int i = 0; i < threadCount; i++) {
		//std::cout << "running thread " << i << std::endl;
		threads.emplace_back(std::thread(GaussianParallelJob, std::ref(mtx), threadCount, i));
	}

	for (int i = 0; i < threadCount; i++) {
		threads[i].join();
		//std::cout << "thread joined " << i << std::endl;
	}
}

void Decomposer::GaussianParallelJob(DataMatrix& mtx, int threadCount, int threadId)
{
	const int n = mtx.size();
	threadId++;

	// 1-d row agglomeration
	for (unsigned int k = 0; k < n - 1; k++)
	{
		for (unsigned int i = k + threadId; i < n; i += threadCount)
		{
			// gaussian forward elimination
			const double lik = mtx.at(i, k) / mtx.at(k, k);

			for (int j = k + 1; j < n; j++)
			{
				mtx.at(i, j) = mtx.at(i, j) - lik * mtx.at(k, j);
			}

			mtx.at(i, k) = lik;
		}

		m_barrier->arrive_and_wait();
	}
}

void Decomposer::PrintMatrix(DataMatrix& mtx)
{
	if (mtx.size() < 100)
	{
		for (unsigned int i = 0; i < mtx.size(); i++)
		{
			for (unsigned int j = 0; j < mtx.size(); j++)
			{
				std::cout << mtx.at(i, j) << (((j + 1) < mtx.size()) ? ' ' : '\n');
			}
		}
	}
}



bool Decomposer::CompareMatrix(DataMatrix& a, DataMatrix& b)
{
	if (a.size() != b.size())
	{
		return false;
	}

	const int n = a.size();

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (a.at(i, j) != b.at(i, j))
			{
				return false;
			}
		}
	}

	return true;
}