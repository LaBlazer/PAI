#pragma once

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <random>
#include <barrier>
#include <thread>

class DataMatrix {
	double * Data;
	uint32_t Size;

public:
	DataMatrix() = delete;

	DataMatrix(const DataMatrix& other) {
		if (this != &other)
		{
			Size = other.Size;
			if (Size != 0) {
				Data = new double[Size * Size];
				//std::cout << "AAAAAAAAAH, i'm COPYINGGGG " << &other << " TO " << this << std::endl;


				std::memcpy(Data, other.Data, Size * Size * sizeof(double));
			}
			else {
				Data = nullptr;
			}
		}
	}


	DataMatrix(uint32_t size) {
		//std::cout << "AAAAAAAAAH, i'm CONSTRUCTING " << this << std::endl;
		Size = size;

		Data = (size != 0) ? new double[size * size] : nullptr;
	}

	~DataMatrix() {
		//std::cout << "AAAAAAAAAH, i'm DESTROOOCTING " << this << std::endl;

		delete[] Data;
	}

	inline uint32_t size() { return Size; }

	inline double& at(int x, int y) { return Data[Size * y + x]; }
};

class Decomposer
{
	static inline std::unique_ptr<std::barrier<>> m_barrier;
	static void GaussianParallelJob(DataMatrix& mtx, int threadCount, int threadId);

public:
	static void PrintMatrix(DataMatrix& mtx);
	static bool CompareMatrix(DataMatrix& a, DataMatrix& b);
	static void GaussianSerial(DataMatrix& mtx);
	static void GaussianParallel(DataMatrix& mtx, int threadCount);
};