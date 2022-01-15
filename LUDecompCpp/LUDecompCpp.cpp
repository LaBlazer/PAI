#include <iostream>
#include <fstream>
#include <random>
#include <chrono>

using namespace std::chrono;

#include "Decomposer.h"

DataMatrix LoadMatrix(const std::string& filename)
{
	std::string line, token;
	std::ifstream file(filename);

	if (file.is_open())
	{
		// get mtx size
		std::getline(file, line);
		const size_t n = std::count(line.begin(), line.end(), ' ') + 1;

		DataMatrix dm(n);

		unsigned int i = 0, j = 0;
		size_t start = 0, end = 0;
		do
		{
			j = 0;
			start = 0;
			end = line.find(' ', start);

			while (end != std::string::npos) {
				dm.at(i, j) = std::stod(line.substr(start, end - start));

				//std::cout << std::stod(line.substr(start, end - start)) << std::endl;
				//std::cout << dm.at(i, j) << std::endl;

				start = end + 1;
				end = line.find(' ', start);

				j++;
			}

			dm.at(i, j) = std::stod(line.substr(start));

			i++;
		} while (std::getline(file, line) && i < n);



		file.close();

		return dm;
	}

	return DataMatrix(0);
}

int main(int argc, char* argv[])
{
	double speedSerial, speedParallel;

	if (argc < 3)
	{
		std::cout << "Bad args!\nUsage: input.txt <proc count>" << std::endl;
		return 0;
	}

	std::cout << "Loading from " << argv[1] << std::endl;
	auto dm = LoadMatrix(std::string(argv[1]));
	
	std::cout << "Input mtx size " << dm.size() << std::endl;
	Decomposer::PrintMatrix(dm);

	std::cout << "===========" << std::endl << "Gaussian serial" << std::endl;
	DataMatrix dmSerial(dm);

	std::cout << std::endl;

	auto start = high_resolution_clock::now();
	Decomposer::GaussianSerial(dmSerial);
	auto stop = high_resolution_clock::now();
	speedSerial = duration_cast<milliseconds>(stop - start).count();

	std::cout << "Time taken: " << speedSerial << "ms" << std::endl;

	std::cout << std::endl;

	std::cout << "===========" << std::endl << "Gaussian parallel" << std::endl;

	start = high_resolution_clock::now();
	Decomposer::GaussianParallel(dm, std::stoi(argv[2]));
	stop = high_resolution_clock::now();
	speedParallel = duration_cast<milliseconds>(stop - start).count();

	Decomposer::PrintMatrix(dm);

	std::cout << "Time taken: " << speedParallel << "ms" << std::endl << std::endl;

	std::cout << "Speedup: " << speedSerial/speedParallel << "x" << std::endl;

	
}
