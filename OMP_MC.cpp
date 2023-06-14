#include <iostream>
#include <stdlib.h>
#include <omp.h>
#include <random>
#include <cmath>

#define PI 3.141592653589793


//Generates n samples (Please enter n a multiple of thread_count :) )
void Box_Muller(int n, char* method, double* vals) {
    int m = omp_get_max_threads();
	int num_iter{ 0 };
	double* nums;
	if ("antithetic" == method) { int num_iter = n / (2 * m); }
	else { int num_iter = n / (2 * m); }
	std::mt19937 rng;
	rng.seed(1234);
	std::uniform_real_distribution<double> gen(0.0, 1.0);
#pragma omp parallel num_threads(m)
	nums = (double*) calloc(num_iter, sizeof(double));
	int tid = omp_get_thread_num();
	#pragma omp for
	for (int i = 0; i < num_iter; i++) {
		double u1{ 0.0 };
		double u2{ 0.0 };
		u1 = gen(rng);
		u2 = gen(rng);
		double r = sqrt(-2.0f * log(u1));
		nums[i] = r * cos(2.0 * PI * u2);
		nums[num_iter + i] = r * sin(2.0 * PI * u2);
		if ("antithetic" == method) { nums[2 * num_iter + i] = -r * cos(2.0f * PI * u2); nums[3 * num_iter + i] = -r * sin(2.0f * PI * u2); }
	}
	#pragma omp barrier
	#pragma omp critical
	{for (int i = 0; i < num_iter; i++) { vals[tid * num_iter + i] = nums[i]; }}
	free(nums);
}

int main()
{
	int maxThreads = omp_get_max_threads();
	int numThreads = omp_get_num_threads();
	std::cout << "Maximum number of threads: " << maxThreads << std::endl;
	std::cout << "Number of threads: " << numThreads << std::endl;

	return 0;
}
