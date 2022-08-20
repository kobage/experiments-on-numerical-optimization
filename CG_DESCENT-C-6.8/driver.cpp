#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "cg_user.h"
#include "cg_descent.h"
#include "cg_blas.h"
#include "testFunctions.h"
#include <iostream>
#include <chrono>
#include "game.h"

void run_cg_descent_Tests(struct testFunction* t, INT number)
{
	makeTestCollection(test);
	int repNumber = 10;
	int running_times[10];
	cg_stats Stats;

	puts("The following tests are available. plz enter test number");
	INT i;
	for (i = 0; i < number; i++)
		printf("%d. %s, size of test: %d\n ", i + 1, test[i].name, test[i].size);

	scanf("%d", &i);
	getchar();
	if (i < 1 || i > 45) return;
	i--;

	INT n = test[i].size;
	printf("Chosen: %s, size of test: %d\n ", test[i].name, test[i].size);

	double *x = (double *)malloc(n * sizeof(double));
	
	std::chrono::high_resolution_clock::time_point  st;  
	std::chrono::high_resolution_clock::duration  diff;  
	int j;
	for (j = 0; j < repNumber; j++)
	{
		test[i].StartingGess(x, n);
		cg_parameter Parm;
		cg_default(&Parm);
		//Parm.PrintLevel = 0 ;

		st = std::chrono::high_resolution_clock::now();

		cg_descent(x, n, &Stats, &Parm, 1.e-6, test[i].value, test[i].grad, test[i].valgrad, NULL);
		//cg_descent (x, n, &Stats, NULL, 1.e-6, test[i].value, test[i].grad, test[i].valgrad, NULL) ;
		diff = std::chrono::high_resolution_clock::now() - st;
		running_times[j] = std::chrono::duration_cast<std::chrono::microseconds>(diff).count();
	}
	iSort(running_times, 10);

	freopen("results.txt", "w", stdout);
	{
		std::cout << "Function: " << test[i].name << std::endl;
		printf("Median time in microseconds:   %d\n", running_times[repNumber / 2]);
		printf("Function's value:   %f\n", Stats.f);
		int quarter = n / 4;
		int ii;
		for (ii = 0; ii < quarter; ii++)
		{
			printf("x[%3d] = %-15.11f     x[%3d] = %15.11f     x[%3d] = %15.11f     x[%3d] = %2.11f  \n",
				ii, x[ii], quarter + ii, x[quarter + ii], 2 * quarter + ii, x[2 * quarter + ii], 3 * quarter + ii, x[3 * quarter + ii]);
		}
	}
	free(x);
}

#include<string>
#include<iostream>
#include<iomanip>
using namespace std;

int digits(long int k)
{
	if (k <= 0) return -1;
	int c = 0;
	while (k /= 10)
		++c;
	return c;
}

template<class T>
void printVector(T* vec, std::string s, long int n) noexcept
{
	long int quarter = n / 4;
	int sp = digits(n - 1) + 1;
	for (size_t i = 0; i < quarter; ++i)
	{
		std::cout << std::left << s << std::left << '[' << std::left << std::setw(sp) << i << std::left << "]="
			<< std::left << std::setw(21) << vec[i];
		std::cout << std::left << s << std::left << '[' << std::left << std::setw(sp) << i + quarter << std::left << "]="
			<< std::left << std::setw(21) << vec[i + quarter];
		std::cout << std::left << s << std::left << '[' << std::left << std::setw(sp) << i + 2 * quarter << std::left << "]="
			<< std::left << std::setw(21) << vec[i + 2 * quarter];
		std::cout << std::left << s << std::left << '[' << std::left << std::setw(sp) << i + 3 * quarter << std::left << "]="
			<< std::left << std::setw(21) << vec[i + 3 * quarter];
		std::cout << std::endl;
	}
	for (int i = 4 * quarter; i < n; i++)
		std::cout << std::left << s << std::left << '[' << std::left << std::setw(sp) << i << std::left << "]="
		<< std::left << std::setw(21) << vec[i];
	std::cout << std::endl;
}

//--------------------------------  Benchmarking ANN  ---------------------------------------
#include<iostream>
#include<iomanip>
#include<vector>
#include<algorithm>
#include<chrono>
#include"log_regression.h"
using namespace std;

void runANNTest()
{
	cg_stats Stats;  
	long int n = a->n;
	double* x = new double[n];
//	freopen("results.txt", "w", stdout);
	
	size_t j = 1;									 //runs for the specific neuron
	{
		auto st = chrono::high_resolution_clock::now();
		a->set_y(j);
		cg_parameter Parm;
		cg_default(&Parm);
		Parm.PrintLevel = 1;       //for less information take: Parm.PrintLevel = 0;
		
		Initial_ANN(x);
		//The precision is set on line 1694 in file "cg_descent.cppp", so 3rd parameter does not matter
		cg_descent(x, n, &Stats, &Parm, 1.e-6, ANNVal, ANNGrad, ANNValGrad, NULL);


		auto diff = chrono::high_resolution_clock::now() - st;
		auto time = chrono::duration_cast<chrono::microseconds>(diff);

		std::cout << "Neuron " << j << ":  " << endl;
		std::cout << "Problem:  " << a->name << ",  neuron " << j << std::endl;
		std::cout << "Time:  " << time.count() << " microseconds" << std::endl;
		printVector(x, "x", n);

		std::cout << std::endl << std::endl;
	}
	delete[] x;
	delete a;
}


int main(void)
{
	//uncommnet line 1695 in file "cg_descent.cppp"
//	run_cg_descent_Tests(test, TestNumber);  //TestNumber - size of the array of the test-functions
	
	//uncommnet line 1694 in file "cg_descent.cppp"
	runANNTest();  

	//uncommnet line 1696 in file "cg_descent.cppp"
//	runSymGameTest();
}