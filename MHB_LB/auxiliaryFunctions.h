#pragma once
#include <iostream>
#include <random>
#include <cmath>

inline double infNorm
(
	double* x, /* vector */
	int     k /* length of vector */
)
{
	int i, n5;
	double t = 0.;
	n5 = k % 5;

	for (i = 0; i < n5; i++) if (t < fabs(x[i])) t = fabs(x[i]);
	for (; i < k; i += 5)
	{
		if (t < fabs(x[i])) t = fabs(x[i]);
		if (t < fabs(x[i + 1])) t = fabs(x[i + 1]);
		if (t < fabs(x[i + 2])) t = fabs(x[i + 2]);
		if (t < fabs(x[i + 3])) t = fabs(x[i + 3]);
		if (t < fabs(x[i + 4])) t = fabs(x[i + 4]);
	}
	return (t);
}

inline void vec2norm(double* ans, double* gr, const int k)
{
	int i, n5;
	*ans = 0.;
	n5 = k % 5;
	for (i = 0; i < n5; i++) *ans += gr[i] * gr[i];
	for (; i < k; i += 5)
	{
		*ans += gr[i] * gr[i] + gr[i + 1] * gr[i + 1] + gr[i + 2] * gr[i + 2]
			+ gr[i + 3] * gr[i + 3] + gr[i + 4] * gr[i + 4];
	}
	*ans = sqrt(*ans);
}

inline double vecProd(double* x, double* y, const int k)
{
	int i, n5;
	double sum(0.0);
	if (k <= 0) return sum;
	n5 = k % 5;
	for (i = 0; i < n5; i++) sum += x[i] * y[i];
	for (; i < k; i += 5)
	{
		sum += x[i] * y[i] + x[i + 1] * y[i + 1] + x[i + 2] * y[i + 2]
			+ x[i + 3] * y[i + 3] + x[i + 4] * y[i + 4];
	}
	return sum;
}

inline void vecCopy(double* x, double* y, const int n)
{
	int i, n5;
	if (n <= 0) return;
	n5 = n % 5;
	for (i = 0; i < n5; i++)  y[i] = x[i];
	while (i < n)
	{
		y[i] = x[i];
		++i;
		y[i] = x[i];
		++i;
		y[i] = x[i];
		++i;
		y[i] = x[i];
		++i;
		y[i] = x[i];
		++i;
	}
}

void printSolution(double value, double* sol, int dim)
{
	printf("  fx = %f\n\n", value);
	int quarter = dim / 4;
	for (int i = 0; i < quarter; i++)
	{
		printf("x[%3d] = %-15.11f     x[%3d] = %15.11f     x[%3d] = %15.11f     x[%3d] = %2.11f  \n",
			i, sol[i], quarter + i, sol[quarter + i], 2 * quarter + i,
			sol[2 * quarter + i], 3 * quarter + i, sol[3 * quarter + i]);
	}
	for (int i = 4 * quarter; i < dim; i++)
		printf("x[%3d] = %-15.11f     ", i, sol[i]);
	printf("\n");
}

void fillRandomMatrix(double** M, int k)
{
	std::default_random_engine dre;
	std::uniform_real_distribution<double> urdi(-1, 1);
	//std::normal_distribution<double> urdi(-1,1);
	//std::cauchy_distribution<double> urdi(-1, 1);

	for (int i = 0; i < k; i++)
	{
		for (int j = 0; j < i; j++)
			M[i][j] = urdi(dre);
		M[i][i] = 0;
	}
	for (int i = 0; i < k; i++)
		for (int j = i + 1; j < k; j++)
			M[i][j] = (-1.0 * M[j][i]);
}