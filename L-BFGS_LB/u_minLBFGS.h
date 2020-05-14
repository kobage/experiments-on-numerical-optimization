#pragma once			
#include<iostream>
using namespace std;
#include "auxiliaryFunctions.h" 

namespace u_min
{
	int n;
	double EPS;
	const int mVALUE = 5;	//a default value for parameter that determines the number of BFGS corrections saved
	double f0(0.0);			//objective function's value at x0 
	double f1(0.0);
	double  *x0, *dir, *tmp;
	double * g1, *g0;
	double * x1;
	double **s, **gamma;		//BFGS corrections to be saved 
	double *alpha_2loop, *ro_2loop;
	//functions
	void setDimentions(int n, double tolerance, int m);
	template<typename Functor, typename Problem, typename BreakType>
	void l_bfgs(Functor, Problem, BreakType, int m) noexcept;
	void L_BFGS_two_loop_recursion(int m);
	void destructLbfgsData();
};

void u_min::setDimentions(int k, double tolerance, int m)
{
	n = k;
	EPS = tolerance;
	dir = (double *)malloc(k * sizeof(double));
	x0 = (double *)malloc(k * sizeof(double));
	x1 = (double *)malloc(k * sizeof(double));
	g1 = (double *)malloc(k * sizeof(double));
	g0 = (double *)malloc(k * sizeof(double));
	s = (double **)malloc(m * sizeof(double*));
	gamma = (double **)malloc(m * sizeof(double*));
	for (int i = 0; i < m; ++i)
	{
		gamma[i] = (double *)malloc(k * sizeof(double));
		s[i] = (double *)malloc(k * sizeof(double));
	}
	alpha_2loop = (double *)malloc(m * sizeof(double));
	ro_2loop = (double *)malloc(m * sizeof(double));
}
void  u_min::destructLbfgsData()
{
	free(dir);
	free(x0);
	free(x1);
	free(g0);
	free(g1);
}

void  u_min::L_BFGS_two_loop_recursion(int m)
{
	vecCopy(g0, dir, n);
	//First loop
	for (int i = m - 1; i >-1; --i)
	{
		alpha_2loop[i] = ro_2loop[i] * vecProd(dir, s[i], n);
		for (int j = 0; j < n; ++j)
			dir[j] -= alpha_2loop[i] * gamma[i][j];
	}

	double coeff = 1 / (vecProd(gamma[m - 1], gamma[m - 1], n)* ro_2loop[m - 1]);
	for (int j = 0; j < n; ++j)
		dir[j] *= coeff;

	//Second loop
	double beta;
	for (int i = 0; i <m; ++i)
	{
		beta = ro_2loop[i] * vecProd(dir, gamma[i], n);
		for (int j = 0; j < n; ++j)
			dir[j] += (alpha_2loop[i] - beta) * s[i][j];
	}
}

template<typename Functor, typename Problem, typename BreakType>
void  u_min::l_bfgs(Functor lnSrchMethod, Problem functionName, BreakType breaking, int m) noexcept
{
	f0 = functionName(x0, g0, n);
	vecCopy(g0, dir, n);

	lnSrchMethod(functionName);

	for (int i = 0; i < n; ++i)
	{
		gamma[0][i] = g1[i] - g0[i];
		s[0][i] = x1[i] - x0[i];
	}
	ro_2loop[0] = 1 / vecProd(gamma[0], s[0], n);
	tmp = x0; x0 = x1; x1 = tmp;	//swap(x0,x1); 
	tmp = g0; g0 = g1; g1 = tmp;	//swap(g0,g1);

	int j = 1;
	while (j < m)
	{
		L_BFGS_two_loop_recursion(j);
		lnSrchMethod(functionName);

		for (int i = 0; i < n; ++i)
		{
			gamma[j][i] = g1[i] - g0[i];
			s[j][i] = x1[i] - x0[i];
		}
		ro_2loop[j] = 1 / vecProd(gamma[j], s[j], n);

		tmp = x0; x0 = x1; x1 = tmp;	//swap(x0,x1);
		tmp = g0; g0 = g1; g1 = tmp;	//swap(g0,g1);
		++j;
	}

	while (!breaking())
	{
		L_BFGS_two_loop_recursion(m);
		lnSrchMethod(functionName);

		tmp = gamma[0];
		for (int i = 0; i < m - 1; ++i)
			gamma[i] = gamma[i + 1];
		gamma[m - 1] = tmp;
		tmp = s[0];
		for (int i = 0; i < m - 1; ++i)
			s[i] = s[i + 1];
		s[m - 1] = tmp;

		for (int i = 0; i < m - 1; ++i)
			ro_2loop[i] = ro_2loop[i + 1];

		for (int i = 0; i < n; ++i)
		{
			gamma[m - 1][i] = g1[i] - g0[i];
			s[m - 1][i] = x1[i] - x0[i];
		}
		ro_2loop[m - 1] = 1 / vecProd(gamma[m - 1], s[m - 1], n);
		tmp = x0; x0 = x1; x1 = tmp;	//swap(x0,x1);
		tmp = g0; g0 = g1; g1 = tmp;	//swap(g0,g1);
	}
}