#pragma once
#include<iostream>
#include "phasePrimitives.h"
#include "lineSearch.h"
#include<algorithm>
#include<numeric>
#include <string>
using namespace std;

class l_bgfs
{
	long int& n;
	double& f0;			//objective function's value at x0 
	double& f1;
	double*& x0;
	double*& x1;
	double*& g1;
	double*& g0;
	double*& dir;
	double*& moment;
	long int& FuncGradEvaluations;
	
	long int LineSearchEvaluations{};
	double* tmp;
	int m;    //a value for parameter that determines the number of BFGS corrections saved
	double** s, ** gamma;		//BFGS corrections to be saved 
	double* alpha_2loop, * ro_2loop;

public:
	l_bgfs(phasePrimitives& prmt, int mValue)
		: n{ prmt.n }, f0{ prmt.f0 }, f1{ prmt.f1 }, x0{ prmt.x0 }, x1{ prmt.x1 }, g0{ prmt.g0 },m{ mValue },
		g1{ prmt.g1 }, dir{ prmt.dir }, moment{ prmt.moment }, FuncGradEvaluations{ prmt.FuncGradEvaluations } 
	{
		tmp = new double[n];
		s = (double**)malloc(m * sizeof(double*));
		gamma = (double**)malloc(m * sizeof(double*));
		for (int i = 0; i < m; ++i)
		{
			gamma[i] = (double*)malloc(n * sizeof(double));
			s[i] = (double*)malloc(n * sizeof(double));
		}
		alpha_2loop = (double*)malloc(m * sizeof(double));
		ro_2loop = (double*)malloc(m * sizeof(double));
	}
	~l_bgfs()
	{
		delete[] tmp;
		for (int i = 0; i < m; ++i)
		{
			gamma[i] = (double*)malloc(n * sizeof(double));
			s[i] = (double*)malloc(n * sizeof(double));
		}
	}


	long int getSize() const noexcept { return n; }
	long int getRestarts() const noexcept { return LineSearchEvaluations; }
	long int getEvaluations() const noexcept { return FuncGradEvaluations; }
	double* getVariables() const noexcept { return x0; }
	double* getGradient() const noexcept { return g0; }
	double getValue() const noexcept { return f0; }

	//solver
	template<class ProblemPtr>
	void solve(lineSearch& lnSrch, ProblemPtr ppr) noexcept;

	void L_BFGS_two_loop_recursion(int);

	void printStats()
	{
		std::cout << "Algorithm: L-BFGS,    Objective function:  " << f0 <<
			",   Value-Grad evaluations:  " << getEvaluations() <<
			",   Restarts (line searches):  " << getRestarts() << std::endl;
	}
};

void  l_bgfs::L_BFGS_two_loop_recursion(int j)
{
	std::copy(g0, g0+n,dir);
	//First loop
	for (int i = j - 1; i > -1; --i)
	{
		alpha_2loop[i] = ro_2loop[i] * std::inner_product(dir, dir + n, s[i], double{});
		for (int j = 0; j < n; ++j)
			dir[j] -= alpha_2loop[i] * gamma[i][j];
	}
	
	double coeff = 1 / (std::inner_product(gamma[j - 1], gamma[j - 1] + n, gamma[j - 1], double{}) * ro_2loop[j - 1]);
	for (int j = 0; j < n; ++j)
		dir[j] *= coeff;

	//Second loop
	double beta;
	for (int i = 0; i < j; ++i)
	{		
		beta = ro_2loop[i] * std::inner_product(dir, dir + n, gamma[i], double{});
		for (int j = 0; j < n; ++j)
			dir[j] += (alpha_2loop[i] - beta) * s[i][j];
	}
} 
template<class ProblemPtr>
void l_bgfs::solve(lineSearch& lnSrch, ProblemPtr ppr) noexcept
{

	f0 = ppr->valGrad(x0, g0);
	++FuncGradEvaluations;
	std::copy(g0, g0+n,dir);
	++LineSearchEvaluations;
	lnSrch(ppr);

	for (int i = 0; i < n; ++i)
	{
		gamma[0][i] = g1[i] - g0[i];
		s[0][i] = x1[i] - x0[i];
	}
	
	ro_2loop[0] = 1 / std::inner_product(s[0], s[0] + n, gamma[0], double{});
	swap(x0,x1); 
	swap(g0,g1);

	int j = 1;
	while (j < m)
	{
		L_BFGS_two_loop_recursion(j);
		++LineSearchEvaluations;/////
		lnSrch(ppr);

		for (int i = 0; i < n; ++i)
		{
			gamma[j][i] = g1[i] - g0[i];
			s[j][i] = x1[i] - x0[i];
		}
		ro_2loop[j] = 1 / std::inner_product(s[j], s[j] + n, gamma[j], double{});

		swap(x0,x1);
		swap(g0,g1);
		++j;
	}
	cout << "The number of nunction and gradient evaluations, before main loop:  " << FuncGradEvaluations << endl;   ///////

	while (!ppr->stoppingCondition(g0))
	{
		L_BFGS_two_loop_recursion(m);
		++LineSearchEvaluations;
		lnSrch(ppr);

		cout << "f=" << f0 << "   " << "|g| =" << infNorm(g0, n) << ",  evaluations:  " << FuncGradEvaluations << endl;///////////////////

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
		
		ro_2loop[m - 1] = 1 / std::inner_product(s[m - 1], s[m - 1] + n, gamma[m - 1], double{});
		swap(x0,x1);
		swap(g0,g1);
	}
}