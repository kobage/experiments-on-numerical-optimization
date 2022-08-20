#pragma once  
#include <string>
#include"problems.h"
#include<algorithm>
#include "phasePrimitives.h"
#include "lineSearch.h"
using namespace std;

class MNAG
{
	long int& n;
	double& f0;
	double& f1;
	double*& x0;
	double*& x1;
	double*& g1;
	double*& g0;
	double*& dir;
	double*& moment;

	long int& FuncGradEvaluations;
	long int counter{};
	long int LineSearchCunter{};

	double alpha{0.001};

public:
	MNAG(phasePrimitives& prmt)
//	: ppr{ prmt.ppr }, n{ prmt.n }, f0{ prmt.f0 }, f1{ prmt.f1 }, x0{ prmt.x0 }, x1{ prmt.x1 }, g0{ prmt.g0 },
	:  n{ prmt.n }, f0{ prmt.f0 }, f1{ prmt.f1 }, x0{ prmt.x0 }, x1{ prmt.x1 }, g0{ prmt.g0 },
		g1{ prmt.g1 }, dir{ prmt.dir }, moment{ prmt.moment }, FuncGradEvaluations{ prmt.FuncGradEvaluations }	{ }
	~MNAG() {}

	long int getSize() const noexcept { return n; }
	long int getRestarts() const noexcept { return LineSearchCunter; }
	long int getEvaluations() const noexcept { return FuncGradEvaluations; }
	double* getVariables() const noexcept { return x0; }
	double* getGradient() const noexcept { return g0; }
	double getValue() const noexcept { return f0; }
	
	//solver
	template<class ProblemPtr>
	void solve(lineSearch& lnSrch, ProblemPtr ppr) noexcept;

	void printStats()
	{
		std::cout << "Algorithm: MNAG,    Objective function:  " << f0 <<
			",   Value-Grad evaluations:  " << getEvaluations() <<
			",   Restarts (line searches):  " << getRestarts() << std::endl;
	}
};
template<class ProblemPtr>
void  MNAG::solve(lineSearch& lnSrch, ProblemPtr ppr) noexcept
{
	FuncGradEvaluations = 0;
	LineSearchCunter = 0;

	f0 = f1 = ppr->valGrad(x0, g0);
	++FuncGradEvaluations;
	if (ppr->stoppingCondition(g0))
		return;

	while (true)
	{
		if (f0 <= f1 || counter == 100000)
		{
			dir = g0;
			++LineSearchCunter;
			alpha = lnSrch(ppr);
			counter = 0;
		}

		for (int i = 0; i < n; ++i)
			moment[i] = x1[i] - x0[i];
		swap(x0, x1);
		swap(g0, g1);

		if (ppr->stoppingCondition(g0))
			return;

		f0 = f1;
		for (int i = 0; i < n; ++i)
			x1[i] = x0[i] + moment[i];

		ppr->valGrad(x1, g1);
		++FuncGradEvaluations;
		if (ppr->stoppingCondition(g1))
		{
			swap(x0, x1);
			swap(g0, g1);
			return;
		}

		for (int i = 0; i < n; ++i)
			x1[i] -=  (alpha * g1[i]);

		cout << "f0=" << f0<<   "   " << "|g| =" <<  infNorm(g0,n)  << ",  evaluations:  " << FuncGradEvaluations << endl;   //

		f1 = ppr->valGrad(x1, g1);
		++FuncGradEvaluations;
		++counter;
	}
}