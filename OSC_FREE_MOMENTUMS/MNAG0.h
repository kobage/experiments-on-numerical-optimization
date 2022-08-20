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
	double* tmp;				// to maintain former x1 temporarily

	const double NAG_OR_HB{ 0.7 };   //from [0.5,1]: 0.5 gives pure MNAG, 1 gives pure MHB

	long int& FuncGradEvaluations;
	long int counter{};
	long int LineSearchCunter{};

	double alpha{ 0.001 };

	bool jumpToLnSrch{};
public:
	MNAG(phasePrimitives& prmt)
		//	: ppr{ prmt.ppr }, n{ prmt.n }, f0{ prmt.f0 }, f1{ prmt.f1 }, x0{ prmt.x0 }, x1{ prmt.x1 }, g0{ prmt.g0 },
		: n{ prmt.n }, f0{ prmt.f0 }, f1{ prmt.f1 }, x0{ prmt.x0 }, x1{ prmt.x1 }, g0{ prmt.g0 },
		g1{ prmt.g1 }, dir{ prmt.dir }, moment{ prmt.moment }, FuncGradEvaluations{ prmt.FuncGradEvaluations }	{tmp = new double[n]; }
	~MNAG() { delete[] tmp; }

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
	jumpToLnSrch = true;

	while (true)
	{
		if (jumpToLnSrch)
		{
			++LineSearchCunter;
			dir = g0;
			alpha = lnSrch(ppr);
			counter = 0;
			jumpToLnSrch = false;
		}

		if (ppr->stoppingCondition(g1))
		{
			swap(x0, x1);
			swap(g0, g1);
			return;
		}

		for (int i = 0; i < n; ++i)
			moment[i] = NAG_OR_HB * (x1[i] - x0[i]);

		for (int i = 0; i < n; ++i)
			x0[i] = NAG_OR_HB * x1[i] + (1 - NAG_OR_HB) * x0[i];
		f0 = f1;
		swap(g0, g1);
		std::copy(x1, x1 + n, tmp);

		//	cout << "f0=" << f0<<   "   " << "f1 =" << f1 << "   " <<  infNorm(g0,n)  << ",  evaluations:  " << FuncGradEvaluations << endl;///////////////////

		for (int i = 0; i < n; ++i)
			x1[i] = x0[i] + moment[i] - alpha * g0[i];
		f1 = ppr->valGrad(x1, g1);
		++FuncGradEvaluations;
		++counter;

		if (f0 <= f1 || counter == 100)
		{
			swap(x0, tmp);
			jumpToLnSrch = true;
		}
	}
}