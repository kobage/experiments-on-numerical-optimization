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
	long int ValueEvaluations{};
	long int counter{};
	long int LineSearchCunter{};
	double alpha{0.001};

	MpSpMtr hessian; 
public:
	MNAG(phasePrimitives& prmt)
//	: ppr{ prmt.ppr }, n{ prmt.n }, f0{ prmt.f0 }, f1{ prmt.f1 }, x0{ prmt.x0 }, x1{ prmt.x1 }, g0{ prmt.g0 },
	:  n{ prmt.n }, f0{ prmt.f0 }, f1{ prmt.f1 }, x0{ prmt.x0 }, x1{ prmt.x1 }, g0{ prmt.g0 },
		g1{ prmt.g1 }, dir{ prmt.dir }, moment{ prmt.moment }, FuncGradEvaluations{ prmt.FuncGradEvaluations }	
	{
		hessian.matrix.resize(n);
		hessian.LLT.resize(n);
		hessian.n = n;
	}

	~MNAG() {}

	long int getSize() const noexcept { return n; }
	long int getRestarts() const noexcept { return LineSearchCunter; }
	long int getEvaluations() const noexcept { return FuncGradEvaluations; }
	double* getVariables() const noexcept { return x0; }
	double* getGradient() const noexcept { return g0; }
	double getValue() const noexcept { return f0; }
	
	//solvers
	template<class ProblemPtr>
	void solve(lineSearch& lnSrch, ProblemPtr ppr) noexcept;

	template<class ProblemPtr>
	void solve_BigBeta(lineSearch& lnSrch, ProblemPtr ppr) noexcept;

	template<class ProblemPtr>
	void solve_RandCorrection(lineSearch& lnSrch, ProblemPtr ppr) noexcept;

	void printStats()
	{
		std::cout << "Algorithm: MNAG,    Objective function:  " << f0 <<
			",   Value-Grad evaluations:  " << getEvaluations() << ",  " << std::endl
			<< "Values evaluations : " << ValueEvaluations <<
			",   Restarts (line searches):  " << getRestarts() << std::endl;
	}
};

#include<random>
#include<iomanip>
template<class ProblemPtr>
void  MNAG::solve(lineSearch& lnSrch, ProblemPtr ppr) noexcept
{
	while (true)
	{
		if (f0 <= f1 || counter == 100000)
		{
			if (1 != counter)
			{
				f0 = ppr->valGrad(x0, g0);
				++FuncGradEvaluations;
			}
			if (ppr->stoppingCondition(g0))
				return;

			dir = g0;
			++LineSearchCunter;
			alpha = lnSrch(ppr);
			counter = 0;
		}

		for (int i = 0; i < n; ++i)
			moment[i] = x1[i] - x0[i];
		swap(x0, x1);
		swap(g0, g1);
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
			x1[i] = x0[i] + moment[i] - alpha * g1[i];

		f1 = ppr->value(x1);

	//  Statistics for ANN
	//	cout << "f=" << setw(10) << left <<  f1 << "   " << "|g| =" << setw(10) << left << infNorm(g1, n)
	//		<< "  evaluations:  " << FuncGradEvaluations << endl;   //

		++ValueEvaluations;
		++counter;
	}
}

#include<random>
#include<iomanip>
template<class ProblemPtr>
void  MNAG::solve_RandCorrection(lineSearch& lnSrch, ProblemPtr ppr) noexcept
 {
	default_random_engine dre(std::random_device{}());
	uniform_real_distribution<double> di(0.5, 1.);
	double correction_rate{};

	while (true)
	{
		if (f0 <= f1 || counter == 100000)
		{
			if (1 != counter)
			{
				f0 = ppr->valGrad(x0, g0);
				++FuncGradEvaluations;
			}
			if (ppr->stoppingCondition(g0))
				return;

			dir = g0;
			++LineSearchCunter;
			alpha = lnSrch(ppr);
			counter = 0;
		}

		for (int i = 0; i < n; ++i)
			moment[i] = x1[i] - x0[i];
		swap(x0, x1);
		swap(g0, g1);
		f0 = f1;

		correction_rate = di(dre);
		for (int i = 0; i < n; ++i)
			x1[i] = x0[i] + correction_rate * moment[i];

		ppr->valGrad(x1, g1);
		++FuncGradEvaluations;
		if (ppr->stoppingCondition(g1))
		{
			swap(x0, x1);
			swap(g0, g1);
			return;
		}

		for (int i = 0; i < n; ++i)
			x1[i] = x0[i] + moment[i] - alpha * g1[i];
		++counter;

		f1 = ppr->value(x1);

		//  Statistics for ANN
		// cout << "f=" << setw(10) << left <<  f1 << "   " << "|g| =" << setw(10) << left << infNorm(g1, n)
		//	<< "  evaluations:  " << FuncGradEvaluations << endl;   //

		++ValueEvaluations;
		++counter;
	}
}

template<class ProblemPtr>
void  MNAG::solve_BigBeta(lineSearch& lnSrch, ProblemPtr ppr) noexcept
{
	while (true)
	{
		if (f0 <= f1 || counter == 100000)
		{
			if (1 != counter)
			{
				f0 = ppr->valGrad(x0, g0);
				++FuncGradEvaluations;
			}
			if (ppr->stoppingCondition(g0))
				return;

			dir = g0;
			++LineSearchCunter;
			alpha = lnSrch(ppr);
			counter = 0;
		}

		for (int i = 0; i < n; ++i)
			moment[i] = x1[i] - x0[i];
		swap(x0, x1);
		swap(g0, g1);
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
			x1[i] = x0[i] + 1.0005 * moment[i] - alpha * g1[i];

		f1 = ppr->value(x1);

		//  Statistics for ANN
		//	cout << "f=" << setw(10) << left <<  f1 << "   " << "|g| =" << setw(10) << left << infNorm(g1, n)
		//		<< "  evaluations:  " << FuncGradEvaluations << endl;   //

		++ValueEvaluations;
		++counter;
	}
}