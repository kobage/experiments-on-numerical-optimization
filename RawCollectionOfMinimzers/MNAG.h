#pragma once  
#include <string>
#include<algorithm>
#include"cm_namespace.h"

using namespace std;

class cm::MNAG
{
	double* moment;
	long int counter{};

public:
	MNAG(problem* p)
	{
		n = p->getSize();
		ppr = p;
		moment = new double[n];
		x0 = new double[n];
		x1 = new double[n];
		g1 = new double[n];
		g0 = new double[n];
	}
	~MNAG()
	{
		delete[]  moment;
		delete[]  x0;
		delete[]  x1;
		delete[]  g0;
		delete[]  g1;
	}

	//MNAG versions:
	template<class line_search>
	void operator()(line_search& , const double ) noexcept;
/*
	void solve(line_search&) noexcept;

	template<class line_search>
	void solve_BigBeta(line_search&) noexcept;

	template<class line_search>
	void solve_RandCorrection(line_search&) noexcept;
*/
	void printStats()
	{
		std::cout << "   Algorithm MNAG;  Objective function: " << f0 << ",  Value-Grad evaluations:  " << ValGradEvaluations << std::endl;
		std::cout << "   Value evaluations:  " << ValueEvaluations << ",  Restarts (line searches):  " << LineSearchCunter << std::endl;
	}


};

template<class line_search>
void cm::MNAG::operator()(line_search& lnSrch, const double BETA) noexcept
{
	while (true)
	{
		if (f0 <= f1 || counter == 100000)
		{
			if (1 != counter)
			{
				f0 = ppr->valGrad(x0, g0);
				++ValGradEvaluations;
			}
			if (ppr->stoppingCondition(g0))
				return;

			dir = g0;
			++LineSearchCunter;
			lnSrch();
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
		++ValGradEvaluations;
		if (ppr->stoppingCondition(g1))
		{
			swap(x0, x1);
			swap(g0, g1);
			return;
		}

		for (int i = 0; i < n; ++i)
			x1[i] = x0[i] + BETA * moment[i] - alpha * g1[i];
	

		f1 = ppr->value(x1);

		//  Statistics for ANN - uncomment of ANN 
		//	cout << "f=" << setw(10) << left <<  f1 << "   " << "|g| =" << setw(10) << left << infNorm(g1, n)
		//		<< "  evaluations:  " << ValGradEvaluations << endl;   //

		++ValueEvaluations;
		++counter;
	}
}