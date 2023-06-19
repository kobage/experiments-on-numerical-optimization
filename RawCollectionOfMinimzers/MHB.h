#pragma once  
#include <string>   
#include<algorithm>
#include"cm_namespace.h"

using namespace std;

class cm::MHB
{
	double* moment;
	long int counter{};
public:
	MHB(problem* p)
	{
		n = p->getSize();
		ppr = p;
		moment = new double[n];
		x0 = new double[n];
		x1 = new double[n];
		g1 = new double[n];
		g0 = new double[n];
	}
	~MHB()
	{
		delete[]  moment;
		delete[]  x0;
		delete[]  x1;
		delete[]  g0;
		delete[]  g1;
	}

	template<typename line_search>
	void operator()(line_search&, const double) noexcept;

	void printStats()
	{
		std::cout << "   Algorithm MHB;  Objective function: " << f0 << ",  Value-Grad evaluations:  " << ValGradEvaluations <<
			"   Restarts (line searches):  " << LineSearchCunter << endl;
	}
};

template<typename line_search>
void cm::MHB::operator()(line_search& lnSrch, const double BETA) noexcept
{
	++ValGradEvaluations;
	f0 = f1 = ppr->valGrad(x0, g0);

	if (ppr->stoppingCondition(g0))
		return;

	while (true)
	{
		if (f0 <= f1 || counter == 100000)
		{
			dir = g0;
			++LineSearchCunter;
			lnSrch();
			counter = 0;
		}
		
		for (int i = 0; i < n; ++i)
			moment[i] = x1[i] - x0[i];
		swap(x0, x1);
		swap(g0, g1);

		if (ppr->stoppingCondition(g0))
			return;

		// uncomment for "runANNTest()"
		// cout << "f0=" << f0<<   "   " << "|g| =" <<  infNorm(g0,n)  << ",  evaluations:  " << ValGradEvaluations << endl;   //

		f0 = f1;
		for (int i = 0; i < n; ++i)
			x1[i] = x0[i] + BETA* moment[i] - alpha * g0[i];

		f1 = ppr->valGrad(x1, g1);
		++ValGradEvaluations;
		++counter;
	}
}