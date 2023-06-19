#pragma once  
#include <string>   
#include<algorithm>
#include"cm_namespace.h"

using namespace std;

class cm::ADAM
{
	double const BETA1{ 0.9 };
	double const BETA2{ 0.999 };
	double const EPSILON{ 1e-8 };

	double t{};
	double* m0{};
	double* v0{};
	double* m_tilda;
	double* v_tilda;
public:
	ADAM(problem* p)
	{
		n = p->getSize();
		ppr = p;
		m0 = new double[n];
		v0 = new double[n];
		m_tilda = new double[n];
		v_tilda = new double[n];
		x0 = new double[n];
		g0 = new double[n];
	}
	~ADAM()
	{
		delete[] m0;
		delete[] v0;
		delete[] m_tilda;
		delete[] v_tilda;
		delete[]  x0;
		delete[]  g0;
	}

	template<typename line_search>
	void operator()(line_search&) noexcept;

	void printStats()
	{
		std::cout << "   Algorithm ADAM;  Objective function: " << f0 << ",  Value-Grad evaluations:  " << ValGradEvaluations << std::endl;
		std::cout << "   Restarts (line searches):  " << LineSearchCunter;
	}
};

template<typename line_search>
void cm::ADAM::operator()(line_search& lnSrch) noexcept
{
	for (int i = 0; i < n; ++i)
	{
		m0[i] = 0;
		v0[i] = 0;
	}

	while (true)
	{
		++t;
		f0 = ppr->valGrad(x0, g0);
	//	    cout << "f0=" << f0 << "   "  << infNorm(g0, n) << ",  evaluations:  " << ValGradEvaluations << endl;///////////////////
		++ValGradEvaluations;
		if (ppr->stoppingCondition(g0))
			return;
		for (int i = 0; i < n; ++i)
		{
			m0[i] = BETA1 * m0[i] + (1 - BETA1) * g0[i];
			v0[i] = BETA2 * v0[i] + (1 - BETA2) * (g0[i] * g0[i]);
			m_tilda[i] = m0[i] / (1 - pow(BETA1, t));
			v_tilda[i] = v0[i] / (1 - pow(BETA2, t));
			x0[i] = x0[i] - (alpha * m_tilda[i]) / (sqrt(v_tilda[i]) + EPSILON);
		}
	}
}
