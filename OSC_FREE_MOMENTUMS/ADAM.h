#pragma once  
#include <string>
#include"problems.h"
#include<algorithm>
#include "phasePrimitives.h"

using namespace std;

class ADAM
{
	long int& n;
	double& f0;
	double& f1;
	double*& x0;
	double*& x1;
	double*& g0;
	double*& g1;
	long int& FuncGradEvaluations;

	long int RestartCunter{};
	double const ALPHA{ 0.001 };
	double const BETA1{ 0.9 };
	double const BETA2{ 0.999 };
	double const EPSILON { 1e-8 };

	double t{};
	double* m0{};
	double* v0{};
	double* m_tilda;
	double* v_tilda;
public:
	ADAM(phasePrimitives& prmt)
		: n{ prmt.n }, f0{ prmt.f0 }, f1{ prmt.f1 }, x0{ prmt.x0 }, x1{ prmt.x1 },
		g0{ prmt.g0 }, g1{ prmt.g1 },  FuncGradEvaluations{ prmt.FuncGradEvaluations }	
	{
		m0 = new double[n];
		v0 = new double[n];
		m_tilda = new double[n];
		v_tilda = new double[n];
	
	}
	~ADAM() 
	{ 
		delete[] m0;
		delete[] v0;
		delete[] m_tilda;
		delete[] v_tilda;
	}

	long int getSize() const noexcept { return n; }
	long int getRestarts() const noexcept { return RestartCunter; }
	long int getEvaluations() const noexcept { return FuncGradEvaluations; }
	double* getVariables() const noexcept { return x0; }
	double* getGradient() const noexcept { return g0; }
	double getValue() const noexcept { return f0; }

	//solvers
	template<class ProblemPtr>
	void solve(ProblemPtr ppr) noexcept;

	void printStats()
	{
		std::cout << "Algorithm: ADAM,    Objective function:  " << f0 <<
			",   Value-Grad evaluations:  " << getEvaluations() <<
			",   Restarts (line searches):  " << getRestarts() << std::endl;
	}
};

template<class ProblemPtr>
void  ADAM::solve(ProblemPtr ppr) noexcept
{
	FuncGradEvaluations = 0;
	t = 0;
	for (int i = 0; i < n; ++i)
	{
		m0[i] = 0;
		v0[i] = 0;
	}

	while (true)
	{
		++t;
		f0 = ppr->valGrad(x0, g0);
		//    cout << "f0=" << f0 << "   "  << infNorm(g0, n) << ",  evaluations:  " << FuncGradEvaluations << endl;///////////////////
		++FuncGradEvaluations;
		if (ppr->stoppingCondition(g0))
			return;		
		for (int i = 0; i < n; ++i)
		{
			m0[i] = BETA1 * m0[i] + (1 - BETA1) * g0[i];
			v0[i] = BETA2 * v0[i] + (1 - BETA2) * (g0[i]* g0[i]);
			m_tilda[i] = m0 [i] / (1 - pow(BETA1, t));
			v_tilda[i] = v0 [i] / (1 - pow(BETA2, t));
			x0[i] = x0[i] - (ALPHA * m_tilda[i]) / (sqrt(v_tilda[i]) + EPSILON);
		}
	}
}
