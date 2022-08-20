#pragma once  
#include <string>
using namespace std;

class phasePrimitives
{
public:
	long int n{};
	long int FuncGradEvaluations{};
	double f0{};		//objective function's value at x0 - მიზნის ფუნქციის მნიშვნელობა x0-ში
	double f1{};
	double* x0, * x1;   //vectors to iterate
	double* moment;		//to be used in 1d minimization
	double* g1, * g0;   //to store gradients, corresponding x0 and x1
	double* dir;		//descent direction. For MHB simply anti-gradient

//	phasePrimitives(problem* problemAddress) :n{ problemAddress->getSize() }
	phasePrimitives(long int k) :n{ k }
	{
		moment = new double[n];
		x0 = new double[n];
		x1 = new double[n];
		g1 = new double[n];
		g0 = new double[n];
		dir = new double[n];
	}

	~phasePrimitives()
	{
		delete[] moment;
		delete[] x0;
		delete[] x1;
		delete[] g0;
		delete[] g1;
		delete[] dir;
	}
};