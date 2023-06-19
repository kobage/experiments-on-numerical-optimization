#pragma once			
#include<iostream>
#include "auxiliaryFunctions.h" 
#include"Problems.h"
using namespace std;

namespace cm		//collection of minimizers
{
	int n;
	double EPS;
	int mVALUE;				//a default value for parameter that determines the number of BFGS corrections saved, usually 11
	double f0(0.0);			//objective function's value at x0 
	double f1(0.0);
	double alpha{ 0.001 };
	double* x0, * dir, * tmp;
	double* g1, * g0;
	double* x1;

	long int ValGradEvaluations{};
	long int ValueEvaluations{};
	long int LineSearchCunter{};

	problem* ppr{};   //pointer to a problem

	//line search
	class lineSearch;

	//Solvers
	class L_BFGS;
 	class MHB;
	class MNAG;
	class ADAM;
};

