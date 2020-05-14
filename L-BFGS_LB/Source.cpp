#include"auxiliaryFunctions.h"
#include"u_minLBFGS.h"
#include"lineSearch.h"
#include"TestsUCON.h"
#include"add-inLB.h"

using namespace u_min;


int main()
{
	//To run unconstrained minimization test problems
	//runUCONTests();

	//To run lower-bound constrained minimization test problems
	runTests_LB();
}