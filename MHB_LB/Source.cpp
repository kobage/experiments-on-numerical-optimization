#include"auxiliaryFunctions.h"
#include"u_minMHB.h"
#include"lineSearch.h"
#include"TestsUCON.h"
#include"add-inLB.h"

using namespace u_min;


int main()
{
	//To run unconstrained minimization test problems;
	//runUCONTests();

	//To run lower-bound constrained minimization test problems
	//solver can be changed on lines a442-1443 and 1445-1446, in file "add-inLB"
	runTests_LB();
}