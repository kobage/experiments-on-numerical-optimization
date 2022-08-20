#include<vector>
#include<chrono>
#include<algorithm>
#include "L_BFGS.h"
#include "lineSearch.h"
#include "testers.h"
using namespace std;

int main()
{
	//To run unconstrained minimization test problems
//	runUCONTests();
// 
	//To run ANN test problems
	runANNTest();
}