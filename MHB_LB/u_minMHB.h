#pragma once  
#include <string>
using namespace std;

namespace u_min
{
	int n;
	long int counter;
	double EPS;
	double alpha;
	double f0;			//objective function's value at x0 - მიზნის ფუნქციის მნიშვნელობა x0-ში
	double f1;
	double  *x0, *dir, *tmp;
	double * g1, *g0;
	double * x1;

	//functions
	void constructData(double tolerance);
	void destructData();

	template<typename Functor, typename Problem, typename BreakType>
	void mhb(Functor, Problem, BreakType) noexcept;

	template<typename Functor, typename Problem, typename BreakType>
	void steepestDescent(Functor, Problem, BreakType) noexcept;
};

void u_min::constructData(double tolerance)
{
	EPS = tolerance;
	alpha = 0.001;
	dir = (double *)malloc(n * sizeof(double));
	x0 = (double *)malloc(n * sizeof(double));
	x1 = (double *)malloc(n * sizeof(double));
	g1 = (double *)malloc(n * sizeof(double));
	g0 = (double *)malloc(n * sizeof(double));
}

void u_min::destructData()
{
	free(dir);
	free(x0);
	free(x1);
	free(g0);
	free(g1);
}

template<typename Functor, typename Problem, typename BreakType>
void u_min::mhb(Functor lnSrchMethod, Problem functionName, BreakType breaking) noexcept
{
	f0 = f1 = functionName(x0, g0, n);
	if (breaking()) 	return;
	while (true)
	{
		if (f0 <= f1 || counter == 100000)
		{
			lnSrchMethod(functionName);
			counter = 0;
		}
		for (int i = 0; i < n; ++i)
			dir[i] = x1[i] - x0[i];
		tmp = x0; x0 = x1; x1 = tmp;
		tmp = g0; g0 = g1; g1 = tmp;

		if (breaking()) 	return;

		f0 = f1;
		for (int i = 0; i < n; ++i)
			x1[i] = x0[i] + dir[i] - alpha*g0[i];
		f1 = functionName(x1, g1, n);
		++counter;
	}
}

template<typename Functor, typename Problem, typename BreakType>
void u_min::steepestDescent(Functor lnSrchMethod, Problem functionName, BreakType breaking) noexcept
{
	f0 = f1 = functionName(x0, g0, n);
	if (breaking()) 	return;
	while (true)
	{

		lnSrchMethod(functionName);
		tmp = x0; x0 = x1; x1 = tmp;
		tmp = g0; g0 = g1; g1 = tmp;
		f0 = f1;
		if (breaking()) 	return;
	}
}
