#include<iostream>
#include<vector>
#include<ranges>
#include<numeric>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include"problems.h"
#include"solver.h"

using boost::multiprecision::cpp_dec_float_50;
using namespace std;


int main()
{
	using
		real = cpp_dec_float_50;
	//	real = double;

//	problem<real>* a = new CURLY10<real>();
	problem<real>* a = new Woods<real>();
	
	MHB<real> minimizer{ 1e-6,a};

	freopen("results.txt", "w", stdout);
	auto st = chrono::high_resolution_clock::now();

	a->initialize(minimizer.getVariables());
	bool what = minimizer.solver();

	auto diff = chrono::high_resolution_clock::now() - st;
	auto time = chrono::duration_cast<chrono::microseconds>(diff);

	std::cout << "Optimal Value:  "<< minimizer.getValue() << endl;
	std::cout << "Time:  "  << time.count() << endl << endl;
	if(what)
		minimizer.printVector(minimizer.getVariables(),"x");

	cout << what << endl;
}