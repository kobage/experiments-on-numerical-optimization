l-bfgs is implemented by us and its design fully meets the requirements that arise when 
integrating the add-in into unconstrained mnimization algorithm. This implementation is 
enough simple. The only modern technic is line search procedure borrowed from CG-DESCENT-C-6.8. 
But its design is object oriented and partly generic. As a result, implementation is effective 
both on unconstrained minimization problems and on problems, 
where variables are constrained with only their lower bounds.  

To evaluate performance characteristic of the presented add-in, integrated into l-bfgs, 
on symmetric matrix game tests, 
the following code can be used. It uses the fastest mathematical 
programming solver  Gurobi optimizer 800 with C++ interface. 
To run with Visual Studio 2017, the following files should be added (Project -> Add Existing Item):
gurobi80.lib, gurobi_c++md2017.lib, gurobi_c++.h.

Varying commented lines, different ways of collecting data are tried.

#include <iostream>
#include <chrono>
#include <random>
#include "C:\gurobi800\win64\include\gurobi_c++.h"
using namespace std;

static bool
dense_optimize(GRBEnv* env,
	int     rows,
	int     cols,
	double** A,     /* constraint matrix */
	char*   sense,  /* constraint senses */
	double* rhs,    /* RHS vector */
	double* lb,     /* variable lower bounds */
	double* ub,     /* variable upper bounds */
	char*   vtype,  /* variable types (continuous, binary, etc.) */
	double* solution,
	double* objvalP)
{
	GRBModel model = GRBModel(*env);
	model.set(GRB_IntParam_OutputFlag, false);

	int i, j;
	bool success = false;

	/* Add variables to the model */
	GRBVar* vars = model.addVars(lb, ub, NULL, vtype, NULL, cols);

	for (i = 0; i < rows; i++) {
		GRBLinExpr lhs = 0;
		for (j = 0; j < cols; j++)
			if (A[i][j] != 0)
				lhs += A[i][j] * vars[j];
		model.addConstr(lhs, sense[i], rhs[i]);
	}
	{
		GRBLinExpr lhs = 0;
		for (i = 0; i < rows; i++)
			lhs += vars[i];
		model.addConstr(lhs, '=', 1.);
	}

	GRBQuadExpr obj = 0;
	model.setObjective(obj);

	model.optimize();

	//	model.write("dense.lp");

	if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
		*objvalP = model.get(GRB_DoubleAttr_ObjVal);
		for (i = 0; i < cols; i++)
			solution[i] = vars[i].get(GRB_DoubleAttr_X);
		success = true;
	}
	delete[] vars;
	return success;
}

int  main()
{
	GRBEnv* env = 0;

	int n;
	cout << "Enter dimension: " << endl;
	cin >> n;

	try {
		auto st = chrono::high_resolution_clock::now();

		env = new GRBEnv();

		double **A;
		char*    sense;
		double*  rhs;
		double*  lb;
		bool    success;
		double  objval;
		double*  sol;

		//fillBasicData(n, A, sense, rhs, lb, sol);
		A = (double **)malloc(sizeof(double *)*n);

		for (int i = 0; i < n; i++)
			A[i] = (double *)malloc(sizeof(double)*n);

		std::default_random_engine dre;
		std::uniform_real_distribution<double> urdi(-1, 1);
		//std::normal_distribution<double> urdi(-1,1);
		//std::cauchy_distribution<double> urdi(-1, 1);

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < i; j++)
				A[i][j] = urdi(dre);
			A[i][i] = 0;
		}
		for (int i = 0; i < n; i++)
			for (int j = i + 1; j < n; j++)
				A[i][j] = (-1.0* A[j][i]);

		lb = (double *)malloc(sizeof(double)*n);
		for (int i = 0; i < n; ++i)
			lb[i] = 0.;

		rhs = (double *)malloc(sizeof(double)*n);
		for (int i = 0; i < n; ++i)
			rhs[i] = 0.;

		sol = (double *)malloc(sizeof(double)*n);

		sense = (char *)malloc(sizeof(char)*n);
		for (int i = 0; i < n; ++i)
			sense[i] = '<';

		freopen("results.txt", "w", stdout);

		success = dense_optimize(env, n, n, A, sense, rhs, lb, NULL, NULL, sol, &objval);

		auto diff = chrono::high_resolution_clock::now() - st;
		auto time = chrono::duration_cast<chrono::microseconds>(diff);
		cout << "Time: " << time.count() << endl;
		cout << "The minimum of Onjective Function: " << objval << endl;
		int quarter = n / 4;
		for (int i = 0; i < quarter; i++)
		{
			printf("x[%3d] = %-15.11f     x[%3d] = %15.11f     x[%3d] = %15.11f     x[%3d] = %2.11f  \n",
				i, sol[i], quarter + i, sol[quarter + i], 2 * quarter + i, sol[2 * quarter + i], 3 * quarter + i, sol[3 * quarter + i]);
		}
	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Exception during optimization" << endl;
	}

	delete env;
}