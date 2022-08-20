#pragma once
#include<numeric>
#include<vector>
using namespace std;
class quadrPenalSMG 
{
	double Dev{};
	int i{}, j{};

public:
	vector<vector<double>> A;
	double EPS{};
	long int n{};
	quadrPenalSMG() {}
	quadrPenalSMG(double const tolerance) { EPS = tolerance; n = 200; }
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "quadrPenalSMG"; }
	double valGrad
	(
		double* x,
		double* g
	)
	{
		double fx = 0.;
		Dev = 0.;
		double prod{};
		double sum{};
		double* p;
		for (i = 0; i < n; i++)
		{
			g[i] = 0;
			sum += x[i];
		}
		for (i = 0; i < n; i++)
		{
			p = A[i].data();
			prod = std::inner_product(p, p + n, x, double{});
			if (prod > 0.)
			{
				if (Dev < prod) Dev = prod;
				fx += prod * prod;
				for (j = 0; j < n; j++)
					g[j] += prod * p[j];
			}
		}
		if (Dev < fabs(sum - 1)) Dev = fabs(sum - 1);
		fx += (sum - 1) * (sum - 1);
		for (i = 0; i < n; i++)
			g[i] += (sum - 1);

		for (i = 0; i < n; i++)
		{
			if (x[i] < 0.)
			{
				if (-x[i] > Dev)
					Dev = -x[i];
				prod = x[i];
				fx += prod * prod;
				g[i] += prod;
			}
		}

		return fx / 2;
	}

	double value
	(
		double* x
	)
	{
		double fx = 0.;
		Dev = 0.;
		double prod{};
		double sum{};
		double* p;
		for (i = 0; i < n; i++)
			sum += x[i];
		
		for (i = 0; i < n; i++)
		{
			p = A[i].data();
			prod = std::inner_product(p, p + n, x, double{});
			if (prod > 0.)
			{
				if (Dev < prod) Dev = prod;
				fx += prod * prod;
			}
		}
		if (Dev < fabs(sum - 1)) Dev = fabs(sum - 1);
		fx += (sum - 1) * (sum - 1);


		for (i = 0; i < n; i++)
		{
			if (x[i] < 0.)
			{
				if (-x[i] > Dev)
					Dev = -x[i];
				prod = x[i];
				fx += prod * prod;
			}
		}
		return fx / 2;
	}

	void Grad
	(
		double* x,
		double* g
	)
	{
		Dev = 0.;
		double prod{};
		double sum{};
		double* p;
		for (i = 0; i < n; i++)
		{
			g[i] = 0;
			sum += x[i];
		}
		for (i = 0; i < n; i++)
		{
			p = A[i].data();
			prod = std::inner_product(p, p + n, x, double{});
			if (prod > 0.)
			{
				if (Dev < prod) Dev = prod;
				for (j = 0; j < n; j++)
					g[j] += prod * p[j];
			}
		}
		if (Dev < fabs(sum - 1)) Dev = fabs(sum - 1);
		for (i = 0; i < n; i++)
			g[i] += (sum - 1);

		for (i = 0; i < n; i++)
		{
			if (x[i] < 0.)
			{
				if (-x[i] > Dev)
					Dev = -x[i];
				prod = x[i];
				g[i] += prod;
			}
		}
	}
	void printVector()
	{
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
				cout << A[i][j] << "  ";
			cout << endl;
		}
		cout << endl << endl;
	}

	void initialize(double* x)
	{
		for (i = 0; i < n; ++i)
		{
			x[i] = 1. / n;
		}
		A.resize(n);
		for (i = 0; i < n; ++i)
		{
			for (j = 0; j < i; ++j)
				A[i].push_back((2.0 * rand()) / RAND_MAX - 1.0);
			A[i].push_back(0.);
		}
		for (i = 0; i < n; i++)
			for (j = i + 1; j < n; j++)
				A[i].push_back(-1.0 * A[j][i]);
	}
	bool  stoppingCondition() const
	{
		return (Dev < EPS);
	}
};

quadrPenalSMG game(1e-6);

double MGValue
(
	double* x,
	INT       n
)
{
	return game.value(x);
}

void MGGrad
(
	double* g,
	double* x,
	INT        n
)
{
	game.Grad(x,g);
}

double MGValGrad
(
	double* g,
	double* x,
	INT        n
)
{
	return game.valGrad(x, g);
}
#include <iostream>
using namespace std;

bool stoppingByDeviation()
{
	return game.stoppingCondition();
}

void runSymGameTest(void)
{
	cg_stats Stats;
	long int n = game.n;
	double* x;
	x = (double*)malloc(n * sizeof(double));
	game.initialize(x);

	cg_parameter Parm;
	cg_default(&Parm);
	Parm.stop = stoppingByDeviation;

	std::chrono::high_resolution_clock::time_point  st;  
	std::chrono::high_resolution_clock::duration  diff;  

	st = chrono::high_resolution_clock::now();
	cg_descent(x, n, &Stats, &Parm, 123123, MGValue, MGGrad, MGValGrad, NULL);/////////////////
	diff = chrono::high_resolution_clock::now() - st;

	freopen("results.txt", "w", stdout);
	printf("Median time in microseconds:   %d\n",std::chrono::duration_cast<std::chrono::microseconds>(diff).count());
	{
		printf("Function's value:   %f\n", Stats.f);
		int quarter = n / 4;
		int ii;
		for (ii = 0; ii < quarter; ii++)
		{
			printf("x[%3d] = %-15.11f     x[%3d] = %15.11f     x[%3d] = %15.11f     x[%3d] = %2.11f  \n",
				ii, x[ii], quarter + ii, x[quarter + ii], 2 * quarter + ii, x[2 * quarter + ii], 3 * quarter + ii, x[3 * quarter + ii]);
		}
	}
	free(x);
}