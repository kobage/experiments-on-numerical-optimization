#pragma once
#include<iostream>
using namespace std;

namespace hb
{
	int n;
	double EPS;
	double alpha(0.001);
	double f0(0.0);  //es x0-shi mnisvneloba			
	double f1(0.0);

	const double delta(0.1);	// used in the Wolfe conditions
	const double sigma(0.9);	// used in the Wolfe conditions
	const double epsilon(1E-6); // used in the approximate Wolfe termination T2
	const double gama(0.66);	// determines when a bisection step s performed in (L2)
	const double ro(5.0);		//expansion factor used in bracket rule
	const double fsi0(0.01);	//small factor used in the starting guess I0
	const double fsi2(2.0);		// factor to multiply pevious step's alpha in I2
	const double teta(0.5);

	double  a, b, c, c1, d;
	double 	fi0, fiA, fiB, fiC, fiD;
	double  fiPrime0, fiPrime_a, fiPrime_b, fiPrime_c, fiPrime_d;
	double gNorm;

	double(*functionName)
		(
			double *x,
			double *g,
			const int n
			);
	double  *x0, *dir, *tmp;
	double  *x2;
	double * g1, *g0;
	double * x1;
	//functions
	void setDimentions(int n, double tolerance);
	void lineSearch();
	void mhb(void lnSrchMethod(void));
};

inline void vec2norm(double *ans, double *gr, const int n)
{
	int i, n5;
	*ans = 0.;
	n5 = n % 5;
	for (i = 0; i < n5; i++) *ans += gr[i] * gr[i];
	for (; i < n; i += 5)
	{
		*ans += gr[i] * gr[i] + gr[i + 1] * gr[i + 1] + gr[i + 2] * gr[i + 2]
			+ gr[i + 3] * gr[i + 3] + gr[i + 4] * gr[i + 4];
	}
	*ans = sqrt(*ans);
}

inline double vecProd(double *x, double *y, const int n)
{
	int i, n5;
	double sum(0.0);
	if (n <= 0) return sum;
	n5 = n % 5;
	for (i = 0; i < n5; i++) sum += x[i] * y[i];
	for (; i < n; i += 5)
	{
		sum += x[i] * y[i] + x[i + 1] * y[i + 1] + x[i + 2] * y[i + 2]
			+ x[i + 3] * y[i + 3] + x[i + 4] * y[i + 4];
	}
	return sum;
}

inline double infNorm
(
	double *x, /* vector */
	int     n /* length of vector */
)
{
	int i, n5;
	double t = 0.;
	n5 = n % 5;

	for (i = 0; i < n5; i++) if (t < fabs(x[i])) t = fabs(x[i]);
	for (; i < n; i += 5)
	{
		if (t < fabs(x[i])) t = fabs(x[i]);
		if (t < fabs(x[i + 1])) t = fabs(x[i + 1]);
		if (t < fabs(x[i + 2])) t = fabs(x[i + 2]);
		if (t < fabs(x[i + 3])) t = fabs(x[i + 3]);
		if (t < fabs(x[i + 4])) t = fabs(x[i + 4]);
	}
	return (t);
}

void hb::setDimentions(int k, double tolerance)
{
	n = k;
	EPS = tolerance;
	dir = (double *)malloc(k * sizeof(double));
	x0 = (double *)malloc(k * sizeof(double));
	x1 = (double *)malloc(k * sizeof(double));
	x2 = (double *)malloc(k * sizeof(double));
	g1 = (double *)malloc(k * sizeof(double));
	g0 = (double *)malloc(k * sizeof(double));
}
//--------Line Search-------------
void hb::lineSearch()
{

	fi0 = f0;
	fiPrime0 = -vecProd(g0, g0, n);
	gNorm = 0.;

	double length, t;
	//I:< Initial gess for b >
	//I.begin() 

	vec2norm(&gNorm, g0, n);
	t = infNorm(x0, n);
	if (t != 0.)
		b = (fsi0*t) / infNorm(g0, n);
	else
		if (fi0 != 0.)
			b = (fsi0*fabs(fi0)) / (gNorm*gNorm);
		else
			b = 1.;

	//I.end()

	// < Generating initial [a,b]>;
	a = 0.;
	fiPrime_a = fiPrime0;
	fiA = fi0;

	while (1)
	{
		//< Termination test in b >;
		//T.begin()
		for (int i = 0; i < n; i++)
			x1[i] = x0[i] - b*g0[i];
		fiB = functionName(x1, g1, n);
		fiPrime_b = -vecProd(g1, g0, n);
		if (fiPrime_b <= 0.   &&   fiPrime_b >= sigma * fiPrime0 	&&   fiB <= fi0 + epsilon)
		{
			alpha = b;
			f1 = fiB;
			return;
		}
		//T.end()
		if (fiPrime_b >= 0)  break;
		if (fiB > fi0 + epsilon)
		{
			// < Bisection on [a,b] >
			//B.begin()
			while (1)
			{
				d = (1 - teta)*a + teta*b;
				//< Termination test in d >;
				//T.begin()
				for (int i = 0; i < n; i++)
					x1[i] = x0[i] - d*g0[i];
				fiD = functionName(x1, g1, n);
				fiPrime_d = -vecProd(g1, g0, n);
				if (fiPrime_d <= 0.   &&   fiPrime_d >= sigma * fiPrime0   &&   fiD <= fi0 + epsilon)
				{
					alpha = d;
					f1 = fiD;
					return;
				}
				//T.end()
				if (fiPrime_d >= 0.)
				{
					b = d;
					break;
				}
				if (fiD < fi0 + epsilon)
				{

					a = d;
					fiA = fiD;
					fiPrime_a = fiPrime_d;
				}
				else
					b = d;
			}
			//B.end()
			break;
		}
		a = b;
		fiA = fiB;
		fiPrime_a = fiPrime_b;
		b *= ro;
	}
	//G.end()
	while (1)
	{
		length = b - a;
		c = (a*fiPrime_b - b*fiPrime_a) / (fiPrime_b - fiPrime_a);
		// < Update  a,b,c >
		// U3.begin()
		if (a < c && c < b)
		{
			// < Termination test in c >
			// T.begin()
			for (int i = 0; i < n; i++)
				x1[i] = x0[i] - c*g0[i];
			fiC = functionName(x1, g1, n);
			fiPrime_c = -vecProd(g1, g0, n);
			if (fiPrime_c <= 0.   &&   fiPrime_c >= sigma * fiPrime0    &&    fiC <= fi0 + epsilon)
			{
				alpha = c;
				f1 = fiC;
				return;
			}
			// T.end()
			if (fiPrime_c < 0.   &&  fiC > fi0 + epsilon)
			{
				b = c;
				// < Bisection on [a,b] >
				//B.begin()
				while (1)
				{
					d = (1 - teta)*a + teta*b;
					//< Termination test in d >;
					//T.begin()
					for (int i = 0; i < n; i++)
						x1[i] = x0[i] - d*g0[i];
					fiD = functionName(x1, g1, n);
					fiPrime_d = -vecProd(g1, g0, n);
					if (fiPrime_d <= 0.   &&   fiPrime_d >= sigma * fiPrime0   &&   fiD <= fi0 + epsilon)
					{
						alpha = d;
						f1 = fiD;
						return;
					}
					//T.end()

					if (fiPrime_d >= 0.)
					{
						b = d;
						break;
					}
					if (fiD < fi0 + epsilon)
					{
						a = d;
						fiA = fiD;
						fiPrime_a = fiPrime_d;
					}
					else
						b = d;
				}
				//B.end()
			}
			else
				if (fiPrime_c >= 0.)
				{
					c1 = (c*fiPrime_b - b*fiPrime_c) / (fiPrime_b - fiPrime_c);
					b = c;
				}
				else
				{
					c1 = (a*fiPrime_c - c*fiPrime_a) / (fiPrime_c - fiPrime_a);
					a = c;
					fiA = fiC;
					fiPrime_a = fiPrime_c;
				}
			c = c1;
			//U3.end()
		}
		// < Update a,b>
		// U2.begin()
		if (a < c && c < b)
		{

			// < Termination test in c >
			// T.begin()
			for (int i = 0; i < n; i++)
				x1[i] = x0[i] - c*g0[i];
			fiC = functionName(x1, g1, n);
			fiPrime_c = -vecProd(g1, g0, n);
			if (fiPrime_c <= 0.   &&   fiPrime_c >= sigma * fiPrime0   &&   fiC <= fi0 + epsilon)
			{
				alpha = c;
				f1 = fiC;
				return;
			}
			// T.end()
			if (fiPrime_c < 0.   &&  fiC > fi0 + epsilon)
			{
				b = c;
				// < Bisection on [a,b] >
				//B.begin()
				while (1)
				{
					d = (1 - teta)*a + teta*b;
					//< Termination test in d >;
					//T.begin()
					for (int i = 0; i < n; i++)
						x1[i] = x0[i] - d*g0[i];
					fiD = functionName(x1, g1, n);
					fiPrime_d = -vecProd(g1, g0, n);
					if (fiPrime_d <= 0.   &&   fiPrime_d >= sigma * fiPrime0   &&  fiD <= fi0 + epsilon)
					{
						alpha = d;
						f1 = fiD;
						return;
					}
					//T.end()
					if (fiPrime_d >= 0.)
					{
						b = d;
						break;
					}
					if (fiD < fi0 + epsilon)
					{
						a = d;
						fiA = fiD;
						fiPrime_a = fiPrime_d;
					}
					else
						b = d;
				}
				//B.end()
			}
			else
				if (fiPrime_c >= 0.)
					b = c;
				else
				{
					a = c;
					fiA = fiC;
					fiPrime_a = fiPrime_c;
				}
			//U2.end()
		}
		if (b - a > length * gama)
		{
			c = (a + b) / 2;
			// < Update >
			// U.begin()
			// < Termination test in c >
			// T.begin()
			for (int i = 0; i < n; i++)
				x1[i] = x0[i] - c*g0[i];
			fiC = functionName(x1, g1, n);
			fiPrime_c = -vecProd(g1, g0, n);
			if (fiPrime_c <= 0.   &&   fiPrime_c >= sigma * fiPrime0 	&&   fiC <= fi0 + epsilon)
			{
				alpha = c;
				f1 = fiC;
				return;
			}
			// T.end()
			if (fiPrime_c < 0.   &&  fiC > fi0 + epsilon)
			{
				b = c;	fiB = fiC;	fiPrime_b = fiPrime_c;
				// < Bisection on [a,b] >
				//B.begin()
				while (1)
				{
					d = (1 - teta)*a + teta*b;
					//< Termination test in d >;
					//T.begin()
					for (int i = 0; i < n; i++)
						x1[i] = x0[i] - d*g0[i];
					fiD = functionName(x1, g1, n);
					fiPrime_d = -vecProd(g1, g0, n);
					if (fiPrime_d <= 0.   &&   fiPrime_d >= sigma * fiPrime0   &&   fiD <= fi0 + epsilon)
					{
						alpha = d;
						f1 = fiD;
						return;
					}
					//T.end()
					if (fiPrime_d >= 0.)
					{
						b = d;
						break;
					}
					if (fiD < fi0 + epsilon)
					{
						a = d;
						fiA = fiD;
						fiPrime_a = fiPrime_d;
					}
					else
						b = d;
				}
				//B.end()
			}
			else
				if (fiPrime_c >= 0.)
					b = c;
				else
				{
					a = c;
					fiA = fiC;
					fiPrime_a = fiPrime_c;
				}
			//U.end()
		}
	}
}
//end of line search 

void hb::mhb(void lnSrchMethod(void))
{
	f0 = functionName(x0, g0, n);
	f1 = f0;
	if (infNorm(g0, n) < EPS) 	return;
	long int counter;
	while (true)
	{
		if (f0 <= f1 || counter == 100000)
		{
			lnSrchMethod();
			counter = 0;
		}
		for (int i = 0; i < n; ++i)
			dir[i] = x1[i] - x0[i];
		tmp = x0; x0 = x1; x1 = tmp;
		tmp = g0; g0 = g1; g1 = tmp;

		if (infNorm(g0, n) < EPS) 	return;

		f0 = f1;
		for (int i = 0; i < n; ++i)
			x1[i] = x0[i] + dir[i] - alpha*g0[i];
		f1 = functionName(x1, g1, n);
		++counter;
	}
}