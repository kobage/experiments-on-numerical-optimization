#pragma once
#include "u_minMHB.h"
using namespace u_min;

struct lineSearch
{
	const double DELTA = 0.1;	// used in the Wolfe conditions
	const double SIGMA = 0.9;	// used in the Wolfe conditions

	const double EPSILON = 1E-6; // used in the approximate Wolfe termination T2
	const double GAMA = 0.66;	// determines when a bisection step s performed in (L2)
	const double RO = 5.0;		//expansion factor used in bracket rule
	const double FSI0 = 0.01;	//small factor used in the starting guess I0
	const double FSI2 = 2.0;		// factor to multiply previous step's alpha in I2
	const double TETA = 0.5;

	double  a, b, c, c1, d;
	double 	fiB, fiC, fiD;
	double  fiPrime0, fiPrime_a, fiPrime_b, fiPrime_c, fiPrime_d;
	lineSearch() {}
	template<typename Problem>
	void operator()(Problem)  noexcept;
};
template<typename Problem>
void lineSearch::operator()(Problem functionName) noexcept
{
	fiPrime0 = -vecProd(g0, g0, n);
	double length, t, fiPrime_Low;;
	fiPrime_Low = SIGMA * fiPrime0;

	//I:< Initial gess for b >
	//I.begin() 
	t = infNorm(x0, n);
	if (t != 0.)
		b = (FSI0 * t) / infNorm(g0, n);
	else
		if (f0 != 0.)
			b = -(FSI0 * fabs(f0)) / fiPrime0;
		else
			b = 1.;
	//I.end()

	// < Generating initial [a,b]>;
	a = 0.;
	fiPrime_a = fiPrime0;

	while (1)
	{
		//< Termination test in b >;
		//T.begin()
		for (int i = 0; i < n; i++)
			x1[i] = x0[i] - b * g0[i];
		fiB = functionName(x1, g1, n);
		fiPrime_b = -vecProd(g1, g0, n);
		if (fiPrime_b <= 0. && fiPrime_b >= fiPrime_Low && fiB <= f0 + EPSILON)
		{
			alpha = b;
			f1 = fiB;
			return;
		}
		//T.end()
		if (fiPrime_b >= 0)  break;
		if (fiB > f0 + EPSILON)
		{
			// < Bisection on [a,b] >
			//B.begin()
			while (1)
			{
				d = (1 - TETA) * a + TETA * b;
				//< Termination test in d >;
				//T.begin()
				for (int i = 0; i < n; i++)
					x1[i] = x0[i] - d * g0[i];
				fiD = functionName(x1, g1, n);
				fiPrime_d = -vecProd(g1, g0, n);
				if (fiPrime_d <= 0. && fiPrime_d >= fiPrime_Low && fiD <= f0 + EPSILON)
				{
					alpha = d;
					f1 = fiD;
					return;
				}
				//T.end()
				if (fiPrime_d >= 0.)
				{
					b = d;
					fiPrime_b = fiPrime_d;
					break;
				}
				if (fiD < f0 + EPSILON)
				{

					a = d;
					fiPrime_a = fiPrime_d;
				}
				else
					b = d;
			}
			//B.end()
			break;
		}
		a = b;
		fiPrime_a = fiPrime_b;
		b *= RO;
	}
	//G.end()
	while (1)
	{
		length = b - a;
		c = (a * fiPrime_b - b * fiPrime_a) / (fiPrime_b - fiPrime_a);
		// < Update  a,b,c >
		// U3.begin()
		if (a < c && c < b)
		{
			// < Termination test in c >
			// T.begin()
			for (int i = 0; i < n; i++)
				x1[i] = x0[i] - c * g0[i];
			fiC = functionName(x1, g1, n);
			fiPrime_c = -vecProd(g1, g0, n);
			if (fiPrime_c <= 0. && fiPrime_c >= fiPrime_Low && fiC <= f0 + EPSILON)
			{
				alpha = c;
				f1 = fiC;
				return;
			}
			// T.end()
			if (fiPrime_c < 0. && fiC > f0 + EPSILON)
			{
				b = c;
				// < Bisection on [a,b] >
				//B.begin()
				while (1)
				{
					d = (1 - TETA) * a + TETA * b;
					//< Termination test in d >;
					//T.begin()
					for (int i = 0; i < n; i++)
						x1[i] = x0[i] - d * g0[i];
					fiD = functionName(x1, g1, n);
					fiPrime_d = -vecProd(g1, g0, n);
					if (fiPrime_d <= 0. && fiPrime_d >= fiPrime_Low && fiD <= f0 + EPSILON)
					{
						alpha = d;
						f1 = fiD;
						return;
					}
					//T.end()

					if (fiPrime_d >= 0.)
					{
						b = d;
						fiPrime_b = fiPrime_d;
						break;
					}
					if (fiD < f0 + EPSILON)
					{
						a = d;
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
					c1 = (c * fiPrime_b - b * fiPrime_c) / (fiPrime_b - fiPrime_c);
					fiPrime_b = fiPrime_c;
					b = c;
				}
				else
				{
					c1 = (a * fiPrime_c - c * fiPrime_a) / (fiPrime_c - fiPrime_a);
					a = c;
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
				x1[i] = x0[i] - c * g0[i];
			fiC = functionName(x1, g1, n);
			fiPrime_c = -vecProd(g1, g0, n);
			if (fiPrime_c <= 0. && fiPrime_c >= fiPrime_Low && fiC <= f0 + EPSILON)
			{
				alpha = c;
				f1 = fiC;
				return;
			}
			// T.end()
			if (fiPrime_c < 0. && fiC > f0 + EPSILON)
			{
				b = c;
				// < Bisection on [a,b] >
				//B.begin()
				while (1)
				{
					d = (1 - TETA) * a + TETA * b;
					//< Termination test in d >;
					//T.begin()
					for (int i = 0; i < n; i++)
						x1[i] = x0[i] - d * g0[i];
					fiD = functionName(x1, g1, n);
					fiPrime_d = -vecProd(g1, g0, n);
					if (fiPrime_d <= 0. && fiPrime_d >= fiPrime_Low && fiD <= f0 + EPSILON)
					{
						alpha = d;
						f1 = fiD;
						return;
					}
					//T.end()
					if (fiPrime_d >= 0.)
					{
						b = d;
						fiPrime_b = fiPrime_d;
						break;
					}
					if (fiD < f0 + EPSILON)
					{
						a = d;
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
					b = c;
					fiPrime_b = fiPrime_c;
				}
				else
				{
					a = c;
					fiPrime_a = fiPrime_c;
				}
			//U2.end()
		}
		if (b - a > length * GAMA)
		{
			c = (a + b) / 2;
			// < Update >
			// U.begin()
			// < Termination test in c >
			// T.begin()
			for (int i = 0; i < n; i++)
				x1[i] = x0[i] - c * g0[i];
			fiC = functionName(x1, g1, n);
			fiPrime_c = -vecProd(g1, g0, n);
			if (fiPrime_c <= 0. && fiPrime_c >= fiPrime_Low && fiC <= f0 + EPSILON)
			{
				alpha = c;
				f1 = fiC;
				return;
			}
			// T.end()
			if (fiPrime_c < 0. && fiC > f0 + EPSILON)
			{
				b = c;
				// < Bisection on [a,b] >
				//B.begin()
				while (1)
				{
					d = (1 - TETA) * a + TETA * b;
					//< Termination test in d >;
					//T.begin()
					for (int i = 0; i < n; i++)
						x1[i] = x0[i] - d * g0[i];
					fiD = functionName(x1, g1, n);
					fiPrime_d = -vecProd(g1, g0, n);
					if (fiPrime_d <= 0. && fiPrime_d >= fiPrime_Low && fiD <= f0 + EPSILON)
					{
						alpha = d;
						f1 = fiD;
						return;
					}
					//T.end()
					if (fiPrime_d >= 0.)
					{
						b = d;
						fiPrime_b = fiPrime_d;
						break;
					}
					if (fiD < f0 + EPSILON)
					{
						a = d;
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
					b = c;
					fiPrime_b = fiPrime_c;
				}
				else
				{
					a = c;
					fiPrime_a = fiPrime_c;
				}
			//U.end()
		}
	}
}