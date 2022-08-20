#pragma once
#include "phasePrimitives.h"
#include"problems.h"

class lineSearch
{
public:
	long int& n;
	double& f0;
	double& f1;
	double*& x0;
	double*& x1;
	double*& g1;
	double*& g0;
	double*& dir;
	long int& FuncGradEvaluations;

	int lineSearchcounter = 0;
	const double DELTA = 0.1;	// used in the Wolfe conditions
	const double SIGMA = 0.9;	// used in the Wolfe conditions
	double alpha = 0.001;		//used in //I:< Initial gess for b >, only
	const double EPSILON = 1E-6; // used in the approximate Wolfe termination T2
	const double GAMA = 0.66;	// determines when a bisection step s performed in (L2)
	const double RO = 5.0;		//expansion factor used in bracket rule
	const double FSI0 = 0.01;	//small factor used in the starting guess I0
	const double FSI2 = 2.0;		// factor to multiply previous step's alpha in I2
	const double TETA = 0.5;
	double  a{}, b{}, c{}, c1{}, d{};
	double  fiB{}, fiC{}, fiD{};
	double  fiPrime0{}, fiPrime_a{}, fiPrime_b{}, fiPrime_c{}, fiPrime_d{};

	lineSearch(phasePrimitives& prmt)
		: n{ prmt.n }, f0{ prmt.f0 }, f1{ prmt.f1 }, x0{ prmt.x0 },
		x1{ prmt.x1 }, g0{ prmt.g0 }, dir{ prmt.dir }, g1{ prmt.g1 }, FuncGradEvaluations{ prmt.FuncGradEvaluations }
	{}
	template<typename Problem> 
	void lineSearchGeneral(Problem)  noexcept;
	template<typename Problem>
	void lineSearchConvex(Problem)  noexcept;

	bool isGeneral{ true };

	template<typename Problem>
	void operator()(Problem p)
	{
		if (isGeneral)
			lineSearchGeneral(p);
		else
			lineSearchConvex(p);
	}
};

template<typename Problem>
void lineSearch::lineSearchGeneral(Problem ppr) noexcept
{
	a = 0.;
	fiPrime_a = fiPrime0 = -std::inner_product(g0, g0 + n, dir, double{});
	double primeLow(SIGMA * fiPrime0), primeUp((2 * DELTA - 1) * fiPrime0);

	double length, t;
	//I:< Initial gess for b >
	//I.begin() 
	if (lineSearchcounter > 0)
		b = FSI2 * alpha;
	else
	{
		t = infNorm(x0, n);
		if (t != 0.)
			b = (FSI0 * t) / infNorm(g0, n);
		else
			if (f0 != 0.)
				b = (FSI0 * fabs(f0)) / std::inner_product(g0, g0 + n, g0, double{});
			else
				b = 1.;
	}
	//I.end()

	// < Generating initial [a,b]>;
	while (1)
	{
		//< Termination test in b >;
		//T.begin()
		for (int i = 0; i < n; i++)
			x1[i] = x0[i] - b * dir[i];
		fiB = ppr->valGrad(x1, g1);
		++FuncGradEvaluations;
		fiPrime_b = -std::inner_product(dir, dir + n, g1, double{});

		if (fiPrime_b <= primeUp && fiPrime_b >= primeLow && fiB <= f0 + EPSILON)
		{
			f0 = fiB;
			++lineSearchcounter;
			return;
		}
		//T.end()

		if (fiPrime_b >= 0)
			break;
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
					x1[i] = x0[i] - d * dir[i];
				fiD = ppr->valGrad(x1, g1);
				++FuncGradEvaluations;
				fiPrime_d = -std::inner_product(dir, dir + n, g1, double{});
				if (fiPrime_d <= primeUp && fiPrime_d >= primeLow && fiD <= f0 + EPSILON)
				{
					f0 = fiD;
					++lineSearchcounter;
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
				x1[i] = x0[i] - c * dir[i];
			fiC = ppr->valGrad(x1, g1);
			++FuncGradEvaluations;
			fiPrime_c = -std::inner_product(dir, dir + n, g1, double{});
			if (fiPrime_c <= primeUp && fiPrime_c >= primeLow && fiC <= f0 + EPSILON)
			{
				f0 = fiC;
				++lineSearchcounter;
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
						x1[i] = x0[i] - d * dir[i];
					fiD = ppr->valGrad(x1, g1);
					++FuncGradEvaluations;
					fiPrime_d = -std::inner_product(dir, dir + n, g1, double{});
					if (fiPrime_d <= primeUp && fiPrime_d >= primeLow && fiD <= f0 + EPSILON)
					{
						f0 = fiD;
						++lineSearchcounter;
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
			{
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
			}
			//U3.end()
			// < Update a,b>
				// U2.begin()
			if (a < c && c < b)
			{
				// < Termination test in c >
				// T.begin()
				for (int i = 0; i < n; i++)
					x1[i] = x0[i] - c * dir[i];
				fiC = ppr->valGrad(x1, g1);
				++FuncGradEvaluations;
				fiPrime_c = -std::inner_product(dir, dir + n, g1, double{});

				if (fiPrime_c <= primeUp && fiPrime_c >= primeLow && fiC <= f0 + EPSILON)
				{
					f0 = fiC;
					++lineSearchcounter;
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
							x1[i] = x0[i] - d * dir[i];
						fiD = ppr->valGrad(x1, g1);
						++FuncGradEvaluations;
						fiPrime_d = -std::inner_product(dir, dir + n, g1, double{});

						if (fiPrime_d <= primeUp && fiPrime_d >= primeLow && fiD <= f0 + EPSILON)
						{
							f0 = fiD;
							++lineSearchcounter;
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
		}
		if (b - a > length * GAMA)
		{
			c = (a + b) / 2;
			// < Update >
			// U.begin()
			// < Termination test in c >
			// T.begin()
			for (int i = 0; i < n; i++)
				x1[i] = x0[i] - c * dir[i];
			fiC = ppr->valGrad(x1, g1);
			++FuncGradEvaluations;
			fiPrime_c = -std::inner_product(dir, dir + n, g1, double{});
			if (fiPrime_c <= primeUp && fiPrime_c >= primeLow && fiC <= f0 + EPSILON)
			{
				f0 = fiC;
				++lineSearchcounter;
				return;
			}
			// T.end()
			if (fiPrime_c < 0. && fiC > f0 + EPSILON)
			{
				b = c;	fiB = fiC;	fiPrime_b = fiPrime_c;
				// < Bisection on [a,b] >
				//B.begin()
				while (1)
				{
					d = (1 - TETA) * a + TETA * b;
					//< Termination test in d >;
					//T.begin()
					for (int i = 0; i < n; i++)
						x1[i] = x0[i] - d * dir[i];
					fiD = ppr->valGrad(x1, g1);
					++FuncGradEvaluations;
					fiPrime_d = -std::inner_product(dir, dir + n, g1, double{});
					if (fiPrime_d <= primeUp && fiPrime_d >= primeLow && fiD <= f0 + EPSILON)
					{
						f0 = fiD;
						++lineSearchcounter;
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

template<typename Problem>
void lineSearch::lineSearchConvex(Problem ppr) noexcept
{
	a = 0.;
	fiPrime_a = fiPrime0 = -std::inner_product(g0, g0 + n, dir, double{});
	double primeLow(SIGMA * fiPrime0), primeUp((2 * DELTA - 1) * fiPrime0);

	double length, t;
	//I:< Initial gess for b >
	//I.begin() 
	if (lineSearchcounter > 0)
		b = FSI2 * alpha;
	else
	{
		t = infNorm(x0, n);
		if (t != 0.)
			b = (FSI0 * t) / infNorm(g0, n);
		else
			if (f0 != 0.)
				b = (FSI0 * fabs(f0)) / std::inner_product(g0, g0 + n, g0, double{});
			else
				b = 1.;
	}
	//I.end()

	// < Generating initial [a,b]>;
	while (1)
	{
		//< Termination test in b >;
		//T.begin()
		for (int i = 0; i < n; i++)
			x1[i] = x0[i] - b * dir[i];
		fiB = ppr->valGrad(x1, g1);
		++FuncGradEvaluations;
		fiPrime_b = -std::inner_product(dir, dir + n, g1, double{});

		if (fiPrime_b <= primeUp && fiPrime_b >= primeLow && fiB <= f0 + EPSILON)
		{
			f0 = fiB;
			++lineSearchcounter;
			return;
		}
		//T.end()

		if (fiPrime_b >= 0)
			break;
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
				x1[i] = x0[i] - c * dir[i];
			fiC = ppr->valGrad(x1, g1);
			++FuncGradEvaluations;
			fiPrime_c = -std::inner_product(dir, dir + n, g1, double{});
			if (fiPrime_c <= primeUp && fiPrime_c >= primeLow && fiC <= f0 + EPSILON)
			{
				f0 = fiC;
				++lineSearchcounter;
				return;
			}
			// T.end()
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
			
			// < Update a,b>
				// U2.begin()
			if (a < c && c < b)
			{
				// < Termination test in c >
				// T.begin()
				for (int i = 0; i < n; i++)
					x1[i] = x0[i] - c * dir[i];
				fiC = ppr->valGrad(x1, g1);
				++FuncGradEvaluations;
				fiPrime_c = -std::inner_product(dir, dir + n, g1, double{});
				if (fiPrime_c <= primeUp && fiPrime_c >= primeLow && fiC <= f0 + EPSILON)
				{
					f0 = fiC;
					++lineSearchcounter;
					return;
				}
				// T.end()
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
		}
		if (b - a > length * GAMA)
		{
			c = (a + b) / 2;
			// < Update >
			// U.begin()
			// < Termination test in c >
			// T.begin()
			for (int i = 0; i < n; i++)
				x1[i] = x0[i] - c * dir[i];
			fiC = ppr->valGrad(x1, g1);
			++FuncGradEvaluations;
			fiPrime_c = -std::inner_product(dir, dir + n, g1, double{});
			if (fiPrime_c <= primeUp && fiPrime_c >= primeLow && fiC <= f0 + EPSILON)
			{
				f0 = fiC;
				++lineSearchcounter;
				return;
			}
			// T.end()
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