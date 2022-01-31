#pragma once
#include<string>
template<class T>
class problem
{
public:
	T virtual valGrad(T*, T*) = 0;
	void virtual  initialize(T* ) = 0;
	long int virtual getSize() const = 0;
	std::string virtual getName() const  = 0;
};

//-----------------------  1  -----------------------
template <class T>
class ARWHEAD : public problem<T>
{
	long int n{5000};
public:
	long int virtual getSize() const {return n;}
	std::string getName() const { return "ARWHEAD"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T fx = 0.0;
		T group1(0.0);
		g[n - 1] = 0;
		for (int i = 0; i < n - 1; i++)
		{
			group1 = x[i] * x[i] + x[n - 1] * x[n - 1];
			fx += group1 * group1 + 3.0 - 4.0 * x[i];
			g[i] = -4.0 + 4.0 * group1 * x[i];
			g[n - 1] += 4.0 * group1 * x[n - 1];
		}
		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = 1.0;
	}
};
//-----------------------  2  -----------------------
template <class T>
class BDQRTIC : public problem<T>
{
	long int n{ 1000 };
public:
	long int virtual getSize()  const {	return n;}
	std::string getName() const { return "BDQRTIC"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T fx(0.0);
		T group1(0.0), group2(0.0);

		for (int i = 0; i < n; i++)
			g[i] = 0;
		T current_squared[4] = { 0.,pow(x[0], 2),pow(x[1], 2),pow(x[2], 2) };
		T last(pow(x[n - 1], 2));

		for (int i = 0; i < n - 4; i++)
		{
			current_squared[0] = current_squared[1];
			current_squared[1] = current_squared[2];
			current_squared[2] = current_squared[3];
			current_squared[3] = x[i + 3] * x[i + 3];
			group1 = 3.0 - 4.0 * x[i];
			group2 = current_squared[0] + 2 * current_squared[1] + 3 * current_squared[2]
				+ 4 * current_squared[3] + 5 * last;
			fx += group1 * group1 + group2 * group2;
			g[i] += -8.0 * group1 + 4.0 * group2 * x[i];
			g[i + 1] += 8.0 * group2 * x[i + 1];
			g[i + 2] += 12.0 * group2 * x[i + 2];
			g[i + 3] += 16.0 * group2 * x[i + 3];
			g[n - 1] += 20.0 * group2 * x[n - 1];
		}
		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = 1.0;
	}
};
//-----------------------  3  -----------------------
//precondition: n >= 4 
template <class T>
class BROYDN7D : public problem<T>
{
	long int n{ 1000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "BROYDN7D"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T first(-2.0 * x[1] + 1 + (3. - 2.0 * x[0]) * x[0]);
		T last(-x[n - 2] + 1 + (3. - 2.0 * x[n - 1]) * x[n - 1]);
		T fx = pow(fabs(first), 7 / 3.0) + pow(fabs(last), 7 / 3.0);
		T powFabsFirst4over3, powFabsLast4over3;
		for (int i = 0; i < n; i++)
			g[i] = 0;

		g[0] = (7.0 / 3) * pow(fabs(first), 4 / 3.0) * ((first > 0) ? (1) : (-1)) * (3. - 4. * x[0]);
		g[1] = -(14.0 / 3) * pow(fabs(first), 4 / 3.0) * ((first > 0) ? (1) : (-1));
		g[n - 2] = -(7.0 / 3) * pow(fabs(last), 4 / 3.0) * ((last > 0) ? (1) : (-1));
		g[n - 1] = (7.0 / 3) * pow(fabs(last), 4 / 3.0) * ((last > 0) ? (1) : (-1)) * (3. - 4. * x[n - 1]);

		last = x[0] + x[n / 2];
		fx += pow(fabs(last), 7 / 3.0);
		g[0] += (7.0 / 3) * pow(fabs(last), 4 / 3.0) * ((last > 0) ? (1) : (-1));
		g[n / 2] += (7.0 / 3) * pow(fabs(last), 4 / 3.0) * ((last > 0) ? (1) : (-1));

		for (int i = 1; i < n / 2; i++)
		{
			first = 1 - x[i - 1] - 2.0 * x[i + 1] + (3. - 2.0 * x[i]) * x[i];
			last = x[i] + x[i + n / 2];
			fx += pow(fabs(first), 7 / 3.0) + pow(fabs(last), 7 / 3.0);

			powFabsFirst4over3 = (7.0 / 3) * pow(fabs(first), 4 / 3.0) * ((first > 0) ? (1) : (-1));
			powFabsLast4over3 = (7.0 / 3) * pow(fabs(last), 4 / 3.0) * ((last > 0) ? (1) : (-1));

			g[i - 1] += -powFabsFirst4over3;
			g[i] += powFabsFirst4over3 * (3 - 4 * x[i])
				+ powFabsLast4over3;
			g[i + 1] += -2 * powFabsFirst4over3;
			g[i + n / 2] += powFabsLast4over3;
		}
		for (int i = n / 2; i < n - 1; i++)
		{
			first = 1 - x[i - 1] - 2.0 * x[i + 1] + (3. - 2.0 * x[i]) * x[i];
			fx += pow(fabs(first), 7 / 3.0);
			powFabsFirst4over3 = (7.0 / 3) * pow(fabs(first), 4 / 3.0) * ((first > 0) ? (1) : (-1));
			g[i - 1] += -powFabsFirst4over3;
			g[i] += powFabsFirst4over3 * (3 - 4 * x[i]);
			g[i + 1] += -2 * powFabsFirst4over3;
		}
		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = 1.0;
	}
};
/*===== 4 =================== BRYBND Function ================ 5000 ===================*/
//i < n - must be
template <class T>
class BRYBND : public problem<T>
{
	long int n{ 5000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "BRYBND"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T fx = 0.0;
		T group(0.0), next(x[1] * (1 + x[1])), secGroup(next), current(x[0] * (1 + x[0]));
		for (int i = 0; i < n; i++)
			g[i] = 0;
		group = x[0] * (2.0 + 5.0 * x[0] * x[0]) + 1 - secGroup;
		fx += group * group;
		g[0] += 2.0 * group * (2.0 + 15.0 * x[0] * x[0]);
		g[1] -= 2.0 * group * (1 + 2.0 * x[1]);
		for (int i = 1; i < 6; i++)
		{
			secGroup -= next;
			secGroup += current;
			current = next;
			next = x[i + 1] * (1 + x[i + 1]);
			secGroup += next;
			group = x[i] * (2.0 + 5.0 * x[i] * x[i]) + 1 - secGroup;
			fx += group * group;
			g[i] += 2.0 * group * (2.0 + 15.0 * x[i] * x[i]);

			int  lo((0 > i - 5) ? (0) : (i - 5));
			for (int j = lo; j < i; j++)
				g[j] -= 2.0 * group * (1 + 2.0 * x[j]);
			g[i + 1] -= 2.0 * group * (1 + 2.0 * x[i + 1]);
		}
		for (int i = 6; i < n - 1; i++)
		{
			secGroup -= next;
			secGroup += current;
			current = next;
			next = x[i + 1] * (1 + x[i + 1]);
			secGroup += next;
			secGroup -= x[i - 6] * (1 + x[i - 6]);
			group = x[i] * (2.0 + 5.0 * x[i] * x[i]) + 1 - secGroup;
			fx += group * group;
			g[i] += 2.0 * group * (2.0 + 15.0 * x[i] * x[i]);

			for (int j = i - 5; j < i; j++)
				g[j] -= 2.0 * group * (1 + 2.0 * x[j]);
			g[i + 1] -= 2.0 * group * (1 + 2.0 * x[i + 1]);
		}
		secGroup -= next;
		secGroup -= x[n - 7] * (1 + x[n - 7]);
		secGroup += current;
		group = x[n - 1] * (2.0 + 5.0 * x[n - 1] * x[n - 1]) + 1 - secGroup;
		fx += group * group;
		g[n - 1] += 2.0 * group * (2.0 + 15.0 * x[n - 1] * x[n - 1]);

		for (int j = n - 6; j < n - 1; j++)
			g[j] -= 2.0 * group * (1 + 2.0 * x[j]);

		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = -1.0;
	}
};

//------------ CHAINWOO -----------  5  -----------------------
//    Precondition: n > 4
template <class T>
class CHAINWOO : public problem<T>
{
	long int n{ 1000 };
	T item1, item2, item3, item4, item5, item6;
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "CHAINWOO"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T fx = 1.0;
		for (int i = 0; i < n; i++)
			g[i] = 0;
		for (int i = 0; i < n - 3; i += 2)
		{
			item1 = (x[i + 1] - x[i] * x[i]);
			item2 = (1 - x[i]);
			item3 = (x[i + 3] - x[i + 2] * x[i + 2]);
			item4 = (1 - x[i + 2]);
			item5 = (x[i + 1] + x[i + 3] - 2.0);
			item6 = (x[i + 1] - x[i + 3]);
			fx += 100 * item1 * item1 + item2 * item2 + 90 * item3 * item3
				+ item4 * item4 + 10.0 * item5 * item5 + 0.1 * item6 * item6;
			g[i] -= (400 * item1 * x[i] + 2 * item2);
			g[i + 1] += 200 * item1 + 20.0 * item5 + 0.2 * item6;
			g[i + 2] -= (360 * item3 * x[i + 2] + 2.0 * item4);
			g[i + 3] += 180 * item3 + 20.0 * item5 - 0.2 * item6;
		}
		return fx;
	}
	void initialize(T* x)
	{
		x[0] = -3.0;
		x[1] = -1.0;
		x[2] = -3.0;
		x[3] = -1.0;
		for (int i = 4; i < n; i++)
			x[i] = -2.0;
	}
};
//----------- COSINE ------------  6  -----------------------
template <class T>
class COSINE : public problem<T>
{
	long int n{ 10000 };
	T tmp{};
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "COSINE"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T item;
		T fx(0.0);

		for (int i = 0; i < n; i++)
			g[i] = 0.0;

		for (int i = 0; i < n - 1; i++)
		{
			item = -0.5 * x[i + 1] + x[i] * x[i];
			tmp = sin(item);
			g[i] -= 2.0 * tmp * x[i];
			g[i + 1] += 0.5 * tmp;
			fx += cos(item);
		}
		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = 1.0;
	}
};
//-----------------------  7  -----------------------
template <class T>
class CRAGGLVY : public problem<T>
{
	long int n{ 5000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "CRAGGLVY"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T fx = 0.0;
		T item1, item2, element, item3, item4;
		T item1Squared, item2Squared, item3Squared, item4Squared;
		T xipow2, xipow4;
		for (int i = 0; i < n; i++)
			g[i] = 0;
		for (int i = 0; i < n - 3; i += 2)
		{
			item1 = exp(x[i]) - x[i + 1];
			item2 = x[i + 1] - x[i + 2];
			element = x[i + 2] - x[i + 3];
			item3 = tan(element) + element;
			item4 = (x[i + 3] - 1);

			item1Squared = item1 * item1;
			item2Squared = item2 * item2;
			item3Squared = item3 * item3;
			item4Squared = item4 * item4;
			xipow2 = x[i] * x[i];
			xipow4 = xipow2 * xipow2;

			fx += item1Squared * item1Squared + 100 * item2Squared * item2Squared * item2Squared
				+ item3Squared * item3Squared + xipow4 * xipow4 + item4Squared;

			g[i] += 4 * item1Squared * item1 * exp(x[i]) + 8 * x[i] * xipow4 * xipow2;
			g[i + 1] += -4 * item1 * item1Squared + 600 * item2 * item2Squared * item2Squared;
			g[i + 2] += -600 * item2 * item2Squared * item2Squared
				+ 4 * item3 * item3Squared * (1 / pow(cos(element), 2.0) + 1);

			g[i + 3] += -4 * item3 * item3Squared * (1 / pow(cos(element), 2.0) + 1)
				+ 2 * item4;

		}
		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = 2.0;
	}
};
//-----------------------  8  -----------------------
template <class T>
class CURLY10 : public problem<T>
{
	long int n{ 500 };
	int k{10};
	T fx{};
	T cube;
	int i, j;
	T q;
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "CURLY10"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		fx = 0.;
		for (int i = 0; i < n; i++)
			g[i] = 0;

		q = 0.0;
		for (int j = 0; j <= k; j++)
			q += x[j];
		fx += pow(q, 4.0) - 20 * q * q - 0.1 * q;
		for (int j = 0; j <= k; j++)
			g[j] += 4 * pow(q, 3.0) - 40 * q - 0.1;

		for (i = 1; i < n - k; i++)
		{
			q = q - x[i - 1] + x[i + k];
			T t0 = q, t1 = t0 * t0, t2 = t1 * t1;
			fx += t2 - 20 * t1 - 0.1 * t0;;
			cube = 4 * t0 * t1 - 40 * t0 - 0.1;
			j = i;
			g[j] += cube; ++j;
			for (; j <= i + k; )
			{
				g[j] += cube; ++j;
				g[j] += cube; ++j;
				g[j] += cube; ++j;
				g[j] += cube; ++j;
				g[j] += cube; ++j;
			}
		}
		for (; i < n; i++)
		{
			q -= x[i - 1];
			T t0 = q, t1 = t0 * t0, t2 = t1 * t1;
			fx += t2 - 20 * t1 - 0.1 * t0;;
			cube = 4 * t0 * t1 - 40 * t0 - 0.1;
			for (int j = i; j <= n - 1; j++)
				g[j] += cube;
		}
		return fx;
	}

	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = 0.0001 / (n + 1);
	}
};
//-----------------------  11  ----------------------- 5000
template <class T>
class DIXMAANA : public problem<T>
{
	long int n{ 5000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "DIXMAANA"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T fx = 0.0;
		T item(0.0);
		int m(n / 3);
		for (int i = 0; i < n; i++)
			g[i] = 0;
		for (int i = 0; i < m; i++)
		{
			item = x[i] * pow(x[i + m], 2);
			fx += pow(x[i], 2) + 0.125 * pow(item, 2) + 0.125 * x[i] * x[i + 2 * m];
			g[i] += 2 * x[i] + 0.25 * item * pow(x[i + m], 2) + 0.125 * x[i + 2 * m];
			g[i + m] += 0.5 * item * x[i] * x[i + m];
			g[i + 2 * m] += 0.125 * x[i];
		}
		for (int i = m; i < 2 * m; i++)
		{
			item = x[i] * pow(x[i + m], 2);
			fx += pow(x[i], 2) + 0.125 * pow(item, 2);
			g[i] += 2 * x[i] + 0.25 * item * pow(x[i + m], 2);
			g[i + m] += 0.5 * item * x[i] * x[i + m];
		}
		for (int i = 2 * m; i < n; i++)
		{
			fx += pow(x[i], 2);
			g[i] += 2 * x[i];
		}
		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = 2.0;
	}
};

//-----------------------  15  -----------------------
template <class T>
class DIXMAANE : public problem<T> 
{
	long int n{ 3000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "DIXMAANE"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T fx = 0.0;
		T item(0.0);
		int m(n / 3);
		for (int i = 0; i < n; i++)
			g[i] = 0;
		for (int i = 0; i < m; i++)
		{
			item = x[i] * pow(x[i + m], 2);
			fx += (pow(x[i], 2) * (i + 1)) / n + 0.125 * pow(item, 2) + (0.125 * x[i] * x[i + 2 * m] * (i + 1)) / n;
			g[i] += (2 * x[i] * (i + 1)) / n + 0.25 * item * pow(x[i + m], 2) + (0.125 * x[i + 2 * m] * (i + 1)) / n;
			g[i + m] += 0.5 * item * x[i] * x[i + m];
			g[i + 2 * m] += (0.125 * x[i] * (i + 1)) / n;
		}
		for (int i = m; i < 2 * m; i++)
		{
			item = x[i] * pow(x[i + m], 2);
			fx += (pow(x[i], 2) * (i + 1)) / n + 0.125 * pow(item, 2);
			g[i] += (2 * x[i] * (i + 1)) / n + 0.25 * item * pow(x[i + m], 2);
			g[i + m] += 0.5 * item * x[i] * x[i + m];
		}
		for (int i = 2 * m; i < n; i++)
		{
			fx += (pow(x[i], 2) * i) / n;
			g[i] += (2 * x[i] * (i + 1)) / n;
		}
		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = 2.0;
	}
};

//-----------------------  19  -----------------------
template <class T>
class DIXMAANI : public problem<T>
{
	long int n{ 3000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "DIXMAANI"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T fx = 0.0;
		T item(0.0);
		int m(n / 3);
		for (int i = 0; i < n; i++)
			g[i] = 0;
		for (int i = 0; i < m; i++)
		{
			item = x[i] * pow(x[i + m], 2);
			fx += pow(x[i], 2) * pow((i + 1) / T(n), 2) + 0.125 * pow(item, 2)
				+ 0.125 * x[i] * x[i + 2 * m] * pow((i + 1) / T(n), 2);
			g[i] += 2 * x[i] * pow((i + 1) / T(n), 2) + 0.25 * item * pow(x[i + m], 2)
				+ 0.125 * x[i + 2 * m] * pow((i + 1) / T(n), 2);
			g[i + m] += 0.5 * item * x[i] * x[i + m];
			g[i + 2 * m] += 0.125 * x[i] * pow((i + 1) / T(n), 2);
		}
		for (int i = m; i < 2 * m; i++)
		{
			item = x[i] * pow(x[i + m], 2);
			fx += pow(x[i], 2) * pow((i + 1) / T(n), 2) + 0.125 * pow(item, 2);
			g[i] += 2 * x[i] * pow((i + 1) / T(n), 2) + 0.25 * item * pow(x[i + m], 2);
			g[i + m] += 0.5 * item * x[i] * x[i + m];
		}
		for (int i = 2 * m; i < n; i++)
		{
			fx += pow(x[i], 2) * pow((i + 1) / T(n), 2);
			g[i] += 2 * x[i] * pow((i + 1) / T(n), 2);
		}
		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = 2.0;
	}
};

//-----------------------  23  ----------------------- 
template <class T>
class DIXON3DQ : public problem<T>
{
	long int n{ 5000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "DIXON3DQ"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T fx = 0.0;
		T item{};

		item = x[0] - 1;
		fx += item * item;
		g[0] = 2.0 * item;

		for (int i = 1; i < n - 1; i++)
			g[i] = 0.0;

		item = x[n - 1] - 1;
		fx += item * item;
		g[n - 1] = 2.0 * item;


		for (int i = 1; i < n - 1; i++)
		{
			item = (x[i] - x[i + 1]);
			fx += item * item;
			g[i] += 2.0 * item;
			g[i + 1] -= 2.0 * item;
		}
		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = -1.0;
	}
};

//-----------------------  24  -----------------------
template <class T>
class DQDRTIC : public problem<T>
{
	long int n{ 5000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "DQDRTIC"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		for (int i = 0; i < n; i++)
			g[i] = 0;
		T t0 = x[0] * x[0], t1 = x[1] * x[1], t2 = x[2] * x[2];
		T fx = t0 + 100 * (t1 + t2);
		g[0] += 2 * x[0];
		g[1] += 200 * x[1];
		g[2] += 200 * x[2];
		for (int i = 1; i < n - 2; i++)
		{
			t0 = t1;
			t1 = t2;
			t2 = x[i + 2] * x[i + 2];
			fx += t0 + 100 * (t1 + t2);
			g[i] += 2 * x[i];
			g[i + 1] += 200 * x[i + 1];
			g[i + 2] += 200 * x[i + 2];
		}
		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = 3.0;
	}
};

//-----------------------  25  -----------------------
template <class T>
class DQRTIC : public problem<T>
{
	long int n{ 5000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "DQRTIC"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T fx = 0.0;
		T item, squaredItem;
		for (int i = 0; i < n; i++)
		{
			item = x[i] - i - 1;
			squaredItem = item * item;
			fx += squaredItem * squaredItem;
			g[i] = 4 * squaredItem * item;
		}
		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = 2.0;
	}
};

//-----------------------  26  -----------------------
template <class T>
class EDENSCH : public problem<T>
{
	long int n{ 2000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "EDENSCH"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T fx = 0.0;
		for (int i = 0; i < n; i++)
			g[i] = 0;
		T item1, item2, item3;
		T squaredItem1;
		for (int i = 0; i < n - 1; i++)
		{
			item1 = x[i] - 2;
			item2 = x[i] * x[i + 1] - 2 * x[i + 1];
			item3 = x[i + 1] + 1;
			squaredItem1 = item1 * item1;
			fx += 16 + squaredItem1 * squaredItem1 + item2 * item2 + item3 * item3;
			g[i] += 4 * squaredItem1 * item1 + 2 * item2 * x[i + 1];
			g[i + 1] += 2 * item2 * (x[i] - 2.0) + 2 * item3;
		}
		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = .0;
	}
};

//-----------------------  27  -----------------------
template <class T>
class EG2 : public problem<T>
{
	long int n{ 5000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "EG2"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		for (int i = 0; i < n; i++)
			g[i] = 0;
		T fx = 0.5 * sin(pow(x[n - 1], 2));
		g[n - 1] = cos(pow(x[n - 1], 2)) * x[n - 1];
		T item;
		for (int i = 0; i < n - 1; i++)
		{
			item = x[0] + x[i] * x[i] - 1;;
			fx += sin(item);
			g[0] += cos(item);
			g[i] += 2 * cos(item) * x[i];
		}
		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = .0;
	}
};

//-----------------------  28  -----------------------
template <class T>
class ENGVAL1 : public problem<T>
{
	long int n{ 5000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "ENGVAL1"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T fx = 0.0;
		T item;
		for (int i = 0; i < n; i++)
			g[i] = 0.0;

		for (int i = 0; i < n - 1; i++)
		{
			item = x[i] * x[i] + x[i + 1] * x[i + 1];
			fx += item * item + (3 - 4.0 * x[i]);
			g[i] += 4.0 * item * x[i] - 4.0;
			g[i + 1] += 4.0 * item * x[i + 1];
		}
		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = 2.0;
	}
};

//-----------------------  29  -----------------------
template <class T>
class EXTROSNB : public problem<T>
{
	long int n{ 1000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "EXTROSNB"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T item(x[0] - 1);
		T fx(pow(item, 2));
		for (int i = 1; i < n; i++)
			g[i] = 0.0;
		g[0] = 2 * item;
		for (int i = 1; i < n; i++) {
			item = 10 * (x[i] - x[i - 1] * x[i - 1]);
			fx += item * item;
			g[i - 1] += -40.0 * item * x[i - 1];
			g[i] += 20.0 * item;
		}
		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = -1.0;
	}
};

//-----------------------  30  -----------------------
template <class T>
class FLETCHR : public problem<T>
{
	long int n{ 1000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "FLETCHR"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T item;
		T fx = 0.0;
		for (int i = 0; i < n; i++)
			g[i] = 0;

		for (int i = 0; i < n - 1; i++)
		{
			item = x[i + 1] - x[i] + 1 - x[i] * x[i];
			fx += item * item;
			g[i] += 20.0 * item * (-2.0 * x[i] - 1.0);
			g[i + 1] += 20.0 * item;
		}
		return 100. * fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = .0;
	}
};


//-----------------------  31  -----------------------
template <class T>
class FREUROTH : public problem<T>
{
	long int n{ 5000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "FREUROTH"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T fx = 0.0;
		T item1, item2;;
		for (int i = 0; i < n; i++)
			g[i] = 0;
		for (int i = 0; i < n - 1; i++)
		{
			item1 = (-13 + x[i] + ((5 - x[i + 1]) * x[i + 1] - 2.0) * x[i + 1]);
			item2 = (-29 + x[i] + ((1 + x[i + 1]) * x[i + 1] - 14.0) * x[i + 1]);
			fx += item1 * item1 + item2 * item2;
			g[i] += 2.0 * item1 + 2.0 * item2;
			g[i + 1] += 2.0 * item1 * (10 * x[i + 1] - 3.0 * x[i + 1] * x[i + 1] - 2.0) +
				2.0 * item2 * (2 * x[i + 1] + 3.0 * x[i + 1] * x[i + 1] - 14.0);
		}
		return fx;
	}
	void initialize(T* x)
	{
		x[0] = 0.5;
		x[1] = -2.0;
		for (int i = 2; i < n; i++)
			x[i] = 0.0;
	}
};

//-----------------------  32  -----------------------
template <class T>
class GENHUMPS : public problem<T>
{
	long int n{ 5000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "GENHUMPS"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		for (int i = 0; i < n; i++)
			g[i] = 0;
		T item1 = sin(2.0 * x[0]);
		T item2 = sin(2.0 * x[1]);
		T item11 = item1 * item1;
		T item22 = item2 * item2;
		T t0 = x[0] * x[0], t1 = x[1] * x[1];
		T fx = item11 * item22 + 0.05 * (t0 + t1);
		g[0] = 4.0 * item1 * cos(2.0 * x[0]) * item22 + 0.1 * x[0];
		g[1] = 4.0 * item2 * cos(2.0 * x[1]) * item11 + 0.1 * x[1];

		for (int i = 1; i < n - 1; i++)
		{

			item1 = item2;
			item2 = sin(2.0 * x[i + 1]);
			item11 = item22;
			item22 = item2 * item2;
			t0 = t1;
			t1 = x[i + 1] * x[i + 1];
			fx += item11 * item22 + 0.05 * (t0 + t1);
			t0 = 4.0 * item1 * item2;
			g[i] += t0 * cos(2.0 * x[i]) * item2 + 0.1 * x[i];
			g[i + 1] += t0 * cos(2.0 * x[i + 1]) * item1 + 0.1 * x[i + 1];
		}
		return fx;
	}
	void initialize(T* x)
	{
		x[0] = -506.0;
		for (int i = 1; i < n; i++)
			x[i] = 506.2;
	}
};

//-----------------------  33  -----------------------
template <class T>
class GENROSE : public problem<T>
{
	long int n{ 1000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "GENROSE"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T fx = 1.0;
		T item1, item2;
		for (int i = 0; i < n; i++)
			g[i] = 0;
		for (int i = 1; i < n; i++)
		{
			item1 = 10.0 * (x[i] - x[i - 1] * x[i - 1]);
			item2 = x[i] - 1.0;
			fx += item1 * item1 + item2 * item2;
			g[i - 1] -= 40.0 * item1 * x[i - 1];
			g[i] += 20 * item1 + 2 * item2;
		}
		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = 1.0 / n + 1;
	}
};

//-----------------------  34  -----------------------
template <class T>
class LIARWDH : public problem<T>
{
	long int n{ 5000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "LIARWDH"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T item1, item2;
		T fx(0.0);
		for (int i = 0; i < n; i++)
		{
			item1 = 2.0 * (x[i] * x[i] - x[0]);
			item2 = (x[i] - 1);
			fx += item1 * item1 + item2 * item2;
			g[i] = 8.0 * item1 * x[i] + 2.0 * item2;
			g[0] -= 4.0 * item1;
		}
		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = 4.0;
	}
};

//-----------------------  35  -----------------------
template <class T>
class MOREBV : public problem<T>
{
	long int n{ 5000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "MOREBV"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T h(1.0 / n);
		T element, item;
		T fx(0.0);
		for (int i = 0; i < n; i++)
			g[i] = 0;
		//first term
		element = x[1] + h + 1;
		item = 2 * x[1] - x[2] + (h * h * pow(element, 3)) / 2;
		fx += item * item;
		g[1] = 2 * item * (2.0 + h * h * 1.5 * pow(element, 2));
		g[2] = -2 * item;
		//last term
		element = x[n - 1] + (n - 1) * h + 1;
		item = 2 * x[n - 1] - x[n - 2] + (h * h * pow(element, 3)) / 2;
		fx += item * item;
		g[n - 1] = 4 * item * (2.0 + h * h * 1.5 * pow(element, 2));
		g[n - 2] = -2 * item;
		for (int i = 2; i < n - 1; i++)
		{
			element = x[i] + i * h + 1;
			item = 2 * x[i] - x[i - 1] - x[i + 1] + (h * h * pow(element, 3)) / 2;
			fx += item * item;
			g[i - 1] -= 2 * item;
			g[i] += 2 * item * (2.0 + h * h * 1.5 * pow(element, 2));
			g[i + 1] -= 2 * item;
		}
		return fx;
	}
	void initialize(T* x)
	{
		T h(1.0 / n);
		x[0] = 0.;
		for (int i = 1; i < n; i++)
			x[i] = i * h * (i * h - 1.0);
	}
};

//-----------------------  36  -----------------------
template <class T>
class NONCVXU2 : public problem<T>
{
	long int n{ 1000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "NONCVXU2"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T fx = 0.0;
		for (int i = 0; i < n; i++)
			g[i] = 0;
		T item1, item2;
		for (int j = 0; j < n; j++)
		{
			int i(j + 1);
			item1 = (x[j] + x[(3 * i - 2) % n] + x[(7 * i - 3) % n]);
			item2 = x[j] + x[(3 * i - 2) % n] + x[(7 * i - 3) % n];
			fx += item1 * item1 + 4 * cos(item2);
			g[j] += 2 * item1 - 4 * sin(item2);
			g[(3 * i - 2) % n] += 2 * item1 - 4 * sin(item2);
			g[(7 * i - 3) % n] += 2 * item1 - 4 * sin(item2);
		}
		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = i + 1;
	}
};

//-----------------------  37  -----------------------
template <class T>
class NONDIA : public problem<T>
{
	long int n{ 10000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "NONDIA"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T item(x[0] - 1.0);
		T fx(item * item);

		g[0] = 2.0 * item;
		item = 10.0 * (x[0] - x[0] * x[0]);
		fx += item * item;
		g[0] += (20.0 - 40.0 * x[0]) * item;

		for (int i = 1; i < n; i++)
		{
			item = 10.0 * (x[0] - x[i] * x[i]);
			fx += item * item;
			g[0] += 20.0 * item;
			g[i] = -40.0 * x[i] * item;
		}
		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = -1.0;
	}
};

//-----------------------  38  -----------------------
template <class T>
class NONDQUAR : public problem<T>
{
	long int n{ 1000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "NONQUAR"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T item, tmp;
		T fx(0.0);

		for (int i = 0; i < n; i++)
			g[i] = 0.0;

		item = x[0] - x[1];
		fx += item * item;
		g[0] = 2.0 * item;
		g[1] = -2.0 * item;

		item = x[n - 2] + x[n - 1];
		fx += item * item;
		g[n - 2] = 2.0 * item;
		g[n - 1] = 2.0 * item;

		for (int i = 0; i < n - 2; i++)
		{
			T t0 = x[i] + x[i + 1] + x[n - 1], t1 = t0 * t0, t2 = t1 * t1, t3 = t0 * t1;

			fx += t2;
			tmp = 4.0 * t3;
			g[i] += tmp;
			g[i + 1] += tmp;
			g[n - 1] += tmp;
		}
		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = -1.0;
	}
};

//-----------------------  39  -----------------------
template <class T>
class PENALTY1 : public problem<T>
{
	long int n{ 1000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "PENALTY1"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T item, tmp;
		T tail(0.0);
		T a(1E-5);

		T fx = 0.0;

		for (int i = 0; i < n; i++)
		{
			item = x[i] - 1;
			fx += item * item;
			g[i] = 2.0 * a * item;
			tail += x[i] * x[i];
		}
		fx *= a;

		fx += pow(tail - 0.25, 2);
		tmp = 4 * (tail - 0.25);
		for (int i = 0; i < n; i++)
			g[i] += tmp * x[i];
		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = i + 1;
	}
};

//-----------------------  40  -----------------------
template <class T>
class PENALTY2 : public problem<T>
{
	long int n{ 100 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "PENALTY2"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T item1, item2;
		T tail(0.0);
		T a(1E-5);
		T a1(0.2 * a);
		const T ExpMinus1by10 = exp(-1 / 10.0);
		for (int i = 0; i < n; i++)
			g[i] = 0.0;
		item1 = x[0] - 0.2;
		T fx(pow(item1, 2));
		g[0] = 2 * item1;
		tail += (n)*pow(x[0], 2);
		T currentTerm = exp(x[0] / 10);
		T prevTer;
		T currentI = exp(1 / 10.0);
		T prevI;
		for (int i = 1; i < n; i++)
		{
			prevTer = currentTerm;
			currentTerm = exp(x[i] / 10);
			prevI = currentI;
			currentI = exp((i + 1) / 10.0);
			item1 = currentTerm + prevTer - currentI - prevI;
			item2 = currentTerm - ExpMinus1by10;
			tail += (n - i) * x[i] * x[i];
			fx += a * (item1 * item1 + item2 * item2);
			g[i - 1] += a1 * item1 * prevTer;
			g[i] += a1 * currentTerm * (item1 + item2);
		}
		fx += (tail - 1) * (tail - 1);
		for (int i = 0; i < n; i++)
			g[i] += 4 * (tail - 1) * (n - i) * x[i];
		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = 0.5;
	}
};

//-----------------------  41  -----------------------
template <class T>
class POWER : public problem<T>
{
	long int n{ 1000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "POWER"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T item;
		T fx(0.0);

		for (int i = 0; i < n; i++)
		{
			item = (i + 1) * x[i];
			fx += item * item;
			g[i] = 2.0 * item * (i + 1);
		}
		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = 1.0;
	}
};
	
//-----------------------  42  -----------------------
template <class T>
class SROSENBR : public problem<T>
{
	long int n{ 10000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "SROSENBR"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		int i;
		T fx = 0.0;

		for (i = 0; i < n; i += 2) {
			T t1 = 1.0 - x[i];
			T t2 = 10.0 * (x[i + 1] - x[i] * x[i]);
			g[i + 1] = 20.0 * t2;
			g[i] = -2.0 * (x[i] * g[i + 1] + t1);
			fx += t1 * t1 + t2 * t2;
		}
		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i += 2)
		{
			x[i] = -1.2;
			x[i + 1] = 1.0;
		}
	}
};

//-----------------------  43  -----------------------
template <class T>
class TRIDIA : public problem<T>
{
	long int n{ 10000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "TRIDIA"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T item(x[0] - 1.0);
		T fx(item * item);
		for (int i = 0; i < n; i++)
			g[i] = 0;
		g[0] = 2.0 * item;

		for (int i = 1; i < n; i++)
		{
			item = (2.0 * x[i] - x[i - 1]);
			fx += (i + 1) * item * item;
			g[i] += 4.0 * item * (i + 1);
			g[i - 1] -= 2.0 * (i + 1) * item;
		}
		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = 1.0;
	}
};

//-----------------------  44  -----------------------
template <class T>
class Woods : public problem<T>
{
	long int n{ 5000 };
public:
	long int virtual getSize()  const { return n; }
	std::string getName() const { return "Woods"; }
	T valGrad
	(
		T* x,
		T* g
	)
	{
		T fx(0.0);
		T item1, item2, item3, item4, item5, item6;

		for (int i = 0; i < n; i += 4)
		{
			item1 = (x[i + 1] - x[i] * x[i]);
			item2 = (1 - x[i]);
			item3 = (x[i + 3] - x[i + 2] * x[i + 2]);
			item4 = (1 - x[i + 2]);
			item5 = (x[i + 1] + x[i + 3] - 2.0);
			item6 = (x[i + 1] - x[i + 3]);
			fx += (100 * item1 * item1 + item2 * item2 + 90 * item3 * item3
				+ item4 * item4 + 10.0 * item5 * item5 + 0.1 * item6 * item6);
			g[i] = -400 * item1 * x[i] - 2 * item2;
			g[i + 1] = 200 * item1 + 20.0 * item5 + 0.2 * item6;
			g[i + 2] = -360 * item3 * x[i + 2] - 2.0 * item4;
			g[i + 3] = 180 * item3 + 20.0 * item5 - 0.2 * item6;
		}
		return fx;
	}
	void initialize(T* x)
	{
		for (int i = 0; i < n; i += 2)
		{
			x[i] = -3;
			x[i + 1] = -1;
		}
	}
};

//--------------------------------  Benchmarking program for long reals  ---------------------------------------

#include<vector>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include"solver.h"
using boost::multiprecision::cpp_dec_float_50;
using namespace std;


//vector < unique_ptr<problem<double>>> d_problems;
vector<unique_ptr<problem<cpp_dec_float_50>>> dec_float_50_problems;

void makeFoat50TestsVector(void)
{
	dec_float_50_problems.emplace_back(make_unique<ARWHEAD<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<BDQRTIC<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<BROYDN7D<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<BRYBND<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<CHAINWOO<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<COSINE<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<CRAGGLVY<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<CURLY10<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<DIXMAANA<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<DIXMAANE<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<DIXMAANI<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<DIXON3DQ<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<DQDRTIC<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<DQRTIC<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<EDENSCH<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<EG2<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<ENGVAL1<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<EXTROSNB<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<FLETCHR<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<FREUROTH<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<GENHUMPS<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<GENROSE<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<LIARWDH<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<MOREBV<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<NONCVXU2<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<NONDIA<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<NONDQUAR<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<PENALTY1<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<PENALTY2<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<POWER<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<SROSENBR<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<TRIDIA<cpp_dec_float_50>>());
	dec_float_50_problems.emplace_back(make_unique<Woods<cpp_dec_float_50>>());
}
 
void runUCONFoat50Tests()
{
	using
		real = cpp_dec_float_50;

	bool what{};
	int repNumber(1);
	makeFoat50TestsVector();
	vector<unique_ptr<problem<cpp_dec_float_50>>>& v = dec_float_50_problems;
	std::cout <<"Test Problems: " << std::endl;
	for(size_t i=0; i<v.size(); ++i)
		std::cout << i+1 << ":   Problem:  " << v[i]->getName() << ",  size: " << v[i]->getSize() << endl;

	std::cout << "Enter test's index between  "
		<< 1 << "  and  " << v.size() << std::endl;
	int i;
	cin >> i;
	if (1 > i || i > v.size())
	{
		std::cout << "Wrong input!" << std::endl;
		return;
	}
	--i;
	vector<_int64> repetitions(repNumber);
	freopen("results.txt", "w", stdout);

	for (int j = 0; j < repNumber; j++)
	{
		auto st = chrono::high_resolution_clock::now();
		MHB<real> minimizer{ 1e-6,v[i].get() };
		v[i]->initialize(minimizer.getVariables());
		what = minimizer.solve();

		auto diff = chrono::high_resolution_clock::now() - st;
		auto time = chrono::duration_cast<chrono::microseconds>(diff);
		repetitions[j] = time.count();

		if (j == repNumber - 1)
		{
			sort(repetitions.begin(), repetitions.end());
			std::cout << "Problem:  " << v[i]->getName() << ",  size: " << v[i]->getSize() << std::endl;
			std::cout << "Average time:  " << repetitions[repNumber / 2] << " microseconds" << std::endl;
			if (what)
				minimizer.printVector(minimizer.getVariables(), "x");
		}
	}
}

//--------------------------------  Benchmarking program for doubles  ---------------------------------------
vector<unique_ptr<problem<double>>> d_prec_problems;

void makeDoubleTestsVector(void)
{
	d_prec_problems.emplace_back(make_unique<ARWHEAD<double>>());
	d_prec_problems.emplace_back(make_unique<BDQRTIC<double>>());
	d_prec_problems.emplace_back(make_unique<BROYDN7D<double>>());
	d_prec_problems.emplace_back(make_unique<BRYBND<double>>());
	d_prec_problems.emplace_back(make_unique<CHAINWOO<double>>());
	d_prec_problems.emplace_back(make_unique<COSINE<double>>());
	d_prec_problems.emplace_back(make_unique<CRAGGLVY<double>>());
	d_prec_problems.emplace_back(make_unique<CURLY10<double>>());
	d_prec_problems.emplace_back(make_unique<DIXMAANA<double>>());
	d_prec_problems.emplace_back(make_unique<DIXMAANE<double>>());
	d_prec_problems.emplace_back(make_unique<DIXMAANI<double>>());
	d_prec_problems.emplace_back(make_unique<DIXON3DQ<double>>());
	d_prec_problems.emplace_back(make_unique<DQDRTIC<double>>());
	d_prec_problems.emplace_back(make_unique<DQRTIC<double>>());
	d_prec_problems.emplace_back(make_unique<EDENSCH<double>>());
	d_prec_problems.emplace_back(make_unique<EG2<double>>());
	d_prec_problems.emplace_back(make_unique<ENGVAL1<double>>());
	d_prec_problems.emplace_back(make_unique<EXTROSNB<double>>());
	d_prec_problems.emplace_back(make_unique<FLETCHR<double>>());
	d_prec_problems.emplace_back(make_unique<FREUROTH<double>>());
	d_prec_problems.emplace_back(make_unique<GENHUMPS<double>>());
	d_prec_problems.emplace_back(make_unique<GENROSE<double>>());
	d_prec_problems.emplace_back(make_unique<LIARWDH<double>>());
	d_prec_problems.emplace_back(make_unique<MOREBV<double>>());
	d_prec_problems.emplace_back(make_unique<NONCVXU2<double>>());
	d_prec_problems.emplace_back(make_unique<NONDIA<double>>());
	d_prec_problems.emplace_back(make_unique<NONDQUAR<double>>());
	d_prec_problems.emplace_back(make_unique<PENALTY1<double>>());
	d_prec_problems.emplace_back(make_unique<PENALTY2<double>>());
	d_prec_problems.emplace_back(make_unique<POWER<double>>());
	d_prec_problems.emplace_back(make_unique<SROSENBR<double>>());
	d_prec_problems.emplace_back(make_unique<TRIDIA<double>>());
	d_prec_problems.emplace_back(make_unique<Woods<double>>());
}

void runUCONDoubleTests()
{
	using
		real = double;

	bool what{};
	int repNumber(1);
	makeDoubleTestsVector();
	vector<unique_ptr<problem<double>>>& v = d_prec_problems;
	std::cout << "Test Problems: " << std::endl;
	for (size_t i = 0; i < v.size(); ++i)
		std::cout << i + 1 << ":   Problem:  " << v[i]->getName() << ",  size: " << v[i]->getSize() << endl;

	std::cout << "Enter test's index between  "
		<< 1 << "  and  " << v.size() << std::endl;
	int i;
	cin >> i;
	if (1 > i || i > v.size())
	{
		std::cout << "Wrong input!" << std::endl;
		return;
	}
	--i;
	vector<_int64> repetitions(repNumber);
	freopen("results.txt", "w", stdout);

	for (int j = 0; j < repNumber; j++)
	{
		auto st = chrono::high_resolution_clock::now();
		MHB<real> minimizer{ 1e-6,v[i].get() };
		v[i]->initialize(minimizer.getVariables());
		what = minimizer.solve();

		auto diff = chrono::high_resolution_clock::now() - st;
		auto time = chrono::duration_cast<chrono::microseconds>(diff);
		repetitions[j] = time.count();

		if (j == repNumber - 1)
		{
			sort(repetitions.begin(), repetitions.end());
			std::cout << "Problem:  " << v[i]->getName() << ",  size: " << v[i]->getSize() << std::endl;
			std::cout << "Average time:  " << repetitions[repNumber / 2] << " microseconds" << std::endl;
			if (!what)
				std::cout << "Solved Patially... gradient inf norm: " 
				          << minimizer.infNorm(minimizer.getGradient()) << std::endl;
			minimizer.printVector(minimizer.getVariables(), "x");
		}
	}
}

//--------------------------------  Benchmarking: the hybrid approach -----------------------------

void runUCONHybridTests()
{
	bool what{};
	int repNumber(1);
	makeDoubleTestsVector();
	vector<unique_ptr<problem<double>>>& v = d_prec_problems;

	makeFoat50TestsVector();
	vector<unique_ptr<problem<cpp_dec_float_50>>>& u = dec_float_50_problems;

	std::cout << "Test Problems: " << std::endl;
	for (size_t i = 0; i < v.size(); ++i)
		std::cout << i + 1 << ":   Problem:  " << v[i]->getName() << ",  size: " << v[i]->getSize() << endl;

	std::cout << "Enter test's index between  "
		<< 1 << "  and  " << v.size() << std::endl;
	int i;
	cin >> i;
	if (1 > i || i > v.size())
	{
		std::cout << "Wrong input!" << std::endl;
		return;
	}
	--i;
	vector<_int64> repetitions(repNumber);
	freopen("results.txt", "w", stdout);

	for (int j = 0; j < repNumber; j++)
	{
		auto st = chrono::high_resolution_clock::now();
		MHB<double> minimizer{ 1e-6,v[i].get() };
		v[i]->initialize(minimizer.getVariables());
		what = minimizer.solve();

		if (what)
		{
			auto diff = chrono::high_resolution_clock::now() - st;
			auto time = chrono::duration_cast<chrono::microseconds>(diff);
			repetitions[j] = time.count();
			if (j == repNumber - 1)
			{
				sort(repetitions.begin(), repetitions.end());
				std::cout << "Problem:  " << v[i]->getName() << ",  size: " << v[i]->getSize() << std::endl;
				std::cout << "Average time:  " << repetitions[repNumber / 2] << " microseconds" << std::endl;
				std::cout << "Gradient inf norm: " << minimizer.infNorm(minimizer.getGradient()) << std::endl;
				minimizer.printVector(minimizer.getVariables(), "x");
			}
			continue;
		}

		if (j == repNumber - 1)
		{
			std::cout << "Precision was changed at.. gradient inf norm: " << minimizer.infNorm(minimizer.getGradient()) << std::endl;
			std::cout << "Objective function:  " << minimizer.getValue() <<
					",   Value-Grad evaluations:  " << minimizer.getEvaluations() <<
					",   Restarts (line searches):  " << minimizer.getRestarts() << std::endl << std::endl;
		}

		MHB<cpp_dec_float_50> longFloatMinimizer{ 1e-6,u[i].get() };
		long int size = minimizer.getSize();
		auto double_var_s = minimizer.getVariables();
		auto long_float_var_s = longFloatMinimizer.getVariables();
		for (long int k = 0; k < size; ++k)
			long_float_var_s[i] = double_var_s[i];

		longFloatMinimizer.solve();
	
		auto diff = chrono::high_resolution_clock::now() - st;
		auto time = chrono::duration_cast<chrono::microseconds>(diff);
		repetitions[j] = time.count();

		if (j == repNumber - 1)
		{
			sort(repetitions.begin(), repetitions.end());
			std::cout << "Problem:  " << v[i]->getName() << ",  size: " << v[i]->getSize() << std::endl;
			std::cout << "Average time:  " << repetitions[repNumber / 2] << " microseconds" << std::endl;
			longFloatMinimizer.printVector(longFloatMinimizer.getVariables(), "x");
		}
	}
}