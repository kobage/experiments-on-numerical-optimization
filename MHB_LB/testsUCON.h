#pragma once
#include<iostream>
#include"chrono"
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include"lineSearch.h"
using namespace std;
using namespace u_min;

/*===== 0 =================== ARWHEAD Function ================ 5000 ===================*/
double ARWHEAD
(
	double* x,
	double* g,
	const int n
)
{
	double fx = 0.0;
	double group1(0.0);
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
void InitializeARWHEAD
(
	double* x,
	const int n
)
{
	for (int i = 0; i < n; i++)
		x[i] = 1.0;
}
/*===== 1 ================= BDQRTIC Function ============= 5000 ==================*/
double BDQRTIC
(
	double* x,
	double* g,
	const int n
)
{
	double fx(0.0);
	double group1(0.0), group2(0.0);

	for (int i = 0; i < n; i++)
		g[i] = 0;
	double current_squared[4] = { 0.,pow(x[0], 2),pow(x[1], 2),pow(x[2], 2) };
	double last(pow(x[n - 1], 2));

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
void InitializeBDQRTIC
(
	double* x,
	const int n
)
{
	for (int i = 0; i < n; i++)
		x[i] = 1.0;
}
/*===== 2 ================= BROYDN7D Function ============= 5000==================*/
//precondition: n>=4
double BROYDN7D
(
	double* x,
	double* g,
	const int n
)
{
	double first(-2.0 * x[1] + 1 + (3. - 2.0 * x[0]) * x[0]);
	double last(-x[n - 2] + 1 + (3. - 2.0 * x[n - 1]) * x[n - 1]);
	double fx = pow(fabs(first), 7 / 3.0) + pow(fabs(last), 7 / 3.0);
	double powFabsFirst4over3, powFabsLast4over3;
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
void InitializeBROYDN7D
(
	double* x,
	const int n
)
{
	for (int i = 0; i < n; i++)
		x[i] = 1.0;
}
/*===== 3 =================== BRYBND Function ================ 5000 ===================*/
//i < n - must be
double BRYBND
(
	double* x,
	double* g,
	const int n
)
{
	double fx = 0.0;
	double group(0.0), next(x[1] * (1 + x[1])), secGroup(next), current(x[0] * (1 + x[0]));
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
void InitializeBRYBND
(
	double* x,
	const int n
)
{
	for (int i = 0; i < n; i++)
		x[i] = -1.0;
}
/*========= 4 ========= Chained Wood - CHAINWOO ======== 4000 ==========*/
double CHAINWOO
(
	double* x,
	double* g,
	const int n
)
{
	double fx = 1.0;
	double item1, item2, item3, item4, item5, item6;
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
void InitialCHAINWOO
(
	double* x,
	const int n
)
{
	x[0] = -3.0;
	x[1] = -1.0;
	x[2] = -3.0;
	x[3] = -1.0;
	for (int i = 4; i < n; i++)
		x[i] = -2.0;
}
/*========= 5 ========== COSINE Function =========== 10000 ===========*/
double COSINE
(
	double* x,
	double* g,
	const int n
)
{
	double item;
	double fx(0.0);

	for (int i = 0; i < n; i++)
		g[i] = 0.0;

	for (int i = 0; i < n - 1; i++)
	{
		item = -0.5 * x[i + 1] + x[i] * x[i];
		fx += cos(item);
		g[i] -= 2.0 * sin(item) * x[i];
		g[i + 1] += 0.5 * sin(item);
	}
	return fx;
}
void InitialCOSINE
(
	double* x,
	const int n
)
{
	for (int i = 0; i < n; i++)
		x[i] = 1.0;
}
/*========= 6 ========= CRAGGLVY  =========== 5000 =============*/
double CRAGGLVY
(
	double* x,
	double* g,
	const int n
)
{
	double fx = 1.0;
	double item1, item2, element, item3, item4;
	double item1Squared, item2Squared, item3Squared, item4Squared;
	double xipow2, xipow4;
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
void InitialCRAGGLVY
(
	double* x,
	const int n
)
{
	x[0] = 1.0;

	for (int i = 1; i < n; i++)
		x[i] = 2.0;
}
/*===== 7 =============== CURLY10 Function ============= 500 ===================*/
double CURLY10
(
	double* x,
	double* g,
	const int n
)
{
	double fx = 0.0;
	double cube;
	int k(10);
	int i, j;
	double q;
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
		double t0 = q, t1 = t0 * t0, t2 = t1 * t1;
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
		double t0 = q, t1 = t0 * t0, t2 = t1 * t1;
		fx += t2 - 20 * t1 - 0.1 * t0;;
		cube = 4 * t0 * t1 - 40 * t0 - 0.1;
		for (int j = i; j <= n - 1; j++)
			g[j] += cube;
	}
	return fx;
}
void InitializeCURLY
(
	double* x,
	const int n
)
{
	for (int i = 0; i < n; i++)
		x[i] = 0.0001 / (n + 1);
}
/*===== 8 ============= CURLY20 Function ============= 500 ===================*/
double CURLY20
(
	double* x,
	double* g,
	const int n
)
{
	double fx = 0.0;
	double cube;
	int k(20);
	int i, j;
	double q;
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
		double t0 = q, t1 = t0 * t0, t2 = t1 * t1;
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
		double t0 = q, t1 = t0 * t0, t2 = t1 * t1;
		fx += t2 - 20 * t1 - 0.1 * t0;;
		cube = 4 * t0 * t1 - 40 * t0 - 0.1;
		for (int j = i; j <= n - 1; j++)
			g[j] += cube;
	}
	return fx;
}
//initialization in curly functions is common
/*===== 9 ============= CURLY30 Function ======== 500 =================*/
double CURLY30
(
	double* x,
	double* g,
	const int n
)
{
	double fx = 0.0;
	double cube;
	int k(30);
	int i, j;
	double q;
	for (int i = 0; i < n; i++)
		g[i] = 0;

	q = 0.0;
	for (int j = 0; j <= k; j++)
		q += x[j];
	fx += pow(q, 4.0) - 20 * q * q - 0.1 * q;
	cube = 4 * pow(q, 3.0) - 40 * q - 0.1;
	for (j = 0; j <= k; j++)
		g[j] += cube;

	for (i = 1; i < n - k; i++)
	{
		q = q - x[i - 1] + x[i + k];
		double t0 = q, t1 = t0 * t0, t2 = t1 * t1;
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
		double t0 = q, t1 = t0 * t0, t2 = t1 * t1;
		fx += t2 - 20 * t1 - 0.1 * t0;;
		cube = 4 * t0 * t1 - 40 * t0 - 0.1;
		for (j = i; j <= n - 1; j++)
			g[j] += cube;
	}
	return fx;
}

/*===== 10 ======================= DIXMAANA Function ================== 5000 ================*/
double DIXMAANA
(
	double* x,
	double* g,
	const int n
)
{
	double fx = 1.0;
	double item(0.0);
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
void InitializeDIXMAAN
(
	double* x,
	const int n
)
{
	for (int i = 0; i < n; i++)
		x[i] = 2.0;
}

/*===== 11 ======================= DIXMAANB Function ================== 5000 ================*/
double DIXMAANB
(
	double* x,
	double* g,
	const int n
)
{
	double fx = 1.0;
	double item1(0.0);
	double item2(0.0);
	int m(n / 3);
	for (int i = 0; i < n; i++)
		g[i] = 0;
	for (int i = 0; i < m; i++)
	{
		item1 = 0.25 * x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = 0.25 * x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) + item1 * item1 + item2 * item2 + (0.0625 * x[i] * x[i + 2 * m]);
		g[i] += 2 * x[i] + 0.5 * item1 * (x[i + 1] + x[i + 1] * x[i + 1]) + 0.5 * item2 * pow(x[i + m], 2) +
			(0.0625 * x[i + 2 * m]);
		g[i + 1] += 0.5 * item1 * x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += item2 * x[i] * x[i + m];
		g[i + 2 * m] += 0.0625 * x[i];
	}
	for (int i = m; i < 2 * m; i++)
	{
		item1 = 0.25 * x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = 0.25 * x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) + item1 * item1 + item2 * item2;
		g[i] += 2 * x[i] + 0.5 * item1 * (x[i + 1] + x[i + 1] * x[i + 1]) + 0.5 * item2 * pow(x[i + m], 2);
		g[i + 1] += 0.5 * item1 * x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += item2 * x[i] * x[i + m];
	}
	for (int i = 2 * m; i < n - 1; i++)
	{
		item1 = 0.25 * x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += pow(x[i], 2) + item1 * item1;
		g[i] += 2 * x[i] + 0.5 * item1 * (x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.5 * item1 * x[i] * (1 + 2 * x[i + 1]);
	}
	fx += pow(x[n - 1], 2);
	g[n - 1] += 2 * x[n - 1];
	return fx;
}
/*=====12  ======================= DIXMAANC Function ================== 5000 ================*/
double DIXMAANC
(
	double* x,
	double* g,
	const int n
)
{
	double fx = 1.0;
	double item1(0.0);
	double item2(0.0);
	int m(n / 3);
	for (int i = 0; i < n; i++)
		g[i] = 0;
	for (int i = 0; i < m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) + 0.125 * item1 * item1 + 0.125 * item2 * item2 + 0.125 * x[i] * x[i + 2 * m];
		g[i] += 2 * x[i] + 0.25 * item1 * (x[i + 1] + x[i + 1] * x[i + 1]) + 0.25 * item2 * pow(x[i + m], 2) + 0.125 * x[i + 2 * m];
		g[i + 1] += 0.25 * item1 * x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 0.5 * item2 * x[i] * x[i + m];
		g[i + 2 * m] += 0.125 * x[i];
	}
	for (int i = m; i < 2 * m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) + 0.125 * item1 * item1 + 0.125 * item2 * item2;
		g[i] += 2 * x[i] + 0.25 * item1 * (x[i + 1] + x[i + 1] * x[i + 1]) + 0.25 * item2 * pow(x[i + m], 2);
		g[i + 1] += 0.25 * item1 * x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 0.5 * item2 * x[i] * x[i + m];
	}
	for (int i = 2 * m; i < n - 1; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += pow(x[i], 2) + 0.125 * item1 * item1;
		g[i] += 2 * x[i] + 0.25 * item1 * (x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.25 * item1 * x[i] * (1 + 2 * x[i + 1]);
	}
	fx += pow(x[n - 1], 2);
	g[n - 1] += 2 * x[n - 1];
	return fx;
}
/*===== 13 ======================= DIXMAAND Function ================== 5000 ================*/
double DIXMAAND
(
	double* x,
	double* g,
	const int n
)
{
	double fx = 1.0;
	double item1(0.0);
	double item2(0.0);
	int m(n / 3);
	for (int i = 0; i < n; i++)
		g[i] = 0;
	for (int i = 0; i < m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) + 0.26 * item1 * item1 + 0.26 * item2 * item2 + 0.26 * x[i] * x[i + 2 * m];
		g[i] += 2 * x[i] + 0.52 * item1 * (x[i + 1] + x[i + 1] * x[i + 1]) + 0.52 * item2 * pow(x[i + m], 2) + 0.26 * x[i + 2 * m];
		g[i + 1] += 0.52 * item1 * x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 1.04 * item2 * x[i] * x[i + m];
		g[i + 2 * m] += 0.26 * x[i];
	}
	for (int i = m; i < 2 * m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) + 0.26 * item1 * item1 + 0.26 * item2 * item2;
		g[i] += 2 * x[i] + 0.52 * item1 * (x[i + 1] + x[i + 1] * x[i + 1]) + 0.52 * item2 * pow(x[i + m], 2);
		g[i + 1] += 0.52 * item1 * x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 1.04 * item2 * x[i] * x[i + m];
	}
	for (int i = 2 * m; i < n - 1; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += pow(x[i], 2) + 0.26 * item1 * item1;
		g[i] += 2 * x[i] + 0.52 * item1 * (x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.52 * item1 * x[i] * (1 + 2 * x[i + 1]);
	}
	fx += pow(x[n - 1], 2);
	g[n - 1] += 2 * x[n - 1];
	return fx;
}
/*===== 14 ======================= DIXMAANE Function ================== 3000 ================*/
double DIXMAANE
(
	double* x,
	double* g,
	const int n
)
{
	double fx = 1.0;
	double item(0.0);
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
/*===== 15 ======================= DIXMAANF Function ================== 3000 ================*/
double DIXMAANF
(
	double* x,
	double* g,
	const int n
)
{
	double fx = 1.0;
	double item1(0.0);
	double item2(0.0);
	int m(n / 3);
	for (int i = 0; i < n; i++)
		g[i] = 0;
	for (int i = 0; i < m; i++)
	{
		item1 = 0.25 * x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = 0.25 * x[i] * pow(x[i + m], 2);
		fx += (pow(x[i], 2) * (i + 1)) / n + item1 * item1 + item2 * item2 + (0.0625 * x[i] * x[i + 2 * m] * (i + 1)) / n;
		g[i] += (2 * x[i] * (i + 1)) / n + 0.5 * item1 * (x[i + 1] + x[i + 1] * x[i + 1]) + 0.5 * item2 * pow(x[i + m], 2) +
			(0.0625 * x[i + 2 * m] * (i + 1)) / n;
		g[i + 1] += 0.5 * item1 * x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += item2 * x[i] * x[i + m];
		g[i + 2 * m] += (0.0625 * x[i] * (i + 1)) / n;
	}
	for (int i = m; i < 2 * m; i++)
	{
		item1 = 0.25 * x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = 0.25 * x[i] * pow(x[i + m], 2);
		fx += (pow(x[i], 2) * (i + 1)) / n + item1 * item1 + item2 * item2;
		g[i] += (2 * x[i] * (i + 1)) / n + 0.5 * item1 * (x[i + 1] + x[i + 1] * x[i + 1]) + 0.5 * item2 * pow(x[i + m], 2);
		g[i + 1] += 0.5 * item1 * x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += item2 * x[i] * x[i + m];
	}
	for (int i = 2 * m; i < n - 1; i++)
	{
		item1 = 0.25 * x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += (pow(x[i], 2) * (i + 1)) / n + item1 * item1;
		g[i] += (2 * x[i] * (i + 1)) / n + 0.5 * item1 * (x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.5 * item1 * x[i] * (1 + 2 * x[i + 1]);
	}
	fx += pow(x[n - 1], 2);
	g[n - 1] += 2 * x[n - 1];
	return fx;
}
/*===== 16 ======================= DIXMAANG Function ================== 5000 ================*/
double DIXMAANG
(
	double* x,
	double* g,
	const int n
)
{
	double fx = 1.0;
	double item1(0.0);
	double item2(0.0);
	int m(n / 3);
	for (int i = 0; i < n; i++)
		g[i] = 0;
	for (int i = 0; i < m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += (pow(x[i], 2) * (i + 1)) / n + 0.125 * item1 * item1 + 0.125 * item2 * item2
			+ (0.125 * x[i] * x[i + 2 * m] * (i + 1)) / n;
		g[i] += (2 * x[i] * (i + 1)) / n + 0.25 * item1 * (x[i + 1] + x[i + 1] * x[i + 1]) + 0.25 * item2 * pow(x[i + m], 2)
			+ (0.125 * x[i + 2 * m] * (i + 1)) / n;
		g[i + 1] += 0.25 * item1 * x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 0.5 * item2 * x[i] * x[i + m];
		g[i + 2 * m] += (0.125 * x[i] * (i + 1)) / n;
	}
	for (int i = m; i < 2 * m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += (pow(x[i], 2) * (i + 1)) / n + 0.125 * item1 * item1 + 0.125 * item2 * item2;
		g[i] += (2 * x[i] * (i + 1)) / n + 0.25 * item1 * (x[i + 1] + x[i + 1] * x[i + 1]) + 0.25 * item2 * pow(x[i + m], 2);
		g[i + 1] += 0.25 * item1 * x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 0.5 * item2 * x[i] * x[i + m];
	}
	for (int i = 2 * m; i < n - 1; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += (pow(x[i], 2) * (i + 1)) / n + 0.125 * item1 * item1;
		g[i] += (2 * x[i] * (i + 1)) / n + 0.25 * item1 * (x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.25 * item1 * x[i] * (1 + 2 * x[i + 1]);
	}
	fx += pow(x[n - 1], 2);
	g[n - 1] += 2 * x[n - 1];
	return fx;
}

/*===== 17 ======================= DIXMAANH Function ================== 5000 ================*/
double DIXMAANH
(
	double* x,
	double* g,
	const int n
)
{
	double fx = 1.0;
	double item1(0.0);
	double item2(0.0);
	int m(n / 3);
	for (int i = 0; i < n; i++)
		g[i] = 0;
	for (int i = 0; i < m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += (pow(x[i], 2) * (i + 1)) / n + 0.26 * item1 * item1 + 0.26 * item2 * item2
			+ (0.26 * x[i] * x[i + 2 * m] * (i + 1)) / n;
		g[i] += (2 * x[i] * (i + 1)) / n + 0.52 * item1 * (x[i + 1] + x[i + 1] * x[i + 1]) + 0.52 * item2 * pow(x[i + m], 2)
			+ (0.26 * x[i + 2 * m] * (i + 1)) / n;
		g[i + 1] += 0.52 * item1 * x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 1.04 * item2 * x[i] * x[i + m];
		g[i + 2 * m] += (0.26 * x[i] * (i + 1)) / n;
	}
	for (int i = m; i < 2 * m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += (pow(x[i], 2) * (i + 1)) / n + 0.26 * item1 * item1 + 0.26 * item2 * item2;
		g[i] += (2 * x[i] * (i + 1)) / n + 0.52 * item1 * (x[i + 1] + x[i + 1] * x[i + 1]) + 0.52 * item2 * pow(x[i + m], 2);
		g[i + 1] += 0.52 * item1 * x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 1.04 * item2 * x[i] * x[i + m];
	}
	for (int i = 2 * m; i < n - 1; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += (pow(x[i], 2) * (i + 1)) / n + 0.26 * item1 * item1;
		g[i] += (2 * x[i] * (i + 1)) / n + 0.52 * item1 * (x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.52 * item1 * x[i] * (1 + 2 * x[i + 1]);
	}
	fx += pow(x[n - 1], 2);
	g[n - 1] += 2 * x[n - 1];
	return fx;
}
/*===== 18 ======================= DIXMAANI Function ================== 3000 ================*/
double DIXMAANI
(
	double* x,
	double* g,
	const int n
)
{
	double fx = 1.0;
	double item(0.0);
	int m(n / 3);
	for (int i = 0; i < n; i++)
		g[i] = 0;
	for (int i = 0; i < m; i++)
	{
		item = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) * pow((i + 1) / double(n), 2) + 0.125 * pow(item, 2)
			+ 0.125 * x[i] * x[i + 2 * m] * pow((i + 1) / double(n), 2);
		g[i] += 2 * x[i] * pow((i + 1) / double(n), 2) + 0.25 * item * pow(x[i + m], 2)
			+ 0.125 * x[i + 2 * m] * pow((i + 1) / double(n), 2);
		g[i + m] += 0.5 * item * x[i] * x[i + m];
		g[i + 2 * m] += 0.125 * x[i] * pow((i + 1) / double(n), 2);
	}
	for (int i = m; i < 2 * m; i++)
	{
		item = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) * pow((i + 1) / double(n), 2) + 0.125 * pow(item, 2);
		g[i] += 2 * x[i] * pow((i + 1) / double(n), 2) + 0.25 * item * pow(x[i + m], 2);
		g[i + m] += 0.5 * item * x[i] * x[i + m];
	}
	for (int i = 2 * m; i < n; i++)
	{
		fx += pow(x[i], 2) * pow((i + 1) / double(n), 2);
		g[i] += 2 * x[i] * pow((i + 1) / double(n), 2);
	}
	return fx;
}
/*===== 19 ======================= DIXMAANJ Function ================== 3000 ================*/
double DIXMAANJ
(
	double* x,
	double* g,
	const int n
)
{
	double fx = 1.0;
	double item1(0.0);
	double item2(0.0);
	int m(n / 3);
	for (int i = 0; i < n; i++)
		g[i] = 0;
	for (int i = 0; i < m; i++)
	{
		item1 = 0.25 * x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = 0.25 * x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) * pow((i + 1) / double(n), 2) + item1 * item1 + item2 * item2
			+ 0.0625 * x[i] * x[i + 2 * m] * pow((i + 1) / double(n), 2);
		g[i] += 2 * x[i] * pow((i + 1) / double(n), 2) + 0.5 * item1 * (x[i + 1] + x[i + 1] * x[i + 1])
			+ 0.5 * item2 * pow(x[i + m], 2) + 0.0625 * x[i + 2 * m] * pow((i + 1) / double(n), 2);
		g[i + 1] += 0.5 * item1 * x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += item2 * x[i] * x[i + m];
		g[i + 2 * m] += 0.0625 * x[i] * pow((i + 1) / double(n), 2);
	}
	for (int i = m; i < 2 * m; i++)
	{
		item1 = 0.25 * x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = 0.25 * x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) * pow((i + 1) / double(n), 2) + item1 * item1 + item2 * item2;
		g[i] += 2 * x[i] * pow((i + 1) / double(n), 2) + 0.5 * item1 * (x[i + 1] + x[i + 1] * x[i + 1])
			+ 0.5 * item2 * pow(x[i + m], 2);
		g[i + 1] += 0.5 * item1 * x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += item2 * x[i] * x[i + m];
	}
	for (int i = 2 * m; i < n - 1; i++)
	{
		item1 = 0.25 * x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += pow(x[i], 2) * pow((i + 1) / double(n), 2) + item1 * item1;
		g[i] += 2 * x[i] * pow((i + 1) / double(n), 2) + 0.5 * item1 * (x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.5 * item1 * x[i] * (1 + 2 * x[i + 1]);
	}
	fx += pow(x[n - 1], 2);
	g[n - 1] += 2 * x[n - 1];
	return fx;
}
/*===== 20 ======================= DIXMAANK Function ================== 5000 ================*/
double DIXMAANK
(

	double* x,
	double* g,
	const int n
)
{
	double fx = 1.0;
	double item1(0.0);
	double item2(0.0);
	int m(n / 3);
	for (int i = 0; i < n; i++)
		g[i] = 0;
	for (int i = 0; i < m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) * pow((i + 1) / double(n), 2) + 0.125 * item1 * item1 + 0.125 * item2 * item2
			+ 0.125 * x[i] * x[i + 2 * m] * pow((i + 1) / double(n), 2);
		g[i] += 2 * x[i] * pow((i + 1) / double(n), 2) + 0.25 * item1 * (x[i + 1] + x[i + 1] * x[i + 1])
			+ 0.25 * item2 * pow(x[i + m], 2) + 0.125 * x[i + 2 * m] * pow((i + 1) / double(n), 2);
		g[i + 1] += 0.25 * item1 * x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 0.5 * item2 * x[i] * x[i + m];
		g[i + 2 * m] += 0.125 * x[i] * pow((i + 1) / double(n), 2);
	}
	for (int i = m; i < 2 * m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) * pow((i + 1) / double(n), 2) + 0.125 * item1 * item1 + 0.125 * item2 * item2;
		g[i] += 2 * x[i] * pow((i + 1) / double(n), 2) + 0.25 * item1 * (x[i + 1] + x[i + 1] * x[i + 1])
			+ 0.25 * item2 * pow(x[i + m], 2);
		g[i + 1] += 0.25 * item1 * x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 0.5 * item2 * x[i] * x[i + m];
	}
	for (int i = 2 * m; i < n - 1; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += pow(x[i], 2) * pow((i + 1) / double(n), 2) + 0.125 * item1 * item1;
		g[i] += 2 * x[i] * pow((i + 1) / double(n), 2) + 0.25 * item1 * (x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.25 * item1 * x[i] * (1 + 2 * x[i + 1]);
	}
	fx += pow(x[n - 1], 2);
	g[n - 1] += 2 * x[n - 1];
	return fx;
}

/*===== 21 =================== DIXMAANL Function ================== 5000 ================*/
double DIXMAANL
(
	double* x,
	double* g,
	const int n
)
{
	double fx = 1.0;
	double item1(0.0);
	double item2(0.0);
	int m(n / 3);
	for (int i = 0; i < n; i++)
		g[i] = 0;
	for (int i = 0; i < m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) * pow((i + 1) / double(n), 2) + 0.26 * item1 * item1 + 0.26 * item2 * item2
			+ 0.26 * x[i] * x[i + 2 * m] * pow((i + 1) / double(n), 2);
		g[i] += 2 * x[i] * pow((i + 1) / double(n), 2) + 0.52 * item1 * (x[i + 1] + x[i + 1] * x[i + 1])
			+ 0.52 * item2 * pow(x[i + m], 2) + 0.26 * x[i + 2 * m] * pow((i + 1) / double(n), 2);
		g[i + 1] += 0.52 * item1 * x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 1.04 * item2 * x[i] * x[i + m];
		g[i + 2 * m] += 0.26 * x[i] * pow((i + 1) / double(n), 2);
	}
	for (int i = m; i < 2 * m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) * pow((i + 1) / double(n), 2) + 0.26 * item1 * item1 + 0.26 * item2 * item2;
		g[i] += 2 * x[i] * pow((i + 1) / double(n), 2) + 0.52 * item1 * (x[i + 1] + x[i + 1] * x[i + 1])
			+ 0.52 * item2 * pow(x[i + m], 2);
		g[i + 1] += 0.52 * item1 * x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 1.04 * item2 * x[i] * x[i + m];
	}
	for (int i = 2 * m; i < n - 1; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += pow(x[i], 2) * pow((i + 1) / double(n), 2) + 0.26 * item1 * item1;
		g[i] += 2 * x[i] * pow((i + 1) / double(n), 2) + 0.52 * item1 * (x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.52 * item1 * x[i] * (1 + 2 * x[i + 1]);
	}
	fx += pow(x[n - 1], 2);
	g[n - 1] += 2 * x[n - 1];
	return fx;
}
/*========= 22 ========== DIXON3DQ Function ============= 5000 ==============*/
double DIXON3DQ
(
	double* x,
	double* g,
	const int n
)
{
	double item;
	double fx(0.0);

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
void InitialDIXON3DQ
(
	double* x,
	const int n
)
{
	for (int i = 0; i < n; i++)
		x[i] = -1.0;
}
/*===== 23 ===================== DQDRTIC Function ===== 10000 ======================*/

double DQDRTIC
(
	double* x,
	double* g,
	const int n
)
{
	for (int i = 0; i < n; i++)
		g[i] = 0;
	double t0 = x[0] * x[0], t1 = x[1] * x[1], t2 = x[2] * x[2];
	double fx = t0 + 100 * (t1 + t2);
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
void InitialDQDRTIC
(
	double* x,
	const int n
)
{
	for (int i = 0; i < n; i++)
		x[i] = 3.0;
}

/*===== 24 ===================== DQRTIC Function ===== 5000 ======================*/
double DQRTIC
(
	double* x,
	double* g,
	const int n
)
{
	double fx = 0.0;
	double item, squaredItem;
	for (int i = 0; i < n; i++)
	{
		item = x[i] - i - 1;
		squaredItem = item * item;
		fx += squaredItem * squaredItem;
		g[i] = 4 * squaredItem * item;
	}
	return fx;
}
void InitialDQRTIC
(
	double* x,
	const int n
)
{
	for (int i = 0; i < n; i++)
		x[i] = 2.0;
}
/*===== 25 ===================== EDENSCH Function ===== 2000 ======================*/
double EDENSCH
(
	double* x,
	double* g,
	const int n
)
{
	double fx = 0.0;
	for (int i = 0; i < n; i++)
		g[i] = 0;
	double item1, item2, item3;
	double squaredItem1;
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
void InitialEDENSCH
(
	double* x,
	const int n
)
{
	for (int i = 0; i < n; i++)
		x[i] = 0.0;
}
/*===== 26 ===================== EG2 Function ===== 5000 ======================*/
double EG2
(
	double* x,
	double* g,
	const int n
)
{
	for (int i = 0; i < n; i++)
		g[i] = 0;
	double fx = 0.5 * sin(pow(x[n - 1], 2));
	g[n - 1] = cos(pow(x[n - 1], 2)) * x[n - 1];
	double item;
	for (int i = 0; i < n - 1; i++)
	{
		item = x[0] + x[i] * x[i] - 1;;
		fx += sin(item);
		g[0] += cos(item);
		g[i] += 2 * cos(item) * x[i];
	}
	return fx;
}
void InitialEG2
(
	double* x,
	const int n
)
{
	for (int i = 0; i < n; i++)
		x[i] = 0.;
}
/*========= 27 ========== ENGVAL1 Function ============= 5000 ==============*/
double ENGVAL1
(
	double* x,
	double* g,
	const int n
)
{
	double item;
	double fx(0.0);

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
void InitialENGVAL1
(
	double* x,
	const int n
)
{
	for (int i = 0; i < n; i++)
		x[i] = 2.0;

}
/*===== 28 ======================= EXTROSNB ================== 1000 ================*/
double EXTROSNB
(
	double* x,
	double* g,
	const int n
)
{
	double item(x[0] - 1);
	double fx(pow(item, 2));
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
void InitialEXTROSNB
(
	double* x,
	const int n
)
{
	for (int i = 0; i < n; i++)
		x[i] = -1.0;
}
/*========= 29 ========== FLETCHR Function ============= 1000 ==============*/
double FLETCHR
(
	double* x,
	double* g,
	const int n
)
{
	double item;
	double fx = 0.0;
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
void InitialFLETCHR
(
	double* x,
	const int n
)
{
	for (int i = 0; i < n; i++)
		x[i] = 0.0;

}
/*========= 30 ==============  Freudenstein and Roth: FREUROTH ====== 5000 =======*/
double FREUROTH
(
	double* x,
	double* g,
	const int n
)
{
	double fx = 0.0;
	double item1, item2;;
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
void InitialFREUROTH
(
	double* x,
	const int n
)
{
	x[0] = 0.5;
	x[1] = -2.0;
	for (int i = 2; i < n; i++)
		x[i] = 0.0;
}

/*========= 31 ===================  GENHUMPS ====== 5000 ======================*/
double GENHUMPS
(
	double* x,
	double* g,
	const int n
)
{

	for (int i = 0; i < n; i++)
		g[i] = 0;
	double item1 = sin(2.0 * x[0]);
	double item2 = sin(2.0 * x[1]);
	double item11 = item1 * item1;
	double item22 = item2 * item2;
	double t0 = x[0] * x[0], t1 = x[1] * x[1];
	double fx = item11 * item22 + 0.05 * (t0 + t1);
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
void InitialGENHUMPS
(
	double* x,
	const int n
)
{
	x[0] = -506.0;
	for (int i = 1; i < n; i++)
		x[i] = 506.2;
}
/*===== 32 ======================= GENROSE ====================== 1000 ===================*/
double GENROSE
(
	double* x,
	double* g,
	const int n
)
{
	double fx = 1.0;
	//double item1, item2;
	double item1, item2;
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
void InitialGENROSE
(
	double* x,
	const int n
)
{
	for (int i = 0; i < n; i++)
		x[i] = 1.0 / n + 1;
}
/*========= 33 ========== LIARWDH Function ============= 5000 ==============*/
double LIARWDH
(
	double* x,
	double* g,
	const int n
)
{
	double item1, item2;
	double fx(0.0);
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
void InitialLIARWDH
(
	double* x,
	const int n
)
{
	for (int i = 0; i < n; i++)
		x[i] = 4.0;
}
/*========= 34 ========== MOREBV ============= 5000 ==============*/
double MOREBV
(
	double* x,
	double* g,
	const int n
)
{
	double h(1.0 / n);
	double element, item;
	double fx(0.0);
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
void InitialMOREBV
(
	double* x,
	const int n
)
{
	double h(1.0 / n);
	x[0] = 0.;
	for (int i = 1; i < n; i++)
		x[i] = i * h * (i * h - 1.0);
}
/*===== 35 ======================= NONCVXU2 ============ 1000 ===============*/

double NONCVXU2
(
	double* x,
	double* g,
	const int n
)
{
	double fx = 0.0;
	for (int i = 0; i < n; i++)
		g[i] = 0;
	double item1, item2;
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
void InitialNONCVXU2
(
	double* x,
	const int n
)
{
	for (int i = 0; i < n; i++)
		x[i] = i + 1;
}
/*========= 36 ========== NONDIA Function ============= 10000 ==============*/
double NONDIA
(
	double* x,
	double* g,
	const int n
)
{
	double item(x[0] - 1.0);
	double fx(item * item);

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
void InitialNONDIA
(
	double* x,
	const int n
)
{
	for (int i = 0; i < n; i++)
		x[i] = -1.0;

}

/*========= 37 ========== NONDQUAR Function ============= 10000 ==============*/
double NONDQUAR
(
	double* x,
	double* g,
	const int n
)
{
	double item, tmp;
	double fx(0.0);

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
		double t0 = x[i] + x[i + 1] + x[n - 1], t1 = t0 * t0, t2 = t1 * t1, t3 = t0 * t1;

		fx += t2;
		tmp = 4.0 * t3;
		g[i] += tmp;
		g[i + 1] += tmp;
		g[n - 1] += tmp;
	}
	return fx;
}
void InitialNONDQUAR
(
	double* x,
	const int n
)
{
	for (int i = 0; i < n; i += 2)
	{
		x[i] = 1.0;
		x[i + 1] = -1.0;
	}

}
/*========= 38 ========== PENALTY1 ============= 1000 ==============*/
double PENALTY1
(
	double* x,
	double* g,
	const int n
)
{
	double item, tmp;
	double tail(0.0);
	double a(1E-5);

	double fx = 0.0;

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
void InitialPENALTY1
(
	double* x,
	const int n
)
{
	for (int i = 0; i < n; i++)
		x[i] = i + 1;

}
/*========= 39 ========== PENALTY2 ============= 100 ==============*/
double PENALTY2
(
	double* x,
	double* g,
	const int n
)
{
	double item1, item2;
	double tail(0.0);
	double a(1E-5);
	double a1(0.2 * a);
	const double ExpMinus1by10 = exp(-1 / 10.0);
	for (int i = 0; i < n; i++)
		g[i] = 0.0;
	item1 = x[0] - 0.2;
	double fx(pow(item1, 2));
	g[0] = 2 * item1;
	tail += (n)*pow(x[0], 2);
	double currentTerm = exp(x[0] / 10);
	double prevTer;
	double currentI = exp(1 / 10.0);
	double prevI;
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

void InitialPENALTY2
(
	double* x,
	const int n
)
{
	for (int i = 0; i < n; i++)
		x[i] = 0.5;

}
/*========= 40 ========== POWER Function ======== 1000 ==========*/
double POWER
(
	double* x,
	double* g,
	const int n
)
{
	double item;
	double fx(0.0);

	for (int i = 0; i < n; i++)
	{
		item = (i + 1) * x[i];
		fx += item * item;
		g[i] = 2.0 * item * (i + 1);
	}
	return fx;
}
void InitialPOWER
(
	double* x,
	const int n
)
{
	for (int i = 0; i < n; i++)
		x[i] = 1.0;

}
/*===== 41 =============== SROSENBR ============= 10000 ================*/
double SROSENBR
(
	double* x,
	double* g,
	const int n
)
{
	int i;
	double fx = 0.0;

	for (i = 0; i < n; i += 2) {
		double t1 = 1.0 - x[i];
		double t2 = 10.0 * (x[i + 1] - x[i] * x[i]);
		g[i + 1] = 20.0 * t2;
		g[i] = -2.0 * (x[i] * g[i + 1] + t1);
		fx += t1 * t1 + t2 * t2;
	}
	return fx;
}
void InitialSROSENBR
(
	double* x,
	const int n
)
{
	for (int i = 0; i < n; i += 2)
	{
		x[i] = -1.2;
		x[i + 1] = 1.0;
	}
}

/*========= 42 ========== TRIDIA Function ============= 10000 ==============*/
double TRIDIA
(
	double* x,
	double* g,
	const int n
)
{
	double item(x[0] - 1.0);
	double fx(item * item);
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
void InitialTRIDIA
(
	double* x,
	const int n
)
{
	for (int i = 0; i < n; i++)
		x[i] = 1.0;
}
/*========= 43 =========  Woods ======== 10000 =====================*/
double WOODS
(
	double* x,
	double* g,
	const int n
)
{
	double fx(0.0);
	double item1, item2, item3, item4, item5, item6;

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
void InitialWOODS
(
	double* x,
	const int n
)
{
	for (int i = 0; i < n; i += 2)
	{
		x[i] = -3;
		x[i + 1] = -1;
	}
}

//Some tools to work with tests
struct testFunction
{
	double(*fName)(double*, double*, const int);
	void(*Initializer)(double*, const int);

	int size;
	string testName;
	void show(void)
	{
		cout << "Name:  " << testName << endl;
		cout << "Size:  " << size << endl;
	}
};
vector<testFunction> testsVector;
void makeTestsVector(void)
{
	testFunction tmpTest;
	//0
	tmpTest.fName = ARWHEAD;
	tmpTest.Initializer = InitializeARWHEAD;
	tmpTest.testName = "ARWHEAD function";
	tmpTest.size = 5000;
	testsVector.push_back(tmpTest);
	//1
	tmpTest.fName = BDQRTIC;
	tmpTest.Initializer = InitializeBDQRTIC;
	tmpTest.testName = "BDQRTIC function";
	tmpTest.size = 5000;
	testsVector.push_back(tmpTest);
	//2
	tmpTest.fName = BROYDN7D;
	tmpTest.Initializer = InitializeBROYDN7D;
	tmpTest.testName = "BROYDN7D function";
	tmpTest.size = 5000;
	testsVector.push_back(tmpTest);
	//3
	tmpTest.fName = BRYBND;
	tmpTest.Initializer = InitializeBRYBND;
	tmpTest.testName = "BRYBND function";
	tmpTest.size = 5000;
	testsVector.push_back(tmpTest);
	//4
	tmpTest.fName = CHAINWOO;
	tmpTest.Initializer = InitialCHAINWOO;
	tmpTest.testName = "CHAINWOO function";
	tmpTest.size = 4000;
	testsVector.push_back(tmpTest);
	//5
	tmpTest.fName = COSINE;
	tmpTest.Initializer = InitialCOSINE;
	tmpTest.testName = "COSINE function";
	tmpTest.size = 10000;
	testsVector.push_back(tmpTest);
	//6
	tmpTest.fName = CRAGGLVY;
	tmpTest.Initializer = InitialCRAGGLVY;
	tmpTest.testName = "CRAGGLVY";
	tmpTest.size = 5000;
	testsVector.push_back(tmpTest);
	//7
	tmpTest.fName = CURLY10;
	tmpTest.Initializer = InitializeCURLY;
	tmpTest.testName = "CURLY10";
	tmpTest.size = 500;
	testsVector.push_back(tmpTest);
	//8
	tmpTest.fName = CURLY20;
	tmpTest.Initializer = InitializeCURLY;
	tmpTest.testName = "CURLY20";
	tmpTest.size = 500;
	testsVector.push_back(tmpTest);
	//9
	tmpTest.fName = CURLY30;
	tmpTest.Initializer = InitializeCURLY;
	tmpTest.testName = "CURLY30";
	tmpTest.size = 500;
	testsVector.push_back(tmpTest);
	//10
	tmpTest.fName = DIXMAANA;
	tmpTest.Initializer = InitializeDIXMAAN;
	tmpTest.testName = "DIXMAANA function";
	tmpTest.size = 5000;
	testsVector.push_back(tmpTest);
	//11
	tmpTest.fName = DIXMAANB;
	tmpTest.Initializer = InitializeDIXMAAN;
	tmpTest.testName = "DIXMAANB  function";
	tmpTest.size = 5000;
	testsVector.push_back(tmpTest);
	//12
	tmpTest.fName = DIXMAANC;
	tmpTest.Initializer = InitializeDIXMAAN;
	tmpTest.testName = "DIXMAANC  function";
	tmpTest.size = 5000;
	testsVector.push_back(tmpTest);
	//13
	tmpTest.fName = DIXMAAND;
	tmpTest.Initializer = InitializeDIXMAAN;
	tmpTest.testName = "DIXMAAND  function";
	tmpTest.size = 5000;
	testsVector.push_back(tmpTest);
	//14
	tmpTest.fName = DIXMAANE;
	tmpTest.Initializer = InitializeDIXMAAN;
	tmpTest.testName = "DIXMAANE  function";
	tmpTest.size = 3000;
	testsVector.push_back(tmpTest);
	//15
	tmpTest.fName = DIXMAANF;
	tmpTest.Initializer = InitializeDIXMAAN;
	tmpTest.testName = "DIXMAANF  function";
	tmpTest.size = 3000;
	testsVector.push_back(tmpTest);
	//16
	tmpTest.fName = DIXMAANG;
	tmpTest.Initializer = InitializeDIXMAAN;
	tmpTest.testName = "DIXMAANG  function";
	tmpTest.size = 3000;
	testsVector.push_back(tmpTest);
	//17
	tmpTest.fName = DIXMAANH;
	tmpTest.Initializer = InitializeDIXMAAN;
	tmpTest.testName = "DIXMAANH function";
	tmpTest.size = 3000;
	testsVector.push_back(tmpTest);
	//18
	tmpTest.fName = DIXMAANI;
	tmpTest.Initializer = InitializeDIXMAAN;
	tmpTest.testName = "DIXMAANI function";
	tmpTest.size = 3000;
	testsVector.push_back(tmpTest);
	//19
	tmpTest.fName = DIXMAANJ;
	tmpTest.Initializer = InitializeDIXMAAN;
	tmpTest.testName = "DIXMAANJ function";
	tmpTest.size = 3000;
	testsVector.push_back(tmpTest);
	//20
	tmpTest.fName = DIXMAANK;
	tmpTest.Initializer = InitializeDIXMAAN;
	tmpTest.testName = "DIXMAANK function";
	tmpTest.size = 3000;
	testsVector.push_back(tmpTest);
	//21
	tmpTest.fName = DIXMAANL;
	tmpTest.Initializer = InitializeDIXMAAN;
	tmpTest.testName = "DIXMAANL function";
	tmpTest.size = 3000;
	testsVector.push_back(tmpTest);
	//22
	tmpTest.fName = DIXON3DQ;
	tmpTest.Initializer = InitialDIXON3DQ;
	tmpTest.testName = "DIXON3DQ function";
	tmpTest.size = 5000;
	testsVector.push_back(tmpTest);
	//23
	tmpTest.fName = DQDRTIC;
	tmpTest.Initializer = InitialDQDRTIC;
	tmpTest.testName = "DQDRTIC function";
	tmpTest.size = 10000;
	testsVector.push_back(tmpTest);
	//24
	tmpTest.fName = DQRTIC;
	tmpTest.Initializer = InitialDQRTIC;
	tmpTest.testName = "DQRTIC function";
	tmpTest.size = 5000;
	testsVector.push_back(tmpTest);
	//25
	tmpTest.fName = EDENSCH;
	tmpTest.Initializer = InitialEDENSCH;
	tmpTest.testName = "EDENSCH function";
	tmpTest.size = 2000;
	testsVector.push_back(tmpTest);
	//26
	tmpTest.fName = EG2;
	tmpTest.Initializer = InitialEG2;
	tmpTest.testName = "EG2 function";
	tmpTest.size = 5000;
	testsVector.push_back(tmpTest);
	//27
	tmpTest.fName = ENGVAL1;
	tmpTest.Initializer = InitialENGVAL1;
	tmpTest.testName = "ENGVAL1 function";
	tmpTest.size = 5000;
	testsVector.push_back(tmpTest);
	//28
	tmpTest.fName = EXTROSNB;
	tmpTest.Initializer = InitialEXTROSNB;
	tmpTest.testName = "EXTROSNB";
	tmpTest.size = 1000;
	testsVector.push_back(tmpTest);
	//29
	tmpTest.fName = FLETCHR;
	tmpTest.Initializer = InitialFLETCHR;
	tmpTest.testName = "FLETCHR function";
	tmpTest.size = 1000;
	testsVector.push_back(tmpTest);
	//30
	tmpTest.fName = FREUROTH;
	tmpTest.Initializer = InitialFREUROTH;
	tmpTest.testName = "FREUROTH";
	tmpTest.size = 5000;
	testsVector.push_back(tmpTest);
	//31
	tmpTest.fName = GENHUMPS;
	tmpTest.Initializer = InitialGENHUMPS;
	tmpTest.testName = "GENHUMPS";
	tmpTest.size = 5000;
	testsVector.push_back(tmpTest);
	//32
	tmpTest.fName = GENROSE;
	tmpTest.Initializer = InitialGENROSE;
	tmpTest.testName = "GENROSE";
	tmpTest.size = 1000;
	testsVector.push_back(tmpTest);
	//33
	tmpTest.fName = LIARWDH;
	tmpTest.Initializer = InitialLIARWDH;
	tmpTest.testName = "LIARWDH function";
	tmpTest.size = 5000;
	testsVector.push_back(tmpTest);
	//34
	tmpTest.fName = MOREBV;
	tmpTest.Initializer = InitialMOREBV;
	tmpTest.testName = "MOREBV";
	tmpTest.size = 10000;
	testsVector.push_back(tmpTest);
	//35
	tmpTest.fName = NONCVXU2;
	tmpTest.Initializer = InitialNONCVXU2;
	tmpTest.testName = "NONCVXU2";
	tmpTest.size = 1000;
	testsVector.push_back(tmpTest);
	//36
	tmpTest.fName = NONDIA;
	tmpTest.Initializer = InitialNONDIA;
	tmpTest.testName = "NONDIA function";
	tmpTest.size = 10000;
	testsVector.push_back(tmpTest);
	//37
	tmpTest.fName = NONDQUAR;
	tmpTest.Initializer = InitialNONDQUAR;
	tmpTest.testName = "NONDQUAR function";
	tmpTest.size = 10000;
	testsVector.push_back(tmpTest);
	//38
	tmpTest.fName = PENALTY1;
	tmpTest.Initializer = InitialPENALTY1;
	tmpTest.testName = "PENALTY1";
	tmpTest.size = 1000;
	testsVector.push_back(tmpTest);
	//39
	tmpTest.fName = PENALTY2;
	tmpTest.Initializer = InitialPENALTY2;
	tmpTest.testName = "PENALTY2";
	tmpTest.size = 150;
	testsVector.push_back(tmpTest);
	//40
	tmpTest.fName = POWER;
	tmpTest.Initializer = InitialPOWER;
	tmpTest.testName = "POWER function";
	tmpTest.size = 1000;
	testsVector.push_back(tmpTest);
	//41
	tmpTest.fName = SROSENBR;
	tmpTest.Initializer = InitialSROSENBR;
	tmpTest.testName = "SROSENBR function";
	tmpTest.size = 10000;
	testsVector.push_back(tmpTest);
	//42
	tmpTest.fName = TRIDIA;
	tmpTest.Initializer = InitialTRIDIA;
	tmpTest.testName = "TRIDIA function";
	tmpTest.size = 10000;
	testsVector.push_back(tmpTest);
	//43
	tmpTest.fName = WOODS;
	tmpTest.Initializer = InitialWOODS;
	tmpTest.testName = "WOODS";
	tmpTest.size = 10000;
	testsVector.push_back(tmpTest);
}

bool stopingByInfNorm()
{
	return (infNorm(g0, n) < EPS);
}

void runUCONTests()
{
	int repNumber(1);
	makeTestsVector();
	vector<testFunction>& v = testsVector;

	std::cout << "Enter test's index between  "
		<< 0 << "  and  " << v.size() - 1 << endl;
	int i;
	cin >> i;
	v[i].show();
	n = v[i].size;							//seting tsize
	u_min::constructData(1E-6);				//seting tolerance 
	vector<_int64> repetitions(repNumber);
	freopen("results.txt", "w", stdout);

	for (int j = 0; j < repNumber; j++)
	{
		v[i].Initializer(x0, n);
		lineSearch  lnSrch;
		auto st = chrono::high_resolution_clock::now();
		u_min::mhb(lnSrch, v[i].fName, stopingByInfNorm);
		auto diff = chrono::high_resolution_clock::now() - st;
		auto time = chrono::duration_cast<chrono::microseconds>(diff);
		repetitions[j] = time.count();
	}
	sort(repetitions.begin(), repetitions.end());
	std::cout << "Average time:  " << repetitions[repNumber / 2] << " microseconds" << endl;
	std::cout << endl;
	{
		v[i].show();
		printSolution(f0, x0, n);
	}
}