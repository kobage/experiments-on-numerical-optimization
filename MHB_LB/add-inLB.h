#pragma once
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

int lbNumbers;				// lbNumbers = number of lower constrained variables
double* lBounds;			// Array of all lower bounds
int* lbVarsInds;			// Indexes of lower bound constrained variables
int* freeVarsInds;			// Indexes of free variables
double* y;					// Decision variable after transformation - we are changing the notion of
							// The old variable y[i] = lb[i] + x0[i]*x0[i]
double **A;
double Dev;

//>>>>>> Component 1: Change of Variables,  their representation and container of tests <<<<< 

//Pointers to functions, visible from functions
double(*substitutionName)(double x);
double(*inversesubstitutionName)(double x);
double(*substitutionPrimeName)(double x);

/*===== 0 =================== Quadratic substitution =================================*/
double quadraticSubstitution(double v)
{
	return v*v;
}
double quadraticSubstitutionPrime(double v)
{
	return 2 * v;
}
double inverseQuadraticSubstitution(double v)
{
	return sqrt(v);
}

/*===== 1 =================== Quartic substitution =================================*/
double quarticSubstitution(double v)
{
	return v*v*v*v;
}
double quarticSubstitutionPrime(double v)
{
	return 4 * v*v*v;
}
double inverseQuarticSubstitution(double v)
{
	//return sqrt(sqrt(v));
	return sqrt(sqrt(v));
}

/*===== 2 =================== Exponential substitution =================================*/
double exponentialSubstitution(double v)
{
	return exp(v);
}
double exponentialcSubstitutionPrime(double v)
{
	return exp(v);
}
double inverseExponentialSubstitution(double v)
{
	return log(v);
}

/*===== 3 ================ Sinusoidal substitution, only for symmetric game ==================*/
double sinusoidalSubstitution(double v)
{
	return 1 + sin(v);
}
double sinusoidalcSubstitutionPrime(double v)
{
	return cos(v);		
}
double inverseSinusoidalSubstitution(double v)
{
	return asin(v - 1);		
}
//--------------------------- End of substitutions ----------------------------------------
//structure for a substitution
struct lb_substiturtion
{
	double(*substitution) (double);
	double(*inverseSubstitution)(double);
	double(*substitutionPrime)(double);
	string substitutionName;
	void show(void)
	{
		cout << "Substitution:  " << substitutionName << endl;
	}
};
//Container of substitutions
vector<lb_substiturtion> lb_substiturtionsVector;
void makeLb_substiturtionsVector(void)
{
	lb_substiturtion tmpSubstitution;
	//0
	tmpSubstitution.substitution = quadraticSubstitution;
	tmpSubstitution.inverseSubstitution = inverseQuadraticSubstitution;
	tmpSubstitution.substitutionPrime = quadraticSubstitutionPrime;
	tmpSubstitution.substitutionName = "Quadratic substitution";
	lb_substiturtionsVector.push_back(tmpSubstitution);
	//1
	tmpSubstitution.substitution = quarticSubstitution;
	tmpSubstitution.inverseSubstitution = inverseQuarticSubstitution;
	tmpSubstitution.substitutionPrime = quarticSubstitutionPrime;
	tmpSubstitution.substitutionName = "Quartic substitution";
	lb_substiturtionsVector.push_back(tmpSubstitution);
	//2
	tmpSubstitution.substitution = exponentialSubstitution;
	tmpSubstitution.inverseSubstitution = inverseExponentialSubstitution;
	tmpSubstitution.substitutionPrime = exponentialcSubstitutionPrime;
	tmpSubstitution.substitutionName = "Exponential substitution";
	lb_substiturtionsVector.push_back(tmpSubstitution);
	//3
	tmpSubstitution.substitution = sinusoidalSubstitution;
	tmpSubstitution.inverseSubstitution = inverseSinusoidalSubstitution;
	tmpSubstitution.substitutionPrime = sinusoidalcSubstitutionPrime;
	tmpSubstitution.substitutionName = "Sinusoidal substitution, only for symmetric game";
	lb_substiturtionsVector.push_back(tmpSubstitution);
}

//>>>>>> Component 2: LB-problems, their representation and container of tests <<<<<<

/*===== 0 =================== hatflda Function ================ 4 =================*/
double hatflda
(
	double *x,
	double *g,
	const int k
)
{
	for (int i = 0; i < k; i++)
		g[i] = 0;

	double group(x[0] - 1);
	double fx = group*group;
	g[0] += 2.0*group;

	for (int i = 1; i < k; i++)
	{
		group = x[i - 1] - sqrt(x[i]);
		fx += group*group;
		g[i - 1] += 2.0*group;
		g[i] -= group / sqrt(x[i]);
	}
	return fx;
}

void InitializeHatflda
(
	double *x,
	const int k
)
{
	for (int i = 0; i < k; i++)
		x[i] = 0.1;

	y = (double *)malloc(k * sizeof(double));
	lBounds = (double *)malloc(k * sizeof(double));
	for (int i = 0; i < k; i++)
		lBounds[i] = 0.0000001;
	for (int i = 0; i < k; i++)
	{
		x[i] = inversesubstitutionName(x[i] - lBounds[i]);
	}
	lbNumbers = k;
	lbVarsInds = (int *)malloc(lbNumbers * sizeof(int));
	for (int i = 0; i < lbNumbers; i++)
		lbVarsInds[i] = i;
}

void recoveringAttributesHatflda
(
	double *x,
	const int k
)
{
	free(y);
	for (int i = 0; i < k; i++)
	{
		x[i] = substitutionName(x[i]) + lBounds[i];
	}
	free(lBounds);
	free(lbVarsInds);
}

double hatfldaComposite
(
	double *x,
	double *g,
	const int k
)
{
	for (int i = 0; i < k; i++)
		y[i] = lBounds[i] + substitutionName(x[i]);

	double fx = hatflda(y, g, k);
	for (int i = 1; i < k; i++)
		g[i] *= substitutionPrimeName(x[i]);
	return fx;
}

/*===== 1 ===================    hs001    Function ================ 2  ===================*/
double hs001
(
	double *x,
	double *g,
	const int k
)
{
	double group = x[1] - x[0] * x[0];
	double fx = 100 * group*group + (1 - x[0])*(1 - x[0]);
	g[0] = -400 * group* x[0] - 2 * (1 - x[0]);
	g[1] = 200 * group;
	return fx;
}

void InitializeHs001
(
	double *x,
	const int k
)
{
	x[0] = -2.;
	x[1] = 1.;

	y = (double *)malloc(k * sizeof(double));
	lbNumbers = 1;
	lBounds = (double *)malloc(lbNumbers * sizeof(double));
	lBounds[0] = -1.5;
	lbVarsInds = (int *)malloc(lbNumbers * sizeof(int));
	freeVarsInds = (int *)malloc((n - lbNumbers) * sizeof(int));
	lbVarsInds[0] = 1;
	freeVarsInds[0] = 0;
	int j;
	for (int i = 0; i < lbNumbers; i++)
	{
		j = lbVarsInds[i];
		x[j] = inversesubstitutionName(x[j] - lBounds[i]);
	}
}

void recoveringAttributesHs001
(
	double *x,
	const int k
)
{
	free(y);
	int j;
	for (int i = 0; i < lbNumbers; i++)
	{
		j = lbVarsInds[i];
		x[j] = substitutionName(x[j]) + lBounds[i];
	}
	free(lBounds);
	free(lbVarsInds);
	free(freeVarsInds);
}

double hs001Composite
(
	double *x,
	double *g,
	const int k
)
{
	y[0] = x[0];
	y[1] = lBounds[0] + substitutionName(x[1]);
	double fx = hs001(y, g, k);
	g[1] *= substitutionPrimeName(x[1]);
	return fx;
}

/*===== 2 ===================   hs002     Function ================ 2 ===================*/
void InitializeHs002
(
	double *x,
	const int k
)
{
	x[0] = -2.;
	x[1] = 1.;

	y = (double *)malloc(k * sizeof(double));
	lbNumbers = 1;
	lBounds = (double *)malloc(lbNumbers * sizeof(double));
	lBounds[0] = 1.5;
	lbVarsInds = (int *)malloc(lbNumbers * sizeof(int));
	freeVarsInds = (int *)malloc((n - lbNumbers) * sizeof(int));
	lbVarsInds[0] = 1;
	freeVarsInds[0] = 0;
	int j;
	for (int i = 0; i < lbNumbers; i++)
	{
		j = lbVarsInds[i];
		x[j] = inversesubstitutionName(fabs(x[j] - lBounds[i]));		//???
	}
}

void recoveringAttributesHs002
(
	double *x,
	const int k
)
{
	free(y);
	int j;
	for (int i = 0; i < lbNumbers; i++)
	{
		j = lbVarsInds[i];
		x[j] = substitutionName(x[j]) + lBounds[i];
	}
	free(lBounds);
	free(lbVarsInds);
	free(freeVarsInds);
}
/*===== 3 ===================   hs003    Function ================ 2 ===================*/
double hs003
(
	double *x,
	double *g,
	const int k
)
{
	double fx = 0.00001 * (x[1] - x[0])*(x[1] - x[0]) + x[1];
	g[0] = -0.00002 * (x[1] - x[0]);
	g[1] = 0.00002 * (x[1] - x[0]) + 1;
	return fx;
}

void InitializeHs003
(
	double *x,
	const int k
)
{
	x[0] = 10.;
	x[1] = 1.;

	y = (double *)malloc(k * sizeof(double));
	lbNumbers = 1;
	lBounds = (double *)malloc(lbNumbers * sizeof(double));
	lBounds[0] = 0.;
	lbVarsInds = (int *)malloc(lbNumbers * sizeof(int));
	lbVarsInds[0] = 1;
	freeVarsInds = (int *)malloc((n - lbNumbers) * sizeof(int));
	freeVarsInds[0] = 0;
	int j;
	for (int i = 0; i < lbNumbers; i++)
	{
		j = lbVarsInds[i];
		x[j] = inversesubstitutionName(x[j] - lBounds[i]);
	}
}

void recoveringAttributesHs003
(
	double *x,
	const int k
)
{
	free(y);
	int j;
	for (int i = 0; i < lbNumbers; i++)
	{
		j = lbVarsInds[i];
		x[j] = substitutionName(x[j]) + lBounds[i];
	}
	free(lBounds);
	free(lbVarsInds);
	free(freeVarsInds);
}

double hs003Composite
(
	double *x,
	double *g,
	const int k
)
{
	y[0] = x[0];
	y[1] = lBounds[0] + substitutionName(x[1]);
	double fx = hs003(y, g, k);
	g[1] *= substitutionPrimeName(x[1]);
	return fx;
}

/*===== 4 ===================   hs004    Function ================ 2 ===================*/
double hs004
(
	double *x,
	double *g,
	const int k
)
{
	double fx = (1 + x[0])*(1 + x[0])*(1 + x[0]) + x[1];
	g[0] = 3 * (1 + x[0]);
	g[1] = 1;
	return fx;
}

void InitializeHs004
(
	double *x,
	const int k
)
{
	x[0] = 1.125;
	x[1] = 0.125;

	y = (double *)malloc(k * sizeof(double));
	lbNumbers = 2;
	lBounds = (double *)malloc(lbNumbers * sizeof(double));
	lBounds[0] = 1.;
	lBounds[1] = 0.;
	lbVarsInds = (int *)malloc(lbNumbers * sizeof(int));
	lbVarsInds[0] = 0;
	lbVarsInds[1] = 1;
	int j;
	for (int i = 0; i < lbNumbers; i++)
	{
		j = lbVarsInds[i];
		x[j] = inversesubstitutionName(x[j] - lBounds[i]);
	}
}

void recoveringAttributesHs004
(
	double *x,
	const int k
)
{
	free(y);
	int j;
	for (int i = 0; i < lbNumbers; i++)
	{
		j = lbVarsInds[i];
		x[j] = substitutionName(x[j]) + lBounds[i];
	}
	free(lBounds);
	free(lbVarsInds);
}

double hs004Composite
(
	double *x,
	double *g,
	const int k
)
{
	y[0] = lBounds[0] + substitutionName(x[0]);
	y[1] = lBounds[1] + substitutionName(x[1]);
	double fx = hs004(y, g, k);
	g[0] *= substitutionPrimeName(x[0]);
	g[1] *= substitutionPrimeName(x[1]);
	return fx;
}

/*===== 5 ==================   palmer1a   Function =============== 6 ================*/
double *X, *Y;
double palmer1a
(
	double *x,
	double *g,
	const int k
)
{
	double fx = 0.;
	for (int i = 0; i < 6; ++i) g[i] = 0.;
	double sq;
	double group;
	for (int m = 0; m < 35; ++m)
	{
		sq = X[m] * X[m];
		group = Y[m] - x[0] - x[1] * sq - x[2] * sq*sq - x[3] * sq*sq*sq - x[4] / (x[5] + sq);

		fx += group*group;

		g[0] -= 2 * group;
		g[1] -= 2 * group*sq;
		g[2] -= 2 * group*sq*sq;
		g[3] -= 2 * group*sq*sq*sq;
		g[4] -= (2 * group) / (x[5] + sq);
		g[5] += (2 * group*x[4]) / ((x[5] + sq)*(x[5] + sq));
	}
	return fx;
}

void Initializepalmer1a
(
	double *x,
	const int k
)
{
	x[0] = 1.0;		//A0
	x[1] = 1.0;		//A2
	x[2] = 1.0;		//A4
	x[3] = 1.0;		//A6
	x[4] = 1.0;		//B>=0
	x[5] = 1.0;		//C>=0

	X = (double *)malloc(35 * sizeof(double));
	Y = (double *)malloc(35 * sizeof(double));
	double tmpX[] = { -1.788963 ,-1.745329, -1.658063, -1.570796,-1.483530,-1.396263,-1.308997,
		-1.218612,-1.134464,-1.047198, -0.872665,-0.698132, -0.523599 ,-0.349066,-0.174533, 0.0000000
		, 1.788963, 1.745329, 1.658063, 1.570796, 1.483530, 1.396263, 1.308997, 1.218612
		, 1.134464, 1.047198, 0.872665, 0.698132, 0.523599, 0.349066, 0.174533,-1.8762289, -1.8325957
		, 1.8762289, 1.8325957 };
	for (int i = 0; i < 35; ++i)	X[i] = tmpX[i];
	double tmpY[] = { 78.596218,65.77963,43.96947, 27.038816, 14.6126, 6.2614, 1.538330, 0.000000
		,1.188045, 4.6841, 16.9321, 33.6988, 52.3664, 70.1630, 83.4221, 88.3995, 78.596218, 65.77963
		, 43.96947, 27.038816, 14.6126, 6.2614, 1.538330, 0.000000, 1.188045, 4.6841, 16.9321, 33.6988
		, 52.3664, 70.1630, 83.4221, 108.18086, 92.733676, 108.18086, 92.733676 };
	for (int i = 0; i < 35; ++i)	Y[i] = tmpY[i];


	y = (double *)malloc(k * sizeof(double));
	lbNumbers = 2;
	lBounds = (double *)malloc(lbNumbers * sizeof(double));
	lBounds[0] = 0.;
	lBounds[1] = 0.;
	lbVarsInds = (int *)malloc(lbNumbers * sizeof(int));
	lbVarsInds[0] = 4;
	lbVarsInds[1] = 5;
	int j;
	for (int i = 0; i < lbNumbers; i++)
	{
		j = lbVarsInds[i];
		x[j] = inversesubstitutionName(x[j] - lBounds[i]);
	}
	freeVarsInds = (int *)malloc((n - lbNumbers) * sizeof(int));
	freeVarsInds[0] = 0;
	freeVarsInds[1] = 1;
	freeVarsInds[2] = 2;
	freeVarsInds[3] = 3;
}

void recoveringAttributespalmer1a
(
	double *x,
	const int k
)
{
	free(y);
	int j;
	for (int i = 0; i < lbNumbers; i++)
	{
		j = lbVarsInds[i];
		x[j] = substitutionName(x[j]) + lBounds[i];
	}
	free(lBounds);
	free(lbVarsInds);
	free(freeVarsInds);
	free(X);
	free(Y);
}

double palmer1aComposite
(
	double *x,
	double *g,
	const int k
)
{
	for (int i = 0; i < 4; ++i) y[i] = x[i];
	y[4] = lBounds[0] + substitutionName(x[4]);
	y[5] = lBounds[1] + substitutionName(x[5]);

	double fx = palmer1a(y, g, k);
	g[4] *= substitutionPrimeName(x[4]);
	g[5] *= substitutionPrimeName(x[5]);
	return fx;
}

/*===== 6 ==================   palmer1   Function =============== 4 ================*/
double palmer1
(
	double *x,
	double *g,
	const int k
)
{
	double fx = 0.;
	for (int i = 0; i < 4; ++i) g[i] = 0.;
	double sq;
	double group;
	for (int m = 0; m < 31; ++m)
	{
		sq = X[m] * X[m];
		group = Y[m] - (x[0] * sq + x[1] / (x[2] + sq / x[3]));

		fx += group*group;

		g[0] -= 2 * group*sq;
		g[1] -= 2 * group / (x[2] + sq / x[3]);
		g[2] += 2 * group*x[1] / ((x[2] + sq / x[3])*(x[2] + sq / x[3]));
		g[3] -= 2 * group*x[1] * sq / ((x[2] + sq / x[3])*(x[2] + sq / x[3])* x[3] * x[3]);
	}
	return fx;
}

void Initializepalmer1
(
	double *x,
	const int k
)
{
	x[0] = 1.0;		//A
	x[1] = 1.0;		//B>=0
	x[2] = 1.0;		//C>=0
	x[3] = 1.0;		//D>=0

	int M = 31;
	X = (double *)malloc(M * sizeof(double));
	Y = (double *)malloc(M * sizeof(double));
	double tmpX[] = { -1.788963 ,-1.745329, -1.658063, -1.570796,-1.483530,-1.396263,-1.308997,
		-1.218612,-1.134464,-1.047198, -0.872665,-0.698132, -0.523599 ,-0.349066,-0.174533, 0.0000000
		, 1.788963, 1.745329, 1.658063, 1.570796, 1.483530, 1.396263, 1.308997, 1.218612
		, 1.134464, 1.047198, 0.872665, 0.698132, 0.523599, 0.349066, 0.174533 };
	for (int i = 0; i < M; ++i)	X[i] = tmpX[i];
	double tmpY[] = { 78.596218,65.77963,43.96947, 27.038816, 14.6126, 6.2614, 1.538330, 0.000000
		,1.188045, 4.6841, 16.9321, 33.6988, 52.3664, 70.1630, 83.4221, 88.3995, 78.596218, 65.77963
		, 43.96947, 27.038816, 14.6126, 6.2614, 1.538330, 0.000000, 1.188045, 4.6841, 16.9321, 33.6988
		, 52.3664, 70.1630, 83.4221 };
	for (int i = 0; i < M; ++i)	Y[i] = tmpY[i];

	y = (double *)malloc(k * sizeof(double));
	lbNumbers = 3;
	lBounds = (double *)malloc(lbNumbers * sizeof(double));
	lBounds[0] = 0.00001;
	lBounds[1] = 0.00001;
	lBounds[2] = 0.00001;
	lbVarsInds = (int *)malloc(lbNumbers * sizeof(int));
	lbVarsInds[0] = 1;
	lbVarsInds[1] = 2;
	lbVarsInds[2] = 3;
	int j;
	for (int i = 0; i < lbNumbers; i++)
	{
		j = lbVarsInds[i];
		x[j] = inversesubstitutionName(x[j] - lBounds[i]);
	}
	freeVarsInds = (int *)malloc((k - lbNumbers) * sizeof(int));
	freeVarsInds[0] = 0;
}

void recoveringAttributespalmer1
(
	double *x,
	const int k
)
{
	free(y);
	int j;
	for (int i = 0; i < lbNumbers; i++)
	{
		j = lbVarsInds[i];
		x[j] = substitutionName(x[j]) + lBounds[i];
	}
	free(lBounds);
	free(lbVarsInds);
	free(freeVarsInds);
}

double palmer1Composite
(
	double *x,
	double *g,
	const int k
)
{
	y[0] = x[0];
	y[1] = lBounds[0] + substitutionName(x[1]);
	y[2] = lBounds[1] + substitutionName(x[2]);
	y[3] = lBounds[2] + substitutionName(x[3]);
	double fx = palmer1(y, g, k);
	g[1] *= substitutionPrimeName(x[1]);
	g[2] *= substitutionPrimeName(x[2]);
	g[3] *= substitutionPrimeName(x[3]);
	return fx;
}
/*===== 7 ==================   pspdoc   Function =============== 4 ================*/
double pspdoc
(
	double *x,
	double *g,
	const int k
)
{
	double group = x[0] * x[0] + (x[1] - x[2])*(x[1] - x[2]) + 1;
	double fx = sqrt(group);

	g[0] = x[0] / sqrt(group);
	g[1] = (x[1] - x[2]) / sqrt(group);
	g[2] = -(x[1] - x[2]) / sqrt(group);

	group = x[1] * x[1] + (x[2] - x[3])*(x[2] - x[3]) + 1;
	fx += sqrt(group);

	g[1] += x[1] / sqrt(group);
	g[2] += (x[2] - x[3]) / sqrt(group);
	g[3] = -(x[2] - x[3]) / sqrt(group);

	return fx;
}

void Initializepspdoc
(
	double *x,
	const int k
)
{
	x[0] = 3.0;			//>=1
	x[1] = 3.0;
	x[2] = 3.0;
	x[3] = 3.0;

	y = (double *)malloc(k * sizeof(double));
	lbNumbers = 1;
	lBounds = (double *)malloc(lbNumbers * sizeof(double));
	lBounds[0] = 1;

	lbVarsInds = (int *)malloc(lbNumbers * sizeof(int));
	lbVarsInds[0] = 0;

	int j;
	for (int i = 0; i < lbNumbers; i++)
	{
		j = lbVarsInds[i];
		x[j] = inversesubstitutionName(x[j] - lBounds[i]);
	}
	freeVarsInds = (int *)malloc((n - lbNumbers) * sizeof(int));
	freeVarsInds[0] = 1;
	freeVarsInds[1] = 2;
	freeVarsInds[2] = 3;
}

void recoveringAttributespspdoc
(
	double *x,
	const int k
)
{
	free(y);
	int j;
	for (int i = 0; i < lbNumbers; i++)
	{
		j = lbVarsInds[i];
		x[j] = substitutionName(x[j]) + lBounds[i];
	}
	x[0] *= -1;
	free(lBounds);
	free(lbVarsInds);
	free(freeVarsInds);
}

double pspdocComposite
(
	double *x,
	double *g,
	const int k
)
{
	y[0] = lBounds[0] + substitutionName(x[0]);
	y[1] = x[1];
	y[2] = x[2];
	y[3] = x[3];

	double fx = pspdoc(y, g, k);
	g[0] *= substitutionPrimeName(x[0]);
	return fx;
}
/*===== 8 ===  Symmetric game with objective function, cubic with respect to original variables ========*/
//The following Tests are using A and Dev
double sgValGrad
(
	double *x,
	double *g,
	const int k
)
{
	double fx = 0.;
	Dev = 0.;
	int i, j;
	double prod, qu;
	double sum = 0.;
	double* p;
	for (i = 0; i < k; i++)
	{
		g[i] = 0;
		sum += x[i];
	}
	for (i = 0; i < k; i++)
	{
		p = A[i];
		prod = vecProd(p, x, k);
		if (prod > 0.)
		{
			if (Dev < prod) Dev = prod;
			qu = prod*prod;
			fx += prod * qu;
			for (j = 0; j < k; j++)
				g[j] += qu*p[j];
		}
	}
	if (Dev < fabs(sum - 1)) Dev = fabs(sum - 1);
	int sigh = (sum < 1) ? (-1) : 1;
	prod = (sum - 1)*(sum - 1)*sigh;
	fx += prod*(sum - 1);
	for (i = 0; i < k; i++)
		g[i] += prod;
	return fx / 3;
}

void Initialize_sg
(
	double *x,
	const int k
)
{
	y = (double *)malloc(k * sizeof(double));
	A = (double **)malloc(sizeof(double *)*k);
	for (int i = 0; i<k; i++)
		A[i] = (double *)malloc(sizeof(double)*k);

	fillRandomMatrix(A, k);

	for (int i = 0; i < k; i++)
		x[i] = 1.0 / k;

	lbNumbers = k;
	lBounds = (double *)malloc(lbNumbers * sizeof(double));
	lbVarsInds = (int *)malloc(lbNumbers * sizeof(int));
	for (int i = 0; i < k; i++)
	{
		lBounds[i] = 0.;
		lbVarsInds[i] = i;
	}
	for (int i = 0; i < k; i++)
		x[i] = inversesubstitutionName(x[i]);
}

void recoveringAttributes_sg
(
	double *x,
	const int k
)
{
	free(y);
	int i;
	for (i = 0; i < k; i++)
	{
		x[i] = substitutionName(x[i]);
		free(A[i]);
	}
	free(lBounds);
	free(lbVarsInds);
	free(A);
}

double sgComposite
(
	double *x,
	double *g,
	const int k
)
{
	for (int i = 0; i < k; i++)
		y[i] = substitutionName(x[i]);

	double fx = sgValGrad(y, g, k);
	for (int i = 0; i < k; i++)
		g[i] *= substitutionPrimeName(x[i]);
	return fx;
}
/*===== 9 ===  Symmetric game with objective function, quadratic with respect to original variables ========*/
double sgValGradQuadratic
(
	double *x,
	double *g,
	const int k
)
{
	double fx = 0.;
	Dev = 0.;
	int i, j;
	double prod;
	double sum = 0.;
	double* p;
	for (i = 0; i < k; i++)
	{
		g[i] = 0;
		sum += x[i];
	}
	for (i = 0; i < k; i++)
	{
		p = A[i];
		prod = vecProd(p, x, k);
		if (prod > 0.)
		{
			if (Dev < prod) Dev = prod;
			fx += prod*prod;
			for (j = 0; j < k; j++)
				g[j] += prod*p[j];
		}
	}
	if (Dev < fabs(sum - 1)) Dev = fabs(sum - 1);
	fx += (sum - 1)*(sum - 1);
	for (i = 0; i < k; i++)
		g[i] += (sum - 1);
	return fx / 2;
}
void recoveringAttributes_sgQuadratic
(
	double *x,
	const int k
)
{
	free(y);
	int i;
	for (i = 0; i < k; i++)
	{
		x[i] = substitutionName(x[i]);
		free(A[i]);
	}
	free(lBounds);
	free(lbVarsInds);
	free(A);
}

double sgQuadraticComposite
(
	double *x,
	double *g,
	const int k
)
{
	for (int i = 0; i < k; i++)
	{
		y[i] = substitutionName(x[i]);
	}
	double fx = sgValGrad(y, g, k);
	for (int i = 0; i < k; i++)
		g[i] *= substitutionPrimeName(x[i]);
	return fx;
}

//Container of tests
struct lb_testFunction
{
	double(*fName)				(double *, double*, const int);
	double(*fNameComposite)		(double*, double*, const int);
	void(*Initializer)			(double*, const int);
	void(*recoveringtAttributes)(double *, const int);

	int size;
	string testName;
	void show(void)
	{
		cout << "Name:  " << testName << endl;
		cout << "Size:  " << size << endl;
	}
};
vector<lb_testFunction> lb_testsVector;
void makeLbTestsVector(void)
{
	lb_testFunction tmpTest;
	//0
	tmpTest.fName = hatflda;
	tmpTest.fNameComposite = hatfldaComposite;
	tmpTest.Initializer = InitializeHatflda;
	tmpTest.testName = "hatflda function";
	tmpTest.recoveringtAttributes = recoveringAttributesHatflda;
	tmpTest.size = 4;
	lb_testsVector.push_back(tmpTest);
	//1
	tmpTest.fName = hs001;
	tmpTest.fNameComposite = hs001Composite;
	tmpTest.Initializer = InitializeHs001;
	tmpTest.testName = "hs001 function";
	tmpTest.recoveringtAttributes = recoveringAttributesHs001;
	tmpTest.size = 2;
	lb_testsVector.push_back(tmpTest);
	//2
	tmpTest.fName = hs001;
	tmpTest.fNameComposite = hs001Composite;
	tmpTest.Initializer = InitializeHs002;
	tmpTest.testName = "hs002 function";
	tmpTest.recoveringtAttributes = recoveringAttributesHs002;
	tmpTest.size = 2;
	lb_testsVector.push_back(tmpTest);
	//3
	tmpTest.fName = hs003;
	tmpTest.fNameComposite = hs003Composite;
	tmpTest.Initializer = InitializeHs003;
	tmpTest.testName = "hs003 function";
	tmpTest.recoveringtAttributes = recoveringAttributesHs003;
	tmpTest.size = 2;
	lb_testsVector.push_back(tmpTest);
	//4
	tmpTest.fName = hs004;
	tmpTest.fNameComposite = hs004Composite;
	tmpTest.Initializer = InitializeHs004;
	tmpTest.testName = "hs004 function";
	tmpTest.recoveringtAttributes = recoveringAttributesHs004;
	tmpTest.size = 2;
	lb_testsVector.push_back(tmpTest);
	//5
	tmpTest.fName = palmer1a;
	tmpTest.fNameComposite = palmer1aComposite;
	tmpTest.Initializer = Initializepalmer1a;
	tmpTest.testName = "palmer1a function";
	tmpTest.recoveringtAttributes = recoveringAttributespalmer1a;
	tmpTest.size = 6;
	lb_testsVector.push_back(tmpTest);
	//6
	tmpTest.fName = palmer1;
	tmpTest.fNameComposite = palmer1Composite;
	tmpTest.Initializer = Initializepalmer1;
	tmpTest.testName = "palmer1 function";
	tmpTest.recoveringtAttributes = recoveringAttributespalmer1;
	tmpTest.size = 4;
	lb_testsVector.push_back(tmpTest);
	//7
	tmpTest.fName = pspdoc;
	tmpTest.fNameComposite = pspdocComposite;
	tmpTest.Initializer = Initializepspdoc;
	tmpTest.testName = "pspdoc function";
	tmpTest.recoveringtAttributes = recoveringAttributespspdoc;
	tmpTest.size = 4;
	lb_testsVector.push_back(tmpTest);
	//8
	tmpTest.fName = sgValGrad;
	tmpTest.fNameComposite = sgComposite;
	tmpTest.Initializer = Initialize_sg;
	tmpTest.testName = "Symmetric game with objective function, cubic with respect to original variables";
	tmpTest.recoveringtAttributes = recoveringAttributes_sg;
	tmpTest.size = 200;
	lb_testsVector.push_back(tmpTest);
	//9
	tmpTest.fName = sgValGradQuadratic;
	tmpTest.fNameComposite = sgQuadraticComposite;
	tmpTest.Initializer = Initialize_sg;
	tmpTest.testName = "Symmetric game with objective function, quadratic with respect to original variables";
	tmpTest.recoveringtAttributes = recoveringAttributes_sgQuadratic;
	tmpTest.size = 400;
	lb_testsVector.push_back(tmpTest);
}

//>>>>>>>>>>>>>>   Component 3:   Stopping conditions  and program-driver  <<<<<<<<<<<<<

//------------------    Specific stopping condition (feasibility)    -------------------
bool stoppingByDeviation()
{
	return (Dev < EPS);
}
//-----------------    Stopping condition with gradient projection     -----------------
bool stoppingByProjGradInfNorm_LB()
{
	int j;
	for (int i = 0; i < n; i++)		g1[i] = 0.;
	for (int i = 0; i < lbNumbers; i++)
	{
		j = lbVarsInds[i];
		y[j] = substitutionName(x0[j]) + lBounds[i];
		g1[j] = g0[j] / substitutionPrimeName(x0[j]);
	}
	for (int i = 0; i < (n - lbNumbers); i++)
	{
		j = freeVarsInds[i];
		g1[j] = g0[j];
		y[j] = x0[j];
	}

	for (int i = 0; i < lbNumbers; i++)
	{
		j = lbVarsInds[i];
		if (y[j] - g1[j] < lBounds[j])
			g1[j] = lBounds[i] - y[j];
	}
	return (infNorm(g1, n) < EPS);
}

//Driver:
void runTests_LB()
{
	int repNumber(1);
	makeLbTestsVector();
	makeLb_substiturtionsVector();

	cout << "Available tests:" << endl;
	for (int i = 0; i < lb_testsVector.size(); ++i)
		cout << i << ".  " << lb_testsVector[i].testName << endl;
	std::cout << "Enter test's index between  "
		<< 0 << "  and  " << lb_testsVector.size() - 1 << endl;
	int i;
	cin >> i;

	cout << "Available transformations of variables::" << endl;
	for (int i = 0; i < lb_substiturtionsVector.size(); ++i)
		cout << i << ".  " << lb_substiturtionsVector[i].substitutionName << endl;
	std::cout << "Enter substitution's index between  "
		<< 0 << "  and  " << lb_substiturtionsVector.size() - 1 << endl;
	int j;
	cin >> j;

	lb_testsVector[i].show();
	lb_substiturtionsVector[j].show();

	substitutionName = lb_substiturtionsVector[j].substitution;
	substitutionPrimeName = lb_substiturtionsVector[j].substitutionPrime;
	inversesubstitutionName = lb_substiturtionsVector[j].inverseSubstitution;

	//-----------Seting tolerance and size-----------------//
	n = lb_testsVector[i].size;
	u_min::constructData(1E-6);
	vector<_int64> repetitions(repNumber);
	freopen("results.txt", "w", stdout);

	//Create expensive objects	
	auto st = chrono::high_resolution_clock::now();
	auto diff = chrono::high_resolution_clock::now() - st;
	auto time = chrono::duration_cast<chrono::microseconds>(diff);

	for (int j = 0; j<repNumber; j++)
	{
		st = chrono::high_resolution_clock::now();
		lineSearch  lnSrch;

		lb_testsVector[i].Initializer(x0, n);
		
		if (lb_testsVector[i].fNameComposite == sgComposite || lb_testsVector[i].fNameComposite == sgQuadraticComposite)
			u_min::MNAG(lnSrch, lb_testsVector[i].fNameComposite, stoppingByDeviation);
		//	u_min::mhb(lnSrch, lb_testsVector[i].fNameComposite, stoppingByDeviation);
		else
			u_min::MNAG(lnSrch, lb_testsVector[i].fNameComposite, stoppingByProjGradInfNorm_LB);
		//	u_min::mhb(lnSrch, lb_testsVector[i].fNameComposite, stoppingByProjGradInfNorm_LB);

		lb_testsVector[i].recoveringtAttributes(x0, n);

		diff = chrono::high_resolution_clock::now() - st;
		time = chrono::duration_cast<chrono::microseconds>(diff);
		repetitions[j] = time.count();
	}
	sort(repetitions.begin(), repetitions.end());
	std::cout << "Average time:  " << repetitions[repNumber / 2] << " microseconds" << endl;
	std::cout << endl;
	{
		lb_testsVector[i].show();
		lb_substiturtionsVector[j].show();
		printSolution(f0, x0, n);
	}
	u_min::destructData();
}