#pragma once  
#include <string>
#include"problems.h"
#include<algorithm>
using namespace std;

template<class T> 
class MHB
{
	problem<T>* ppr{};
	long int n{};
	long int FuncGradEvaluations;
	long int counter{};
	long int LineSearchCunter{};
	T EPS{};
	T alpha{};
	T f0;			//objective function's value at x0 - მიზნის ფუნქციის მნიშვნელობა x0-ში
	T f1;
	T* x0, * x1, * dir, * tmp;
	T* g1, * g0;

	//line search
	bool lineSearch() noexcept;

	//used in line search
	const T RO_Right = 10;
	const T RO_Left = 1e-3;
	T  t{1e-3};
	T  fi_alpha, fi_t;
public:
	MHB(const T tolerance, problem<T>* );
	~MHB();

	//infinity norm of the vector
	T infNorm(T*) const noexcept;

	//vector name is supplied as a string object
	void printVector(T*, std::string) const noexcept;
	long int getSize() const noexcept {	return n; }
	long int getRestarts() const noexcept { return LineSearchCunter; }
	long int getEvaluations() const noexcept { return FuncGradEvaluations; }
	T* getVariables() const noexcept {	return x0;}
	T* getGradient() const noexcept { return g0; }
	T getValue() const noexcept { return f0; }
	//solver
	bool solve() noexcept;
};

//---------------------------Implemetations------------------------------------
/*****************************************************************************/

template<class T>
MHB<T>::MHB(const T tolerance, problem<T>* problemAddress)
{
	ppr = problemAddress;
	n = ppr->getSize();
	EPS = tolerance;
	alpha = 0.001;
	dir = new T[n];
	x0 = new T[n];
	x1 = new T[n];
	g1 = new T[n];
	g0 = new T[n];
	tmp = new T[n];
}
template<class T>
MHB<T>::~MHB()
{
	delete[] dir;
	delete[] x0;
	delete[] x1;
	delete[] g0;
	delete[] g1;
	delete[] tmp;
}

template<class T>
T MHB<T>::infNorm(T* x) const noexcept
{
	T mx{};
	T tmp{};
	for (size_t i = 0; i < n; ++i)
	{
		tmp = x[i];
		if (mx < fabs(tmp))
			mx = fabs(tmp);
	}
	return mx;
}

template<class T>
void MHB<T>::printVector(T* vec, std::string s) const noexcept
{
	std::cout << "Objective function:  " << f0 << 
		",   Value-Grad evaluations:  " << getEvaluations() << 
		",   Restarts (line searches):  " << getRestarts() << std::endl << std::endl;
	long int quarter = n / 4;
	uint8_t indDigits{};
	if (n < 1001) indDigits = 3;
	else
		if(n<10001)indDigits = 4;
		else if (n < 100001)indDigits = 5;
	for (size_t i = 0; i < quarter; ++i)
	{
		std::cout << std::left << s<< std::left << '[' << std::left << std::setw(indDigits)<< i << std::left << "]="
			<< std::left << std::setw(21) << vec[i];
		std::cout << std::left << s << std::left << '[' << std::left << std::setw(indDigits) << i + quarter << std::left << "]="
			<< std::left << std::setw(21) << vec[i + quarter];
		std::cout << std::left << s << std::left << '[' << std::left << std::setw(indDigits) << i + 2 * quarter << std::left << "]="
			<< std::left << std::setw(21) << vec[i + 2 * quarter];
		std::cout << std::left << s << std::left << '[' << std::left << std::setw(indDigits) << i + 3* quarter << std::left << "]="
			<< std::left << std::setw(21) << vec[i + 3 * quarter];
		std::cout << std::endl;
	}
	for (int i = 4 * quarter; i < n; i++)
		std::cout << std::left << std::setw(16) << s << '[' << i << "]=" << vec[i];
	std::cout << std::endl;
}

template<class T>
bool MHB<T>::lineSearch() noexcept
{
	for (int i = 0; i < n; i++)
		x1[i] = x0[i] - t * g0[i];
	fi_t = ppr->valGrad(x1, g1);
	++FuncGradEvaluations;

	while (fi_t >= f0)
	{
		t *= RO_Left;
		if (t < std::numeric_limits<T>::epsilon())
			return false;
		for (int i = 0; i < n; i++)
			x1[i] = x0[i] - t * g0[i];
		fi_t = ppr->valGrad(x1, g1);
		++FuncGradEvaluations;
	}

	alpha = t;
	fi_alpha = fi_t;
	while (true)
	{
		t *= RO_Right;
		for (int i = 0; i < n; i++)
			tmp[i] = x0[i] - t * g0[i];
		fi_t = ppr->valGrad(tmp, dir);
		if (fi_t >= fi_alpha)
			break;
		alpha = t;
		fi_alpha = fi_t;
		swap(x1, tmp);
		swap(dir, g1);
	}

	f1 = fi_alpha;
	return true;
}

template<class T>
bool MHB<T>::solve() noexcept
{
	f0 = f1 = ppr->valGrad(x0, g0); 
	++FuncGradEvaluations;
	if (infNorm(g0)<EPS) 	return true;
	while (true)
	{
		if (f0 <= f1 || counter == 100000)
		{
			++LineSearchCunter;
			if (lineSearch() == false) return false;
			counter = 0;
		}
		///	cout << "f0=" << f0<<   "   " << infNorm(g0)  <<  endl;   ///////

		for (int i = 0; i < n; ++i)
			dir[i] = x1[i] - x0[i];
		swap(x0, x1);
		swap(g0, g1);

		if (infNorm(g0) < EPS) 	return true;

		f0 = f1;
		for (int i = 0; i < n; ++i)
			x1[i] = x0[i] + dir[i] - alpha * g0[i];
		f1 = ppr->valGrad(x1, g1);
		++FuncGradEvaluations;
		++counter;
	}
}