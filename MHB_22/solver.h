#pragma once  
#include <string>
#include"problems.h"
using namespace std;

template<class T> 
class MHB
{
	problem<T>* ppr{};
	long int n{};
	long int counter{};
	T EPS{};
	T alpha{};
	T f0;			//objective function's value at x0 - მიზნის ფუნქციის მნიშვნელობა x0-ში
	T f1;
	T* x0, * x1, * dir, * tmp;
	T* g1, * g0;

	//infinity norm of the vector
	T infNorm(T*) const noexcept;
	
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

	//vector name is supplied as a string object
	void printVector(T*, std::string) const noexcept;
	long int getSize() const
	{
		return n;
	}
	T* getVariables()
	{
		return x0;
	}
	T getValue() const
	{
		return f0;
	}
	//solver
	bool solver() noexcept;
};
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
}
template<class T>
MHB<T>::~MHB()
{
	delete[] dir;
	delete[] x0;
	delete[] x1;
	delete[] g0;
	delete[] g1;
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
	for (size_t i = 0; i < n; ++i)
		cout << s << '[' << i<< "]=" << vec[i] << endl;
}

template<class T>
bool MHB<T>::lineSearch() noexcept
{
	for (int i = 0; i < n; i++)
		x1[i] = x0[i] - t * g0[i];
	fi_t = ppr->valGrad(x1, g1);

	while (fi_t >= f0)
	{
		t *= RO_Left;
		if (t < std::numeric_limits<T>::epsilon())
			return false;
		for (int i = 0; i < n; i++)
			x1[i] = x0[i] - t * g0[i];
		fi_t = ppr->valGrad(x1, g1);
		//	++FuncGradEvaluations;
	}

	alpha = 0.;
	fi_alpha = f0;

	while (fi_t < fi_alpha)
	{
		alpha = t;
		fi_alpha = fi_t;
		t *= RO_Right;
		for (int i = 0; i < n; i++)
			x1[i] = x0[i] - t * g0[i];
		fi_t = ppr->valGrad(x1, g1);
		//	++FuncGradEvaluations;
	}
	//აქ უბრალოდ f1-ის მნიშვნელობის მოქექვა უკვე მოძებნილებიდან არ შეიძლება?
	for (int i = 0; i < n; i++)
		x1[i] = x0[i] - alpha * g0[i];
	f1 = ppr->valGrad(x1, g1);
	return true;
}

template<class T>
bool MHB<T>::solver() noexcept
{
	f0 = f1 = ppr->valGrad(x0, g0); 
	if (infNorm(g0)<EPS) 	return true;
	while (true)
	{
		if (f0 <= f1 || counter == 100000)
		{
			if (lineSearch() == false) return false;
			counter = 0;
		}
		///	cout << "f0=" << f0<<   "   " << infNorm(g0)  <<  endl;   ///////

		for (int i = 0; i < n; ++i)
			dir[i] = x1[i] - x0[i];
		tmp = x0; x0 = x1; x1 = tmp;
		tmp = g0; g0 = g1; g1 = tmp; 

		if (infNorm(g0) < EPS) 	return true;

		f0 = f1;
		for (int i = 0; i < n; ++i)
			x1[i] = x0[i] + dir[i] - alpha * g0[i];
		f1 = ppr->valGrad(x1, g1);
		++counter;
	}
}