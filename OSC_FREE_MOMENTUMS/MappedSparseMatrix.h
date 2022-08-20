#pragma once
#include<map>
#include <vector>
#include<string>
#include<iostream>
#include <algorithm>

using namespace std;

class MpSpMtr
{
public:
	uint32_t n{};
	uint32_t nnzz;
	uint8_t p{ 10 };   //for limited memory incomplete Cholesky factorization
	vector<map<int32_t, double>> matrix, LLT;
	std::string matrixName;
	MpSpMtr() {};

	void matrixByVector(double*, double*);
	void printMapData(vector<map<int32_t, double>>);
	void printVector(double*);

	//lower triangular part is L; upper -  - - LT
	bool iCholesky();
	bool l_iCholesky();
	void LLTinverseByVector(double*, double*);
	void diagonal(double*);
};

void MpSpMtr::matrixByVector(double* x, double* res)
{
	for (size_t i = 0; i < n; ++i)
	{
		res[i] = 0.;
		for (auto& e : matrix[i])
			res[i] += x[e.first] * e.second;
	}
}

void MpSpMtr::printMapData(vector<map<int32_t, double>> m)
{
	for (auto& e : m)
	{
		for (auto& el : e)
			cout << '(' << &e - &m[0] << ',' << el.first << ',' << el.second << ")   ";
		cout << endl;
	}
	cout << endl;
}

void MpSpMtr::printVector(double* x)
{
	for (size_t i = 0; i < n; ++i)
		cout << i << ": " << x[i] << "  " << endl;
	cout << endl;
}

bool MpSpMtr::iCholesky()
{
	LLT = matrix;
	double ajj{}, ajk{}, aij{};
	std::map<int32_t, double>::iterator kit, iit;
	int i{}, k{};

	for (size_t j = 0; j < n; ++j)
	{
		auto it = LLT[j].find(j);
		if (!(it->second > 0))
			if (it->second + 0.001 > 0)
				it->second = it->second + 0.001 > 0;
			else
				return false;
		it->second = ajj = sqrt(it->second);

		//loop for k:  for k = 1:j-1 - - - k is iterator here
		for (kit = prev(it); kit != LLT[j].end(); --kit)
		{
			k = kit->first;
			ajk = kit->second;
			//loop for i:  for i = j+1:n - - - i is iterator here
			for (iit = next(it); iit != LLT[j].end(); ++iit)
			{
				i = iit->first;
				if (LLT[i].find(k) != LLT[i].end())
					LLT[i][j] -= LLT[i][k] * ajk;
			}
		}
		for (iit = next(it); iit != LLT[j].end(); ++iit)
		{
			i = iit->first;
			iit->second = aij = LLT[i][j] /= ajj;
			LLT[i][i] -= aij * aij;
		}
	}
	return true;
}
//-----------------------------------------------------
bool MpSpMtr::l_iCholesky()
{
	double ajj{};
	std::map<int32_t, double>::iterator kit, iit;
	int i{}, k{};

	for (size_t j = 0; j < n; ++j)
	{
		std::map<int, double> tmpMap;
		vector<pair<int32_t, double>> to_manipulate;
		auto it = matrix[j].find(j);
		if (!(it->second > 0))
			if (it->second + 1e-07 > 0)
				it->second = it->second + 1e-07 > 0;
			else
			{
				cout << "Oops!" << endl;   //////////////////////////
				return false;
			}
		LLT[j][j] = ajj = sqrt(it->second);
		//Identifying col_len and saving aij - values in temporary map<> object
		int col_len{ 0 }, counter{ 0 };
		for (iit = next(it); iit != matrix[j].end(); ++iit)
		{
			tmpMap[iit->first] = iit->second;
			++col_len;
		}

		//a(i,k)*a(j,k) as a(k,j)*a(k,i)
		for (kit = prev(LLT[j].find(j)); kit != LLT[j].end(); --kit)
		{
			k = kit->first;
			for (iit = next(LLT[k].find(j)); iit != LLT[k].end(); ++iit)
				tmpMap[iit->first] -= iit->second * kit->second;
		}
		//----------------------------------------------------------------------------------------//
		//aij 
		for (auto iter = tmpMap.begin(); iter != tmpMap.end(); ++iter)
		{
			iter->second /= ajj;
			to_manipulate.push_back(*iter);
		}

		//Find largest col_len+p elements in a(:,j) and store; odifying aii's
		if (tmpMap.size() <= col_len + p)
			for (auto iter = tmpMap.begin(); iter != tmpMap.end(); ++iter)
			{
				i = iter->first;
				LLT[i][j] = LLT[j][i] = iter->second;
				LLT[i][i] -= pow(iter->second, 2);
			}
		else
		{
			auto cmp = [](const pair<int32_t, double>& a, const pair<int32_t, double>& b)
			{return (std::fabs(a.second) < std::fabs(b.second)); };
			//	sort(to_manipulate.rbegin(), to_manipulate.rend(), cmp);
			nth_element(to_manipulate.rbegin(), to_manipulate.rbegin() + to_manipulate.size() - col_len - p, to_manipulate.rend(), cmp);
			for (auto ii = to_manipulate.begin(); ii < to_manipulate.begin() + col_len + p; ++ii)
			{
				i = ii->first;
				LLT[i][j] = LLT[j][i] = ii->second;
				LLT[i][i] -= pow(ii->second, 2);
			}
		}
	}
	return true;
}
/*     */
void MpSpMtr::LLTinverseByVector(double* b, double* x)
{
	map<int32_t, double>::iterator it, kit;
	double dif{ 0. };
	x[0] = b[0] / LLT[0][0];

	for (size_t i = 1; i < n; ++i)
	{
		it = LLT[i].find(i);
		dif = b[i];
		for (kit = prev(it); kit != LLT[i].end(); --kit)
			dif -= x[kit->first] * kit->second;
		x[i] = dif / it->second;
	}

	x[n - 1] /= LLT[n - 1][n - 1];
	for (int i = n - 2; i >= 0; --i)
	{
		it = LLT[i].find(i);
		dif = x[i];
		for (kit = next(it); kit != LLT[i].end(); ++kit)
			dif -= x[kit->first] * kit->second;
		x[i] = dif / it->second;
	}
}
void MpSpMtr::diagonal(double* d)
{
	for (size_t i = 0; i < n; ++i)
		d[i] = matrix[i][i];
}