#pragma once
#include"Problems.h"
#include<vector>
#include<chrono>
#include<algorithm>
#include"cm_namespace.h"
#include"L_BFGS.h"
using namespace std;

vector<unique_ptr<problem>> problems_container;

void makeDoubleTestsVector(void)
{
	problems_container.emplace_back(make_unique<ARWHEAD>(1e-6));
	problems_container.emplace_back(make_unique<BDQRTIC>(1e-6));
	problems_container.emplace_back(make_unique<BROYDN7D>(1e-6));
	problems_container.emplace_back(make_unique<BRYBND>(1e-6));
	problems_container.emplace_back(make_unique<CHAINWOO>(1e-6));
	problems_container.emplace_back(make_unique<COSINE>(1e-6));
	problems_container.emplace_back(make_unique<CRAGGLVY>(1e-6));
	problems_container.emplace_back(make_unique<CURLY10>(1e-6));
	problems_container.emplace_back(make_unique<CURLY20>(1e-6));
	problems_container.emplace_back(make_unique<CURLY30>(1e-6));
	problems_container.emplace_back(make_unique<DIXMAANA>(1e-6));
	problems_container.emplace_back(make_unique<DIXMAANB>(1e-6));
	problems_container.emplace_back(make_unique<DIXMAANC>(1e-6));
	problems_container.emplace_back(make_unique<DIXMAAND>(1e-6));
	problems_container.emplace_back(make_unique<DIXMAANE>(1e-6));
	problems_container.emplace_back(make_unique<DIXMAANF>(1e-6));
	problems_container.emplace_back(make_unique<DIXMAANG>(1e-6));
	problems_container.emplace_back(make_unique<DIXMAANH>(1e-6));
	problems_container.emplace_back(make_unique<DIXMAANI>(1e-6));
	problems_container.emplace_back(make_unique<DIXMAANJ>(1e-6));
	problems_container.emplace_back(make_unique<DIXMAANK>(1e-6));
	problems_container.emplace_back(make_unique<DIXMAANL>(1e-6));
	problems_container.emplace_back(make_unique<DIXON3DQ>(1e-6));
	problems_container.emplace_back(make_unique<DQDRTIC>(1e-6));
	problems_container.emplace_back(make_unique<DQRTIC>(1e-6));
	problems_container.emplace_back(make_unique<EDENSCH>(1e-6));
	problems_container.emplace_back(make_unique<EG2>(1e-6));
	problems_container.emplace_back(make_unique<ENGVAL1>(1e-6));
	problems_container.emplace_back(make_unique<EXTROSNB>(1e-6));
	problems_container.emplace_back(make_unique<FLETCHR>(1e-6));
	problems_container.emplace_back(make_unique<FREUROTH>(1e-6));
	problems_container.emplace_back(make_unique<GENHUMPS>(1e-6));
	problems_container.emplace_back(make_unique<GENROSE>(1e-6));
	problems_container.emplace_back(make_unique<LIARWDH>(1e-6));
	problems_container.emplace_back(make_unique<MOREBV>(1e-6));
	problems_container.emplace_back(make_unique<NONCVXU2>(1e-6));
	problems_container.emplace_back(make_unique<NONDIA>(1e-6));
	problems_container.emplace_back(make_unique<NONDQUAR>(1e-6));
	problems_container.emplace_back(make_unique<PENALTY1>(1e-6));
	problems_container.emplace_back(make_unique<PENALTY2>(1e-6));
	problems_container.emplace_back(make_unique<POWER>(1e-6));
	problems_container.emplace_back(make_unique<SROSENBR>(1e-6));
	problems_container.emplace_back(make_unique<TRIDIA>(1e-6));
	problems_container.emplace_back(make_unique<Woods>(1e-6));
	problems_container.emplace_back(make_unique<quadrPenalSMG>(1e-6));
}

#include"lineSearch.h"
#include"L_BFGS.h"
#include"MHB.h"
#include"MNAG.h"
#include"ADAM.h"
using namespace cm;

void runUCONTests()
{
	int repNumber(10);
	makeDoubleTestsVector();

	vector<unique_ptr<problem>>& v = problems_container;
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

	std::chrono::high_resolution_clock::time_point  st;  //
	std::chrono::high_resolution_clock::duration  diff;  //
	for (int j = 0; j < repNumber; j++)
	{
		st = chrono::high_resolution_clock::now();
		cm::lineSearch lnSrch;


	//	cm::L_BFGS solver{ v[i]->getSize() , 100 }; 
	//	cm::L_BFGS solver{ v[i].get() };			//default, with: the number of BFGS corrections saved = 11
	//	cm::MHB solver{ v[i].get() };
		cm::MNAG solver{ v[i].get() };
	//	cm::ADAM solver{ v[i].get() };

		v[i]->initialize(x0);
	//	solver(lnSrch);

		//solver(lnSrch,1.001);			//for MNAG amd MHB
		solver(lnSrch, 1);			//for MNAG

		diff = chrono::high_resolution_clock::now() - st;
		//	auto time = chrono::duration_cast<chrono::microseconds>(diff);
		repetitions[j] = std::chrono::duration_cast<std::chrono::microseconds>(diff).count();

		if (j == repNumber - 1)
		{
			std::sort(repetitions.begin(), repetitions.end());
			std::cout << i + 1 << ". Problem:  " << v[i]->getName() << ",  size: " << v[i]->getSize() << std::endl;
			solver.printStats();
			std::cout << "   Average time:  " << repetitions[repNumber / 2] << " microseconds"  
			 << ", Gradient inf norm: " << infNorm(g0, v[i]->getSize()) << std::endl;
			std::cout << std::endl;
			printVector(x0, "x", n);
		}
	}
}
