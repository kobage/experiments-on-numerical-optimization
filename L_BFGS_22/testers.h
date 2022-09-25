#include<vector>
#include<chrono>
#include<algorithm>
#include <functional>
#include "L_BFGS.h"
#include "lineSearch.h"
#include"problems.h"
#include"log_regression.h"
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
	problems_container.emplace_back(make_unique<DIXMAANA>(1e-6));
	problems_container.emplace_back(make_unique<DIXMAANE>(1e-6));
	problems_container.emplace_back(make_unique<DIXMAANI>(1e-6));
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
	problems_container.emplace_back(make_unique<quadrPenalSMG>(1e-5));
}

void runUCONTests()
{
	int repNumber(1);
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

	for (int j = 0; j < repNumber; j++)
	{
		auto st = chrono::high_resolution_clock::now();
		phasePrimitives prmtv{ v[i]->getSize() };
		lineSearch lnSrch(prmtv);
		l_bgfs minimizer{ prmtv, 10 };
		v[i]->initialize(minimizer.getVariables());

		//	if(i==34)
		//		lnSrch.isGeneral = false;
		minimizer.solve(lnSrch, v[i].get());

		auto diff = chrono::high_resolution_clock::now() - st;
		auto time = chrono::duration_cast<chrono::microseconds>(diff);
		repetitions[j] = time.count();

		if (j == repNumber - 1)
		{
			std::sort(repetitions.begin(), repetitions.end());
			std::cout << "Problem:  " << v[i]->getName() << ",  size: " << v[i]->getSize() << std::endl;
			minimizer.printStats();
			std::cout << "Average time:  " << repetitions[repNumber / 2] << " microseconds" << std::endl;
			std::cout << "Gradient inf norm: " << infNorm(minimizer.getGradient(), v[i]->getSize()) << std::endl;
			std::cout << std::endl;
			printVector(minimizer.getVariables(), "x", minimizer.getSize());
		}
	}
}

//--------------------------------  Benchmarking ANN  ---------------------------------------
#include<iostream>
#include<string>
#include <filesystem>
using namespace std::filesystem;
#include <string>
#include <vector>

path findPathTo(string folder)
{
	path curPath = current_path();
	path tmp{};
	vector<path> subPaths;
	for (auto e : curPath)
	{
		tmp /= e;
		subPaths.push_back(tmp);
	}

	reverse(subPaths.begin(), subPaths.end());

	path found{};
	for (auto e : subPaths)
	{
		for (auto& p : filesystem::directory_iterator(e))
		{
			if (p.path().string().substr(p.path().string().size() - 4) == "Data")
			{
				found = p.path();
				break;
			}
		}
		if (found.string().size() > 3)
			break;
	}
	return found;
}

void runANNTest()
{
	//path to folder "Data"
	path dataPath = findPathTo("Data");

	//mnist or iris
	one_layer_log_regr<double>* a = mnistLoader<double>
		(
			dataPath.string() + "/mnist/train-images.idx3-ubyte",
			dataPath.string() + "/mnist/train-labels.idx1-ubyte",
			dataPath.string() + "/mnist/t10k-images.idx3-ubyte",
			dataPath.string() + "/mnist/t10k-labels.idx1-ubyte"
			);
	/*
	one_layer_log_regr<double>* a = IrisLoader<double>(dataPath.string() +"/Iris/Iris.txt");
	*/

	a->EPS = 1e-0;      //Tolerance

	//for given "j" neuron, j=0,1,...,a->n_neurons
	size_t j = 0;
	{
		auto st = chrono::high_resolution_clock::now();

		problem* p = a;
		a->set_y(j);
		phasePrimitives prmtv{ p->getSize() };
		l_bgfs minimizer{ prmtv, 5 };
		lineSearch lnSrch(prmtv);
		p->initialize(minimizer.getVariables());
		//	lnSrch.isGeneral = false;
		minimizer.solve(lnSrch, p);

		auto diff = chrono::high_resolution_clock::now() - st;
		auto time = chrono::duration_cast<chrono::microseconds>(diff);

		std::cout << "Problem:  " << p->getName() << ",  size: " << p->getSize() << std::endl;
		minimizer.printStats();
		std::cout << "Time:  " << time.count() << " microseconds" << std::endl;
		std::cout << "Gradient inf norm: " << infNorm(minimizer.getGradient(), p->getSize()) << std::endl;
		std::cout << std::endl;
		printVector(minimizer.getVariables(), "x", minimizer.getSize());
		std::cout << std::endl << std::endl;
	}
	delete a;
}