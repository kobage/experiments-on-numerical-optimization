#pragma once
#include"problems.h"
#include<numeric>
#include<chrono>
#include<random>
#include<string>
#include<fstream>
#include<array>
double sigmoid(double x)
{
	return 1. / (1. + exp(-x));
}

class one_layer_log_regr : public problem
{
public:
	double EPS{};
	long int n;
	size_t n_neurons{};
	size_t n_train{};
	size_t n_test{};
	double** X_train{};
	double** X_test{};
	double* p{};				//np.matmul(X_train, z)) - y
	uint8_t* labels_train;
	uint8_t* labels_test;
	uint8_t* y{};
	std::string name{};

	long int  getSize() const { return n; }
	std::string  getName() const { return name; }

	double valGrad
	(
		double* x,
		double* g
	)
	{
		double fx(0.0);
		double tmp{};
		double* z{};

		for (size_t i = 0; i < n_train; ++i)
		{
			z = X_train[i];
			tmp = std::inner_product(z, z + n, x, double{});

			if (0 == y[i])
				fx += log(sigmoid(-tmp));
			else
				fx += log(sigmoid(tmp));
			p[i] = sigmoid(tmp) - y[i];
		}

		for (size_t i = 0; i < n; ++i)
			g[i] = 0.;

		for (size_t j = 0; j < n_train; ++j)
		{
			tmp = p[j];
			for (size_t i = 0; i < n; ++i)
				g[i] += X_train[j][i] * tmp;
		}
		return -fx;
	}
	double value
	(
		double* x
	)
	{
		double fx(0.0);
		double tmp{};
		double* z{};

		for (size_t i = 0; i < n_train; ++i)
		{
			z = X_train[i];
			tmp = std::inner_product(z, z + n, x, double{});

			if (0 == y[i])
				fx += log(sigmoid(-tmp));
			else
				fx += log(sigmoid(tmp));
		}
		return -fx;
	}
	void initialize(double* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = .0;
	}

	bool stoppingCondition(double* g) const
	{
		return (infNorm(g, n) < EPS);
	}

	void set_y(size_t currentNeuron)
	{
		for (size_t j = 0; j < n_train; ++j)
		{
			if (labels_train[j] == currentNeuron) y[j] = 1;
			else y[j] = 0;
		}
	}

	~one_layer_log_regr()
	{
		delete[] p;
		delete[] y;

		for (int i = 0; i < n_train; ++i)
		{
			delete[] X_train[i];
		}
		delete[] X_train;

		for (int i = 0; i < n_test; ++i)
		{
			delete[] X_test[i];
		}
		delete[] X_test;

		delete[] labels_train;
		delete[] labels_test;
	}
};


// Iris Loadig 
one_layer_log_regr* IrisLoader(std::string s)
{
	one_layer_log_regr* a = new one_layer_log_regr();

	size_t nTrains{ 120 };
	size_t nTests{ 30 };
	size_t nNeurons{ 3 };

	a->n = 5;
	a->p = new double[nTrains];
	a->y = new uint8_t[nTrains];
	a->name = "Iris";

	a->X_train = new double* [nTrains];
	for (int i = 0; i < nTrains; ++i)
	{
		a->X_train[i] = new double[5];
	}
	a->X_test = new double* [nTests];
	for (int i = 0; i < nTests; ++i)
	{
		a->X_test[i] = new double[5];
	}

	a->labels_train = new uint8_t[nTrains];
	a->labels_test = new uint8_t[nTests];

	a->n_train = nTrains;
	a->n_neurons = nNeurons;
	a->n_test = nTests;

	std::array<int, 150> inds;
	for (size_t i = 0; i < 150; ++i)
		inds[i] = i;

	// obtain a time-based seed:
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	shuffle(inds.begin(), inds.end(), std::default_random_engine(seed));
	//	shuffle(inds.begin(), inds.end(), std::default_random_engine(0));   //may be used when debugging

	std::ifstream ifs{ s };
	double d{};
	char c{};
	std::string str{};

	//fill
	for (size_t i = 0; i < 150; ++i)
	{
		if (inds[i] < a->n_train)
		{

			for (size_t j = 0; j < 4; ++j)
			{
				ifs >> d >> c;
				a->X_train[inds[i]][j] = d;
			}
			a->X_train[inds[i]][4] = 1.;
			ifs >> str;
			if (str == "setosa")
				a->labels_train[inds[i]] = 0;
			else
				if (str == "versicolor")
					a->labels_train[inds[i]] = 1;
				else
					a->labels_train[inds[i]] = 2;
		}
		else
		{
			for (size_t j = 0; j < 4; ++j)
			{
				ifs >> d >> c;
				a->X_test[inds[i] - a->n_train][j] = d;
			}
			a->X_test[inds[i] - a->n_train][4] = 1.;

			ifs >> str;
			if (str == "setosa")
				a->labels_test[inds[i] - a->n_train] = 0;
			else
				if (str == "versicolor")
					a->labels_test[inds[i] - a->n_train] = 1;
				else
					a->labels_test[inds[i] - a->n_train] = 2;
		}
	}
	return a;
}

// mnist Loader
one_layer_log_regr* mnistLoader
(
	std::string pathToTrainData,
	std::string pathToTrainLabels,
	std::string pathToTestnData,
	std::string pathToTestLabels
)
{
	one_layer_log_regr* a = new one_layer_log_regr();

	size_t nTrains{ 60000 };
	size_t nTests{ 10000 };
	size_t nNeurons{ 10 };
	a->n = 785;
	a->p = new double[nTrains];
	a->y = new uint8_t[nTrains];
	a->name = "mnist";

	a->X_train = new double* [nTrains];
	for (int i = 0; i < nTrains; ++i)
	{
		a->X_train[i] = new double[785];
	}
	a->X_test = new double* [nTests];
	for (int i = 0; i < nTests; ++i)
	{
		a->X_test[i] = new double[785];
	}

	a->labels_train = new uint8_t[nTrains];
	a->labels_test = new uint8_t[nTests];

	a->n_train = nTrains;
	a->n_neurons = nNeurons;
	a->n_test = nTests;

	//In blocks, to reuse notion "ifs"
	//Fill X_train
	{
		std::string image_file = pathToTrainData;
		std::ifstream ifs(image_file.c_str(), std::ios::in | std::ios::binary);
		char p[4];
		ifs.read(p, 4); ifs.read(p, 4); ifs.read(p, 4); ifs.read(p, 4);
		//as we already have parameter values, so just read with zero effect
		char* q = new char[784];
		for (int i = 0; i < 60000; ++i) {
			ifs.read(q, 784);
			for (int j = 0; j < 784; ++j) {
				a->X_train[i][j] = q[j] / 255.0;
			}
			a->X_train[i][784] = 1.;
		}
		delete[] q;
	}
	//Filling "labels_train"
	{
		std::string image_file = pathToTrainLabels;
		std::ifstream ifs(image_file.c_str(), std::ios::in | std::ios::binary);
		char p[4];
		ifs.read(p, 4); ifs.read(p, 4);
		for (int i = 0; i < 60000; ++i) {
			ifs.read(p, 1);
			int label = p[0];
			a->labels_train[i] = label;   //directly p[0]? later?
		}
	}
	//Filling "X_tests"
	{
		std::string image_file = pathToTestnData;
		std::ifstream ifs(image_file.c_str(), std::ios::in | std::ios::binary);
		char p[4];
		ifs.read(p, 4); ifs.read(p, 4); ifs.read(p, 4); ifs.read(p, 4);
		//as we already have parameter values, so just read with zero effect
		char* q = new char[784];
		for (int i = 0; i < 10000; ++i) {
			ifs.read(q, 784);
			for (int j = 0; j < 784; ++j) {
				a->X_test[i][j] = q[j] / 255.0;
			}
			a->X_test[i][784] = 1.;
		}
		delete[] q;
	}
	//Filling "labels_test"
	{
		std::string image_file = pathToTestLabels;
		std::ifstream ifs(image_file.c_str(), std::ios::in | std::ios::binary);
		char p[4];
		ifs.read(p, 4); ifs.read(p, 4);
		for (int i = 0; i < 10000; ++i) {
			ifs.read(p, 1);
			int label = p[0];
			a->labels_test[i] = label;   //directly p[0]? later?
		}
	}

	return a;
}

//--------------------------------  Benchmarking ANN  ---------------------------------------
#include<iostream>
#include <filesystem>
#include"lineSearch.h"
#include"L_BFGS.h"
#include"MHB.h"
#include"MNAG.h"
#include"ADAM.h"
	using namespace cm;
using namespace std::filesystem;

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
	one_layer_log_regr* a = mnistLoader
		(
			dataPath.string() + "/mnist/train-images.idx3-ubyte",
			dataPath.string() + "/mnist/train-labels.idx1-ubyte",
			dataPath.string() + "/mnist/t10k-images.idx3-ubyte",
			dataPath.string() + "/mnist/t10k-labels.idx1-ubyte"
			);
	/*
		one_layer_log_regr<double>* a = IrisLoader<double>(dataPath.string() +"/Iris/Iris.txt");
	*/
	a->EPS = 1e-0;

	size_t j = 7;
	{
		auto st = chrono::high_resolution_clock::now();
		a->set_y(j);
		
		MHB solver(a);
		//MNAG solver{ a};
		//L_BFGS solver{ a };
		//ADAM minimizer{ a };

		cm::lineSearch lnSrch;

		a->initialize(x0);
	//	solver(lnSrch, 1);
		solver(lnSrch,1.001);
		//	minimizer.solve(a);			//only for "ADAM"
		auto diff = chrono::high_resolution_clock::now() - st;
		auto time = chrono::duration_cast<chrono::microseconds>(diff);

		std::cout << std::endl <<  "Problem:  " << a->getName() << ",  size: " << a->getSize() << ".  neuron:  " << j << std::endl;
		solver.printStats();
		std::cout << std::endl << "   time:  " << time.count() << " microseconds,  Gradient inf norm: " << infNorm(g0, n) << std::endl;
		std::cout << std::endl;
		printVector(x0, "x", n);
		std::cout << std::endl << std::endl;
	}
	delete a;
}
