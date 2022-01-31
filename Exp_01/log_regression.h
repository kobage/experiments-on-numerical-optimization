#pragma once
#include"problems.h"
#include<numeric>
#include<chrono>
#include<random>
#include<string>
#include<fstream>

double sigmoid(double x)
{
	return 1. / (1. + exp(-x));
}


template <class VariablesT, class DataT>
class one_layer_log_regr : public problem<VariablesT> 
{
	long int n{ 5000 };
	DataT** X_train{};
	DataT** X_test{};
	VariablesT* p{};     //np.matmul(X_train, z)) - y
	uint8_t** onehot_test{};
	uint8_t** onehot_train{};
	uint8_t* y{};

	size_t n_neurons{};
	size_t n_weiths{};
	size_t n_train{};
	size_t n_test{};

	std::string name{};
public:
	//For Iris
	one_layer_log_regr(std::string s);
	~one_layer_log_regr();
	long int virtual getSize() const { return n; }
	std::string virtual getName() const { return name; }
	
	VariablesT valGrad
	(
		VariablesT* x,
		VariablesT* g
	)
	{
		VariablesT fx(0.0);
		VariablesT tmp{};
		VariablesT* z{};

		for (size_t i = 0; i < n_train; ++i)
		{
			z = X_train[i];
			tmp = std::inner_product(z, z + n_train, x, VariablesT{});
			if (0 == y[i])
				fx += log(sigmoid(-tmp));
			else
				fx += log(sigmoid(tmp));
			p[i] = sigmoid(tmp) - y[i];
		}

		for (size_t i = 0; i < n_weiths; ++i)
		{
			g[i] = 0.;
			for (size_t j = 0; j < n_train; ++j)
				g[i] += X_train[j][i] * p[j];
		}
		return -fx;
	}
	void initialize(VariablesT* x)
	{
		for (int i = 0; i < n_train; i++)
			x[i] = .0;
	}
	void set_y(size_t currentNeuron)
	{
		std::inner_product(onehot_train[currentNeuron], onehot_train[currentNeuron] + n_train, y, VariablesT{});
	}
};

template <class VariablesT, class DataT>
one_layer_log_regr<VariablesT, DataT>::~one_layer_log_regr()
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

	for (int i = 0; i < n_neurons; ++i)
	{
		delete[] onehot_train[i];
	}
	delete[] onehot_train;

	for (int i = 0; i < n_neurons; ++i)
	{
		delete[] onehot_test[i];
	}
	delete[] onehot_test;
}
//For Iris
template <class VariablesT, class DataT>
one_layer_log_regr<VariablesT, DataT>::one_layer_log_regr(std::string s)
{
	p = new VariablesT[n_train];
	y = new uint8_t[n_train];

	X_train = new VariablesT * [n_train];
	for (int i = 0; i < n_train; ++i)
	{
		X_train[i] = new VariablesT[5];
	}
	X_test = new VariablesT * [n_test];
	for (int i = 0; i < n_test; ++i)
	{
		X_test[i] = new VariablesT[5];
	}

	//For each neuron
	onehot_train = new uint8_t* [n_neurons];
	for (int i = 0; i < n_neurons; ++i)
	{
		onehot_train[i] = new uint8_t[n_train];
	}
	onehot_test = new uint8_t * [n_neurons];
	for (int i = 0; i < n_neurons; ++i)
	{
		onehot_test[i] = new uint8_t[n_test];
	}

	std::array<int, 150> inds;
	for (size_t i = 0; i < 150; ++i)
		inds[i] = i;
	// obtain a time-based seed:
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	shuffle(inds.begin(), inds.end(), std::default_random_engine(seed));
	//	shuffle(inds.begin(), inds.end(), std::default_random_engine(0));

	std::ifstream ifs{ s };
	VariablesT d{};
	char c{};
	std::string str{};
	//At first, all data are X_train
	for (size_t i = 0; i < 150; ++i)
	{
		if (inds[i] < n_train)
		{
			for (size_t j = 0; j < 4; ++j)
			{
				ifs >> d >> c;
				X_train[inds[i]][j] = d;
			}
			X_train[inds[i]][4] = 1.;
			std::cout << X_train[inds[i]][0] << " " << X_train[inds[i]][1] << " " << X_train[inds[i]][2] << " " << X_train[inds[i]][3] << " " << std::endl;//////////////
			ifs >> str;
			if (str == "setosa")
			{
				onehot_train[0][inds[i]] = 1;
				onehot_train[1][inds[i]] = 0;
				onehot_train[2][inds[i]] = 0;
			}
			if (s == "versicolor")
			{
				onehot_train[0][inds[i]] = 0;
				onehot_train[1][inds[i]] = 1;
				onehot_train[2][inds[i]] = 0;
			}
			if (s == "virginica")
			{
				onehot_train[0][inds[i]] = 0;
				onehot_train[1][inds[i]] = 0;
				onehot_train[2][inds[i]] = 1;
			}
		}
		else
		{
			for (size_t j = 0; j < 4; ++j)
			{
				ifs >> d >> c;
				X_test[inds[i]-n_train][j] = d;
			}
			X_test[inds[i] - n_train][4] = 1.;

			ifs >> str;
			if (str == "setosa")
			{
				onehot_test[0][inds[i] - n_train] = 1;
				onehot_test[1][inds[i] - n_train] = 0;
				onehot_test[2][inds[i] - n_train] = 0;
			}
			if (s == "versicolor")
			{
				onehot_test[0][inds[i] - n_train] = 0;
				onehot_test[1][inds[i] - n_train] = 1;
				onehot_test[2][inds[i] - n_train] = 0;
			}
			if (s == "virginica")
			{
				onehot_test[0][inds[i] - n_train] = 0;
				onehot_test[1][inds[i] - n_train] = 0;
				onehot_test[2][inds[i] - n_train] = 1;
			}
		}
	}
}

//--------------------------------  Benchmarking ANN  ---------------------------------------
#include<iostream>
void runANNDoubleTest()
{
//	freopen("results.txt", "w", stdout);
//	cout << 
	auto st = chrono::high_resolution_clock::now();
	one_layer_log_regr<double, double>* a = new one_layer_log_regr<double, double>("Data/Iris/Iris.txt");
	problem<double>* p = a;
	a->set_y(0);
	
	MHB<double> minimizer{ 1e-5,p };
	p->initialize(minimizer.getVariables());

	bool what = minimizer.solve();

	auto diff = chrono::high_resolution_clock::now() - st;
	auto time = chrono::duration_cast<chrono::microseconds>(diff);

	std::cout << "Problem:  " << p->getName() << ",  size: " << p->getSize() << std::endl;
	std::cout << "Average time:  " << time.count() << " microseconds" << std::endl;
	minimizer.printVector(minimizer.getVariables(), "x");
}