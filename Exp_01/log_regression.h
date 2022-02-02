#pragma once
#include"problems.h"
#include"mnist_loader.h"
#include<numeric>
#include<chrono>
#include<random>
#include<string>
#include<fstream>

template<class T>
T sigmoid(T x)
{
	return 1. / (1. + exp(-x));
}

template <class VariablesT, class DataT>
class one_layer_log_regr : public problem<VariablesT> 
{
public:
	long int n;
	size_t n_neurons{};
	size_t n_train{};
	size_t n_test{};
	DataT** X_train{};
	DataT** X_test{};
	VariablesT* p{};     //np.matmul(X_train, z)) - y
	uint8_t** onehot_test{};
	uint8_t** onehot_train{};
	uint8_t* y{};

	std::string name{};

	long int  getSize() const { return n; }
	std::string  getName() const { return name; }

	VariablesT valGrad
	(
		VariablesT* x,
		VariablesT* g
	)
	{
		VariablesT fx(0.0);
		VariablesT tmp{};
		DataT* z{};

		for (size_t i = 0; i < n_train; ++i)
		{
			z = X_train[i];
			tmp = std::inner_product(z, z + n, x, VariablesT{});

			if (0 == y[i])
				fx += log(sigmoid(-tmp));
			else
				fx += log(sigmoid(tmp));
			p[i] = sigmoid(tmp) - y[i];
		}

		for (size_t i = 0; i < n; ++i)
		{
			g[i] = 0.;
			for (size_t j = 0; j < n_train; ++j)
				g[i] += X_train[j][i] * p[j];
		}
	
		return -fx;
	}


	void initialize(VariablesT* x)
	{
		for (int i = 0; i < n; i++)
			x[i] = .0;
	}


	void set_y(size_t currentNeuron)
	{
		std::copy(onehot_train[currentNeuron], onehot_train[currentNeuron] + n_train, y);
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
};


// Iris LOader
template <class VariablesT, class DataT>
one_layer_log_regr<VariablesT, DataT>* IrisLoader(std::string s)
{
	one_layer_log_regr<VariablesT, DataT>* a = new one_layer_log_regr<VariablesT, DataT>();
	
	size_t nTrains{ 120 };
	size_t nTests{ 30 };
	size_t nNeurons{ 3 };

	a->n = 5; 
	a->p = new VariablesT[nTrains];
	a->y = new uint8_t[nTrains];
	a->name = "Iris";
	
	a->X_train = new VariablesT* [nTrains];
	for (int i = 0; i < nTrains; ++i)
	{
		a->X_train[i] = new VariablesT[5];
	}
	a->X_test = new VariablesT* [nTests];
	for (int i = 0; i < nTests; ++i)
	{
		a->X_test[i] = new VariablesT[5];
	}

	//For each neuron
	a->onehot_train = new uint8_t* [nNeurons];

	for (int i = 0; i < nNeurons; ++i)
	{
		a->onehot_train[i] = new uint8_t[nTrains];
	}

	a->onehot_test = new uint8_t * [nNeurons];
	for (int i = 0; i < nNeurons; ++i)
	{
		a->onehot_test[i] = new uint8_t[nTests];
	}

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
	VariablesT d{};
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
			{
				a->onehot_train[0][inds[i]] = 1;
				a->onehot_train[1][inds[i]] = 0;
				a->onehot_train[2][inds[i]] = 0;
			}
			if (str == "versicolor")
			{
				a->onehot_train[0][inds[i]] = 0;
				a->onehot_train[1][inds[i]] = 1;
				a->onehot_train[2][inds[i]] = 0;
			}
			if (str == "virginica")
			{
				a->onehot_train[0][inds[i]] = 0;
				a->onehot_train[1][inds[i]] = 0;
				a->onehot_train[2][inds[i]] = 1;
			}
		}
		else
		{
			for (size_t j = 0; j < 4; ++j)
			{
				ifs >> d >> c;
				a->X_test[inds[i]-a->n_train][j] = d;
			}
			a->X_test[inds[i] - a->n_train][4] = 1.;

			ifs >> str;
			if (str == "setosa")
			{
				a->onehot_test[0][inds[i] - a->n_train] = 1;
				a->onehot_test[1][inds[i] - a->n_train] = 0;
				a->onehot_test[2][inds[i] - a->n_train] = 0;
			}
			if (str == "versicolor")
			{
				a->onehot_test[0][inds[i] - a->n_train] = 0;
				a->onehot_test[1][inds[i] - a->n_train] = 1;
				a->onehot_test[2][inds[i] - a->n_train] = 0;
			}
			if (str == "virginica")
			{
				a->onehot_test[0][inds[i] - a->n_train] = 0;
				a->onehot_test[1][inds[i] - a->n_train] = 0;
				a->onehot_test[2][inds[i] - a->n_train] = 1;
			}
		}
	}
	return a;
}

// mnist Loader
template <class VariablesT, class DataT>
one_layer_log_regr<VariablesT, DataT>* mnistLoader
(
	std::string pathToTrainData,
	std::string pathToTrainLabels,
	std::string pathToTestnData,
	std::string pathToTestLabels
)
{
	one_layer_log_regr<VariablesT, DataT>* a = new one_layer_log_regr<VariablesT, DataT>();

	size_t nTrains{ 60000 };
	size_t nTests{ 10000 };
	size_t nNeurons{ 10 };
	a->n = 785;
	a->p = new VariablesT[nTrains];
	a->y = new uint8_t[nTrains];
	a->name = "mnist";

	a->X_train = new VariablesT * [nTrains];
	for (int i = 0; i < nTrains; ++i)
	{
		a->X_train[i] = new VariablesT[785];
	}
	a->X_test = new VariablesT * [nTests];
	for (int i = 0; i < nTests; ++i)
	{
		a->X_test[i] = new VariablesT[785];
	}

	//For each neuron
	a->onehot_train = new uint8_t * [nNeurons];

	for (int i = 0; i < nNeurons; ++i)
	{
		a->onehot_train[i] = new uint8_t[nTrains];
	}

	a->onehot_test = new uint8_t * [nNeurons];
	for (int i = 0; i < nNeurons; ++i)
	{
		a->onehot_test[i] = new uint8_t[nTests];
	}

	a->n_train = nTrains;
	a->n_neurons = nNeurons;
	a->n_test = nTests;

	mnist_loader train = mnist_loader(pathToTrainData, pathToTrainLabels);
	mnist_loader test = mnist_loader(pathToTestnData, pathToTestLabels);

	//Fill X_train
	for (size_t i = 0; i < a->n_train; ++i)
	{
		train.m_images[i].push_back(1.);
		std::copy(train.m_images[i].begin(), train.m_images[i].end(), a->X_train[i]);
	}
	//Fill X_test
	for (size_t i = 0; i < a->n_test; ++i)
	{
		test.m_images[i].push_back(1.);
		std::copy(test.m_images[i].begin(), test.m_images[i].end(), a->X_test[i]);
	}

	//Fill onehot_train
	for (size_t j = 0; j < a->n_train; ++j)
	{
		for (size_t i = 0; i < nNeurons; ++i)
		{
			if (train.m_labels.at(j) == i) a->onehot_train[i][j] = 1;
			else a->onehot_train[i][j] = 0;
		}
	}

	//Fill onehot_test
	for (size_t j = 0; j < a->n_test; ++j)
	{
		for (size_t i = 0; i < nNeurons; ++i)
		{
			if (test.m_labels.at(j) == i) a->onehot_test[i][j] = 1;
			else a->onehot_test[i][j] = 0;
		}
	}

	return a;
}

//--------------------------------  Benchmarking ANN  ---------------------------------------
#include<iostream>

#include <filesystem>
using namespace std::filesystem;


void runANNDoubleTest()
{
//	freopen("results.txt", "w", stdout);
	/*
	auto a = mnistLoader<double, double>
		(
			"Data/mnist/train-images.idx3-ubyte",
			"Data/mnist/train-labels.idx1-ubyte",
			"Data/mnist/t10k-images.idx3-ubyte",
			"Data/mnist/t10k-labels.idx1-ubyte"
			);
	a->set_y(6);
	for (auto i = 0; i < a->n_train; ++i)
		cout << i<< "   " << int( a->y[i]) << endl;


		*/
//	one_layer_log_regr<double,double>* a = IrisLoader<double, double>("Data/Iris/Iris.txt");
	one_layer_log_regr<double, double>* a = 
		mnistLoader<double, double>
		(
			"Data/mnist/train-images.idx3-ubyte", 
			"Data/mnist/train-labels.idx1-ubyte", 
			"Data/mnist/t10k-images.idx3-ubyte", 
			"Data/mnist/t10k-labels.idx1-ubyte"
		);

//	for (size_t j = 0; j < a->n_neurons; ++j)
	size_t j = 0;
	{
		auto st = chrono::high_resolution_clock::now();

		problem<double>* p = a;
		a->set_y(j);

		MHB<double> minimizer{ 0.,p };
		p->initialize(minimizer.getVariables());

		bool what = minimizer.solve();
		
		auto diff = chrono::high_resolution_clock::now() - st;
		auto time = chrono::duration_cast<chrono::microseconds>(diff);

		if(!what)
			std::cout << "Solved Patially..." << std::endl;
		std::cout << "gradient inf norm: " << minimizer.infNorm(minimizer.getGradient()) << std::endl;
		std::cout << "Problem:  " << p->getName() << ",  neuron " << j << ", size: " << p->getSize() << std::endl;
		std::cout << "Average time:  " << time.count() << " microseconds" << std::endl;
		minimizer.printVector(minimizer.getVariables(), "x");

		std::cout << std::endl << std::endl;
	}

	delete a;
}





/*
void runANNLongFloatTest()
{
	using
		real = cpp_dec_float_50;

	//	freopen("results.txt", "w", stdout);

	//	one_layer_log_regr<double,double>* a = IrisLoader<double, double>("Data/Iris/Iris.txt");
	one_layer_log_regr<real, real>* a =
		mnistLoader<real, real>
		(
			"Data/mnist/train-images.idx3-ubyte",
			"Data/mnist/train-labels.idx1-ubyte",
			"Data/mnist/t10k-images.idx3-ubyte",
			"Data/mnist/t10k-labels.idx1-ubyte"
			);

	//	for (size_t j = 0; j < a->n_neurons; ++j)
	size_t j = 9;
	{
		auto st = chrono::high_resolution_clock::now();

		problem<real>* p = a;
		a->set_y(j);

		MHB<real> minimizer{ 1e-1,p };
		p->initialize(minimizer.getVariables());

		bool what = minimizer.solve();

		auto diff = chrono::high_resolution_clock::now() - st;
		auto time = chrono::duration_cast<chrono::microseconds>(diff);

		if (!what)
			std::cout << "Solved Patially..." << std::endl;
		std::cout << "gradient inf norm: " << minimizer.infNorm(minimizer.getGradient()) << std::endl;
		std::cout << "Problem:  " << p->getName() << ",  neuron " << j << ", size: " << p->getSize() << std::endl;
		std::cout << "Average time:  " << time.count() << " microseconds" << std::endl;
		minimizer.printVector(minimizer.getVariables(), "x");

		std::cout << std::endl << std::endl;
	}

	delete a;
}
*/