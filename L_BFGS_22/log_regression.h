#pragma once
#include"problems.h"
#include"mnist_loader.h"
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

template <class DataT>
class one_layer_log_regr : public problem
{
public:
	double EPS{};
	long int n;
	size_t n_neurons{};
	size_t n_train{};
	size_t n_test{};
	DataT** X_train{};
	DataT** X_test{};
	double* p{};					 //np.matmul(X_train, z)) - y
	uint8_t** onehot_test{};
	uint8_t** onehot_train{};
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
		DataT* z{};

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


// Iris Loadig 
template <class DataT>
one_layer_log_regr<DataT>* IrisLoader(std::string s)
{
	one_layer_log_regr<DataT>* a = new one_layer_log_regr<DataT>();

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
				a->X_test[inds[i] - a->n_train][j] = d;
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
template <class DataT>
one_layer_log_regr<DataT>* mnistLoader
(
	std::string pathToTrainData,
	std::string pathToTrainLabels,
	std::string pathToTestnData,
	std::string pathToTestLabels
)
{
	one_layer_log_regr<DataT>* a = new one_layer_log_regr<DataT>();

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
