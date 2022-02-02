#pragma once
#include"problems.h"
#include"log_regression.h"
#include<chrono>
#include<random>
#include<string>
#include<fstream>

// Iris Loader
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

	a->X_train = new VariablesT * [nTrains];
	for (int i = 0; i < nTrains; ++i)
	{
		a->X_train[i] = new VariablesT[5];
	}
	a->X_test = new VariablesT * [nTests];
	for (int i = 0; i < nTests; ++i)
	{
		a->X_test[i] = new VariablesT[5];
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