#include <iostream>
#include<vector>
#include <tuple>
#include <cmath>
#include <string> 
#include <algorithm>
#include <fstream>
#include <random>
#include <chrono>
#include <iterator>
#include "viterbi.h"
#include <boost\dynamic_bitset.hpp>

using namespace std;



void RandomCodes(int num)
{
	std::mt19937_64 rng;
	// initialize the random number generator with time-dependent seed
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	// initialize a uniform distribution between 0 and 1
	std::uniform_real_distribution<double> unif(0, 1);
	
	// ready to generate random numbers
	const int nSimulations = 10;
	ofstream myfile;
	myfile.open("C:\\dev\\EE4504\\Data\\randomBits.txt");

	for (int i = 0; i < num; i++)
	{
		double currentRandomNumber = unif(rng);
		if (currentRandomNumber > 0.5)
		{
			myfile << "0" << endl;
		}
		else
		{
			myfile << "1" << endl;
		}
	}
	myfile.close();
}

vector<int> getbits(long int size)
{
	vector<int> output;
	output.resize(size);


	for (int i = 0; i < size; i++)
	{
		double currentRandomNumber = ((double)rand() / (RAND_MAX));

		if (currentRandomNumber > 0.5)
		{
			output[i] = 0;
		}
		else
		{
			output[i] = 1;
		}
	}
	return output;
}

vector<int> channelNoise(vector<int> codded, float frac)
{
	vector<int> output = codded;


	int numberOFerror = (int)(output.size()*frac);
	vector<int> indexToChange;

	for (int i = 0; i < numberOFerror; i++)
	{
		int currentRandomNumber = rand() % output.size() - 1;
		indexToChange.push_back(currentRandomNumber);
	}

	for (int i = 0; i < numberOFerror; i++)
	{
		if (output[indexToChange[i]] == 1)
		{
			output[indexToChange[i]] = 0;
		}
		else
		{
			output[indexToChange[i]] = 1;
		}
	}

	return output;
}


double test(double snr)
{
	long num_bits = 10000;
	int depth = 10;
	double snrDb = snr;
	double bitErrorFrac;


	bitErrorFrac = 1.0 / (pow(10, (snrDb / 10.0)));
	bitErrorFrac = erfc(bitErrorFrac);
	vector<float> errorRate;
	vector<int> bits = getbits(num_bits);
	vector<int> noisybits(num_bits *2);
	vector<int> decoded_bits(num_bits);

	tuple<string, string> gens = { "10000101", "10101011" };
	//tuple<string, string> gens = { "101", "111" };

	tuple<vector<int>, string> conv_bits;

	cout << "start contructor" << endl;
	fec_conv test(gens, depth);
	cout << "initialised " << endl;

	try
	{
		conv_bits = test.conv_encoder(bits, "0000000");

	}
	catch (const std::exception&)
	{
		cout << "coding failed" << endl;

	};
	
	noisybits = channelNoise(get<0>(conv_bits), snr);

	try
	{
		decoded_bits = test.viterbi_decoder(noisybits, false);

	}
	catch (const std::exception& g)
	{
		cout << "decoding failed: " << g.what() << endl;

	}
	//vector<int> uncoded = channelNoise(get<0>(conv_bits), bitErrorFrac[i]);

	int errors = 0;
	for (int k = 0; k < decoded_bits.size(); k++)
	{
		if (bits[k] != decoded_bits[k])
		{
			errors++;
		}

	}


	//double biterror = (double)errors / (double)(num_bits - depth);

	return (double)errors / (double)(num_bits - depth);
}


void print(vector<int> vec)
{
	for (int i = 0; i < vec.size(); i++)
	{
		cout << vec[i] << endl;
	}
}

int main() {
	srand(time(0));
	
	vector<double> snr = { 1,2,3,4,5,6,7,8,9 };
	vector<double> BERT; 
	double t = test(0.4);
	cout << t << endl;
	//for (int i = 0; i < snr.size(); i++) {
	//	BERT.push_back();
	//}

	tuple<string, string> gens = { "10000101", "10101011" };
	tuple<string, string> ppPat = { "011", "101" };

	fec_conv test(gens,10);
	vector<int> x = { 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0,
		 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1,
		 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0,
		 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1,
		 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1,
		 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1 };
	//vector<int> lol = test.puncture(x, ppPat);

	vector<int> y = {0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,1,1,0,1,1,1,1,0,0,1,1,1,1,1,1,0,1,0,0,0,0,0,1,1,1,1,1,1,0,1,0,1,1,1,0,1,1,0,0,0,1,0,1,1,0,1,0,1,1,0,1,1,0,1,0,0,0,0,0,0,0,1,1,0,1,1,0,0,1,0,1,1,0,0,1,1,0,1,1,1,1,1,1,1,0,0,1,1,0,0,0,1,0,0,1,1,1,0,0,0,0,0,1,0,0,1,0,0,0,1,1,1,0,0,0,1,0,1,0,1,0,1,0,0,1,0,0,1,0,1,1,1,0,0,0,0,1,1,1,0,0,1,1,0,0,1,0,1,1,0,1,0,1,1,1,1,0,0,1,0,1,0,1,0,0,0,0,1,0,1,0,0,1,1,1,1,1,0,1,1,1,0,1,1,1,1,0,1,0,0,0,1,1,0,1,0,1,0,0,1,1,1,1,0,0,0,0,1,0,0,1,0,0,0,1,1,0,0,0,0,1,0,0,0,1,0,1,0,1,1,1,0,1,0,0,1 };
	/*
	vector<int> tt  ={ 0, 0, 1, 0, 1, 0, 1, 0, 2, 0, 1, 0, 0, 0, 1, 0, 1, 0 };
	vector<int>yy = { 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1 };

	vector<int> itr = { 0, 1,  2,  3,  4, 5,  6,  7,  8 };
	int in = 0;
	vector<int> testout;
	for (int i = 0; i < yy.size(); i+=2)
	{
		vector<int> t2 = { yy[i], yy[i + 1] };
		in;
		testout.push_back(test.bm_calc(itr[in], vector<int>{yy[i], yy[i+1]},true));
		in++;
	}

	print(testout);*/

	auto out = test.conv_encoder(x, "0000000");
	auto oout = test.viterbi_decoder(get<0>(out),false);
	cout << y.size() << endl;
	cout << oout.size() << endl;


	//RandomCodes(pow(10,7));
	//cout << "done" << endl;

	for (int i = 0; i < y.size(); i++)
	{
		if (y[i] == get<0>(out)[i])
		{
			cout << "0" << endl;
			//cout << i << endl;
		}
		else
		{
			cout << "1" << endl;
		}
	}

	/*
	int diff = 0;
	for (int i = 0; i < get<0>(out).size(); i++)
	{
		if (get<0>(out)[i] != y[i])
		{
			diff++;
		}
	}

	cout << diff << endl;
	*/
	return 0;
}