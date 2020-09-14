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
#include <boost\dynamic_bitset.hpp>



class trellis_nodes {

public:
	int Ns;
	std::vector<int> fn;
	std::vector<int> tn;
	std::vector<int> out_bits;

	trellis_nodes(int NS = 10) {
		Ns = NS;
		fn.resize(NS);
		tn.resize(NS);
		out_bits.resize(NS);
	}
};


class trellis_branches {
public:
	int Ns;
	std::vector<int>states1; 
	std::vector<int>states2;
	std::vector<int>bits1;
	std::vector<int>bits2;
	std::vector<int>input1;
	std::vector<int>input2;

	trellis_branches(int NS) {
		Ns = NS;
		states1.resize(Ns);
		states2.resize(Ns);
		bits1.resize(Ns);
		bits2.resize(Ns);
		input1.resize(Ns);
		input2.resize(Ns);
	}
};

class trellis_paths {
public:
	int Ns;
	int decision_depth;
	std::vector<std::vector<int>> traceback_states;
	std::vector<std::vector<int>> cumulative_metric;
	std::vector<std::vector<int>> traceback_bits;

	trellis_paths(int NS = 10,int D = 10) {
		Ns = NS;
		decision_depth = D;
		std::vector<std::vector<int>> ntemp (NS, std::vector<int>(decision_depth, 0));
		traceback_states = ntemp;
		cumulative_metric = ntemp;
		traceback_bits = ntemp;
	}

};

/*
this is the class for doing the decoidng and endcoidng using the viterbi encoder and decoder
*/
class fec_conv {
private:
	//returns a vector containng all the indexes where vec[i] = key
	std::vector<int> Where(std::vector<int> vec, int key)
	{
		std::vector<int> outPut;
		for (int i = 0; i < vec.size(); i++)
		{
			if (vec[i] == key)
			{
				outPut.push_back(i);
			}
		}
		return outPut;
	}

	// returns a collomn of a 2d vector by collomn index
	std::vector<int> GetColomn(std::vector<std::vector<int>> vec, int col) {
		std::vector<int> outPut;

		for (int i = 0; i < vec.size(); i++)
		{
			outPut.push_back(vec[i][col]);
		}
		return outPut;
	}

	// prints vector 
	void printVector(std::vector<int> vec)
	{

		for (int j = 0; j < vec.size(); j++)
		{
			std::cout << vec[j] << ", ";
		}
		std::cout << std::endl;
	}

	// prints 2d vector 
	void print2DVector(std::vector<std::vector<int>> vec)
	{
		for (int i = 0; i < vec.size(); i++)
		{
			for (int j = 0; j < vec[i].size(); j++)
			{
				std::cout << vec[i][j] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
		std::cout << std::endl;
	}

	// sets a collumn of a 2d vector to a set value
	std::vector<std::vector<int>>SetColomn(std::vector<std::vector<int>> vec, std::vector<int> inp, int col) {
		std::vector<std::vector<int>> outPut;

		for (int i = 0; i < vec.size(); i++)
		{
			vec[i][col] = inp[i];
		}

		outPut = vec;

		return outPut;
	}

	// removes a colomn from a 2d vector 
	std::vector<std::vector<int>> RemoveColomn(std::vector<std::vector<int>> vec, int col)
	{
		for (int i = 0; i < vec.size(); i++)
		{
			vec[i].erase(vec[i].begin() + col);
		}
		return vec;
	}

	// pops a colomn from a 2d vector 
	std::vector<std::vector<int>> popColumn(std::vector<std::vector<int>> vec)
	{
		std::vector<std::vector<int>> output = vec;
		for (int i = 0; i < output.size(); i++)
		{
			output[i].pop_back();
		}
		return output;
	}

	// converts a vector representing a 2 bit number to an integer
	int vecToInt(std::vector<int> vec)
	{
		std::reverse(vec.begin(), vec.end());

		int result = 0;

		for (int i = 0; i < vec.size(); i++)
		{
			result += (pow(10, i) * vec[i]);
		}

		return(result);
	}

	// calbulated the hamming distance bitween 2 two-bit numbers 
	int hammingDistance(int int1, int int2, int size)
	{

		//return bitset<T>(int1^int2).count();
		int distance = boost::dynamic_bitset<>(size, int1^int2).count();
		return distance;

	}

	// index of for vetors 
	std::vector<int> IndexOf(std::vector<int> indexes, std::vector<int> vec)
	{
		std::vector<int> outPut;
		for (int i = 0; i < indexes.size(); i++)
		{
			outPut.push_back(vec[indexes[i]]);
		}
		return outPut;
	}

	// convers a int to a string of bits in 2 bit form 
	std::string convert_to_binary_string(const unsigned long long int value,
		bool skip_leading_zeroes = false, int ContsLength = 8)
	{
		std::string str;
		bool found_first_one = false;
		const int bits = sizeof(unsigned long long) * 8;  // Number of bits in the type

		for (int current_bit = bits - 1; current_bit >= 0; current_bit--)
		{
			if ((value & (1ULL << current_bit)) != 0)
			{
				if (!found_first_one)
					found_first_one = true;
				str += '1';
			}
			else
			{
				if (!skip_leading_zeroes || found_first_one)
					str += '0';
			}
		}

		if (str.length() < ContsLength)
		{
			while(str.length() != ContsLength)
			{
				str = '0' + str;
			}
			
		}

		return str;
	}
public:
	std::tuple<std::string, std::string> G_polys;
	int constraint_length;
	int Nstates;
	int decision_depth;
	trellis_nodes* input_zero;
	trellis_nodes* input_one;
	trellis_paths* paths;
	trellis_branches* branches;

	// constuctor initialises the data structures 
	fec_conv(std::tuple<std::string, std::string> G, int Depth = 10) {
		G_polys = G;
		constraint_length = std::get<0>(G).length();
		Nstates = pow(2, constraint_length - 1);
		decision_depth = Depth;
		input_zero = new trellis_nodes(Nstates);
		input_one = new trellis_nodes(Nstates);
		paths = new trellis_paths(Nstates, decision_depth);
		branches = new trellis_branches(Nstates);

		for (int i = 0; i < Nstates; i++)
		{
			input_zero->fn[i] = i;
			input_one->fn[i] = i;
			std::string binstring = "";
		
			binstring = convert_to_binary_string(i, true, constraint_length - 1);
	

			auto conv0 = conv_encoder(std::vector<int>{0}, binstring);
			auto conv1 = conv_encoder(std::vector<int>{1}, binstring);
			std::vector<int> output0 = std::get<0>(conv0);
			std::string state0 = std::get<1>(conv0);
			std::vector<int> output1 = std::get<0>(conv1);
			std::string state1 = std::get<1>(conv1);

			input_zero->tn[i] = stoi(state0, 0, 2);
			input_one->tn[i] = stoi(state1, 0, 2);

			input_zero->out_bits[i] = 2 * output0[0] + output0[1];
			input_one->out_bits[i] = 2 * output1[0] + output1[1];
		}

		for (int i = 0; i < Nstates; i++)
		{
			std::vector<int> match_zero_idx = Where(input_zero->tn, i);
			std::vector<int> match_one_idx = Where(input_one->tn, i);
			if (match_zero_idx.size() != 0)
			{
				int t1 = IndexOf(match_zero_idx, input_zero->fn)[0];
				int t2 = IndexOf(match_zero_idx, input_zero->fn)[1];
				int t3 = IndexOf(match_zero_idx, input_zero->out_bits)[0];
				int t4 = IndexOf(match_zero_idx, input_zero->out_bits)[1];
				branches->states1[i] = t1;
				branches->states2[i] = t2;
				branches->bits1[i] = t3;
				branches->bits2[i] = t4;
				branches->input1[i] = 0;
				branches->input2[i] = 0;
			}
			else if (match_one_idx.size() != 0)
			{
				int t1 = IndexOf(match_one_idx, input_one->fn)[0];
				int t2 = IndexOf(match_one_idx, input_one->fn)[1];
				int t3 = IndexOf(match_one_idx, input_one->out_bits)[0];
				int t4 = IndexOf(match_one_idx, input_one->out_bits)[1];
				branches->states1[i] = t1;
				branches->states2[i] = t2;
				branches->bits1[i] = t3;
				branches->bits2[i] = t4;
				branches->input1[i] = 1;
				branches->input2[i] = 1;
			}
		}
	}
	// calculates the branch metric 
	int bm_calc(int ref_code_bits, std::vector<int> rec_code_bits, bool metric_type)
	{
		std::string bits;
		int ref_MSB;
		int ref_LSB;
		int distance;

		if (metric_type)
		{
			bits = convert_to_binary_string(ref_code_bits, true, constraint_length - 1);
			if (ref_code_bits == 0)
			{
				for (int i = 0; i < constraint_length - 2; i++)
					bits = "0" + bits;
			}
			std::string  s(1, bits[0]);
			std::string  s1(1, bits[1]);
			ref_MSB = stoi(s, 0, 2);
			ref_LSB = stoi(s1, 0, 2);
			distance = pow((rec_code_bits[0] - ref_MSB), 2);
			distance = distance + pow((rec_code_bits[1] - ref_LSB), 2);

		}
		else {
			unsigned long s = Nstates - 1;
			const size_t size = sizeof(s) * s;
			//(ref_code_bits^vecToInt(rec_code_bits)).count();

			distance = hammingDistance(ref_code_bits, vecToInt(rec_code_bits), s);

		}

		return distance;
	}

	// decodes the bit stream
	std::vector<int> viterbi_decoder(std::vector<int> x, bool metric_type = true)
	{

		std::vector<int> cm_present;
		cm_present.resize(Nstates);
		int NS = x.size();
		std::vector<int> y(NS - decision_depth);
		y.resize(NS - decision_depth);
		//std::cout << y.size() << std::endl;
		int d1, d2, k = 0;
		std::vector<int> test_lol;
		std::vector<std::vector<int>>  test;
		std::vector<int> cm_past;
		std::vector<std::vector<int>> old_tb_states_temp = paths->traceback_states;
		std::vector<std::vector<int>> old_tb_bits_temp = paths->traceback_bits;

		for (int i = 0; i < NS; i += 2)
		{
			//for (int w = 0; w < ((float)i / (float)NS)*100.0; w++) {
			//	std::cout << "|";
			//}
			//std::cout << ((float)i/ (float)NS)*100.0 << std::endl;
			//std::cout<<std::endl;
			cm_past = GetColomn(paths->cumulative_metric, 0);
			std::vector<std::vector<int>> tb_states_temp = old_tb_states_temp;
			std::vector<std::vector<int>> tb_bits_temp = old_tb_bits_temp;
			tb_states_temp = popColumn(paths->traceback_states);
			tb_bits_temp = popColumn(paths->traceback_bits);
			//std::cout << "1" << std::endl;
			for (int j = 0; j < Nstates; j++)
			{

				d1 = bm_calc(branches->bits1.at(j), std::vector<int>{x.at(i), x.at(i + 1)}, metric_type);
				d1 = d1 + cm_past.at(branches->states1.at(j));

				d2 = bm_calc(branches->bits2.at(j), std::vector<int>{x.at(i), x.at(i + 1)}, metric_type);
				d2 = d2 + cm_past.at(branches->states2.at(j));

				if (d1 <= d2)
				{
					//cout << "d1" << endl;
					cm_present[j] = (d1);
					std::vector<int> tempTBS = tb_states_temp.at(branches->states1.at(j));
					tempTBS.insert(tempTBS.begin(), branches->states1.at(j));
					paths->traceback_states.at(j) = tempTBS;


					std::vector<int> tempTBB = tb_bits_temp.at(branches->states1.at(j));
					tempTBB.insert(tempTBB.begin(), branches->input1.at(j));
					paths->traceback_bits.at(j) = tempTBB;

					//cout << branches->states1[j] << endl;;
					//cout << branches->bits1[j] << endl;;

				}
				else
				{
					//cout << "d2" << endl;

					cm_present[j] = (d2);
					std::vector<int> tempTBS = tb_states_temp.at(branches->states2.at(j));
					tempTBS.insert(tempTBS.begin(), branches->states2.at(j));
					paths->traceback_states.at(j) = tempTBS;


					std::vector<int> tempTBB = tb_bits_temp[branches->states2[j]];
					tempTBB.insert(tempTBB.begin(), branches->input2[j]);
					paths->traceback_bits.at(j) = tempTBB;

					//cout << branches->states2[j] << endl;;
					//cout << branches->bits2[j] << endl;;
				}
			}
			std::vector<std::vector<int>> tempCumPopedBack = popColumn(paths->cumulative_metric);

			for (int l = 0; l < Nstates; l++)
			{
				std::vector<int> temp = tempCumPopedBack.at(l);
				temp.insert(temp.begin(), cm_present.at(l));
				paths->cumulative_metric.at(l) = temp;
			}

			std::vector<int> col = GetColomn(paths->cumulative_metric, 0);
			auto min_metric = min_element(col.begin(), col.end());
			std::vector<int> min_indx = Where(col, *min_metric);
			if (i >= 2 * decision_depth - 2)
			{
				//int shi = paths->traceback_bits[min_metric[0]][paths->traceback_bits[min_metric[0]].size() - 1];

				//cout << shi << ", " << i << endl;
				int min = min_indx[0];
				/*if (min > constraint_length)
				{
					min = std::rand()% (constraint_length);
				}*/
				int val = paths->traceback_bits[min].size() - 1;
		
				y[k] = paths->traceback_bits[min][val];
				k++;
			}

			//print2DVector(paths->traceback_bits);
			//std::cout << i << std::endl;

		}
		y.resize(k);
		//printVector(test_lol);
		return y;

	}

	// encode the bit stream
	std::tuple<std::vector<int>, std::string> conv_encoder(std::vector<int> input, std::string state)
	{
		std::vector<int> output;
		std::string State;
		int k = 0, l = 1;
		for (int i = 0; i < input.size(); i++)
		{
			int u1 = input[i];
			int u2 = input[i];
			for (int j = 1; j < constraint_length - 1; j++) {
				int  gpoly0 = (int)std::get<0>(G_polys)[j] - 48;
				int  gpoly1 = (int)std::get<1>(G_polys)[j] - 48;
				int stateNum = ((int)state[j - 1] - 48);
				if (gpoly0 == 1)
					u1 = u1 ^ stateNum;
				if (gpoly1 == 1)
					u2 = u2 ^ stateNum;
			}
			output.push_back(u1);
			output.push_back(u2);

			std::string newSate = state.substr(0, state.size() - 1);
			state = std::to_string(input[i]);
			state = state + newSate;

		}
		return std::tuple<std::vector<int>, std::string> {output, state};
	}

	// puncturing for the pupose of increase the bit rate 
	std::vector<int> puncture(std::vector<int> code_bits, std::tuple<std::string, std::string> puncture_pattern = {"101", "111"})
	{
		std::vector<int> output;
		std::string pp1 = std::get<0>(puncture_pattern);
		std::string pp2 = std::get<1>(puncture_pattern);

		int lenPP = std::get<0>(puncture_pattern).length();
		int num_codewords = floor(code_bits.size() / 2.0);
		std::vector<int> bits = code_bits;
		bits.resize(num_codewords * 2);
		std::vector<int> x_g1;
		std::vector<int> x_g2;
		for (int i = 0; i < bits.size(); i++)
		{
			if (i % 2 == 0)
			{
				x_g1.push_back(bits[i]);
			}
			else
			{
				x_g2.push_back(bits[i]);
			}
		}


		int num_pp = (int)floor(num_codewords / (float)lenPP);
		x_g1.resize(lenPP*num_pp);
		x_g2.resize(lenPP*num_pp);

		std::vector<int> punx_g1;
		std::vector<int> punx_g2;

		for (int i = 0; i < x_g1.size(); i += lenPP)
		{
			for (int k = 0; k < pp1.size(); k++)
			{
				if (pp1[k] == '1')
				{
					punx_g1.push_back(x_g1[i + k]);
				}
			}
			for (int k = 0; k < pp2.size(); k++)
			{
				if (pp2[k] == '1')
				{
					punx_g2.push_back(x_g2[i + k]);
				}
			}
		}
		int k = 0;
		int j = 0;
		for (int i = 0; i < punx_g1.size() + punx_g2.size(); i++)
		{
			if (i % 2 == 0)
			{
				output.push_back(punx_g1[k]);
				k++;
			}
			else
			{
				output.push_back(punx_g2[j]);
				j++;
			}
		}
		return output;
	}


};