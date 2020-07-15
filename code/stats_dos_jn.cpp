#include<iostream>
#include<fstream>
#include<algorithm>
#include<vector>
#include<random>
#include<chrono>
#include<string>
#include<omp.h>

#include "parameters.hpp"
#include "io_data.hpp"

using namespace std;

/*-----------Parameters-------------*/
const int nblocks = 1e2;											//# Jackknife blocks
const int kmeas = params::nmeas/(double)nblocks;		//# measurements in one block														// V/beta ratio

const string FOLDER_BASE = "Nt"+to_string(params::Nt)+"_Ns"+to_string(params::Ns)+"_beta"+to_string((int)(params::beta*10))+"/"; 												//# intervals
const string FOLDER_IN = "../sim/data/"+FOLDER_BASE;
const string FOLDER_OUT = "../sim/stats/"+FOLDER_BASE;
const string FILE_NAME_IN = "_RE_Q_lambda.dat";
const string FILE_NAME_lvals = "_lambda_vals.dat";
const string FILE_NAME_OUT1 = "_RE_Q_lambda_jnstats.dat";
const string FILE_NAME_OUT2 = "_RE_Q_lambda_jnsamples.dat";

string FOLDER_IN_Q_PRC = "../prc/start_configs/Nt"+to_string(params::Nt)+"_Ns"+
                    to_string(params::Ns)+"_betastart"+to_string((int)(params::beta_start*10))+"/";
string FOLDER_IN_Q_SIM = "../sim/start_configs/Nt"+to_string(params::Nt)+"_Ns"+
                    to_string(params::Ns)+"_betastart"+to_string((int)(params::beta_start*10))+"/";
string FOLDER_IN_DATA = "../sim/data/Nt"+to_string(params::Nt)+"_Ns"+
                  	to_string(params::Ns)+"_beta"+to_string((int)(params::beta*10))+"/";
string FOLDER_OUT_DATA = "../sim/stats/Nt"+to_string(params::Nt)+"_Ns"+
							      to_string(params::Ns)+"_beta"+to_string((int)(params::beta*10))+"/";
/*----------------------------------*/

/*-------Functions-------*/
void read_data(const string& FILE, vector<double>& data);
double get_avg(const vector<double>& data);
void write_data(ofstream& outfile_dat, const vector<double>& lambda_vals, const vector<double>& data);
/*-----------------------*/

int main()
{
	cout << "Hello World - This is a statistics program for OBC SU(2) DoS!" << endl << endl;
	cout << "Method: Jackknife with binning" << endl << endl;
	cout << "The program generates two files per interval: One file contains the rescaled restricted expectation value (<<X>>_n + x_n)/delta + 0.5 & Jackknife errors. \nAdditionally, the program generates one file where mean function & Jackknife samples are stored. Subsequently, these files should be passed to a fit routine ..." << endl << endl;

	vector<double> Q_steps = load_Qsteps_data(FOLDER_IN_Q_PRC+"Qmin.dat");//load fom prc
	//Loop over all intervals: This can be naively parallelized!

	//#pragma omp parallel
	{
	//#pragma omp for
	for(size_t i = 0; i < Q_steps.size()-1; i++)
	{
		vector<double> Q_steps_i = load_Qsteps_data(FOLDER_IN_Q_SIM+to_string(i)+"_Qmin.dat");//load fom sim

		for(size_t j = 0; j < Q_steps_i.size()-1; j++)
		{
			cout << "Hi" << endl;
			vector<double> data;
			read_data(FOLDER_IN_DATA+to_string(i)+"_"+to_string(j)+FILE_NAME_IN, data);

			cout << "Hi" << endl;
			//read lambda vals for the j-th subinterval in the interval i
			vector<double> lambda_vals;
			cout << FOLDER_IN_DATA+to_string(i)+"_"+to_string(j)+FILE_NAME_lvals << endl;
			read_data(FOLDER_IN_DATA+to_string(i)+"_"+to_string(j)+FILE_NAME_lvals, lambda_vals); //+to_string(i)

			cout << "Hi" << endl;
			//Here we write the avg  (<<X>>_n + x_n)/delta + 0.5 and the JN error
			ofstream outfile1;
			outfile1.open(FOLDER_OUT + to_string(i)+"_"+to_string(j) + FILE_NAME_OUT1);

			cout << FOLDER_OUT + to_string(i)+"_"+to_string(j) + FILE_NAME_OUT1 << endl;


			//Here we write the avg + avg of JN samples to file!
			ofstream outfile2;
			outfile2.open(FOLDER_OUT + to_string(i)+"_"+to_string(j) + FILE_NAME_OUT2);

			double avg = 0.0;

			cout << "-->" << " " << lambda_vals.size() << " " << data.size() << " " << kmeas << endl;

			for(size_t k = 0; k < lambda_vals.size(); k++)
			{
				//get all measurements for one value of lambda
	//			vector<double> data_lam(data.begin()+k*nmeas, data.begin()+(k+1)*nmeas);
				vector<double> data_lam;

				for(size_t m = k*params::nmeas; m < (k+1)*params::nmeas; m++) data_lam.push_back(data.at(m));

				//get the mean wrt nmeas
				avg = get_avg(data_lam);

				//get the mean in the blocks
				vector<double> avg_blocks;

				for(int n = 0; n < nblocks; n++)
				{
	//				vector<double> data_lam_blocks(data_lam.begin()+n*kmeas, data_lam.begin()+(n+1)*kmeas);
					vector<double> data_lam_blocks;

					for(int m = n*kmeas; m < (n+1)*kmeas; m++) data_lam_blocks.push_back(data_lam.at(m));

					avg_blocks.push_back(get_avg(data_lam_blocks));
	//				cout << "######> " << data_lam.at(n) << endl;
				}

				//compute the Jackknife samples
				for(auto&& it: avg_blocks)
				{
					it = (params::nmeas*avg - kmeas*it)/((double)(params::nmeas - kmeas));
				}

				//compute the Jackknife error (for a specific interval and lambda value)
				double err = 0.0;
				for(auto it: avg_blocks) err += (avg - it)*(avg - it);
				err = (nblocks - 1.0)*err/nblocks;
				err = sqrt(err);

				//store the avg and Jackknife samples to file
	//			data_lam_jn_samples.at(k).at(0) = avg;
	//			copy(avg_blocks.begin(), avg_blocks.end(), data_lam_jn_samples.at(k).begin()+1);

				//write avg + JN errors to file
				outfile1 << lambda_vals.at(k) << " " << avg << " " << err << endl;

				//write avg + jn samples to file
				outfile2 << lambda_vals.at(k) << " " << avg << " ";
				for(size_t m = 0; m < avg_blocks.size() - 1; m++) outfile2 << avg_blocks.at(m) << " ";
				outfile2 << avg_blocks.back() << endl;
			}
			outfile1.close(); outfile2.close();

		}

	}
	}

	return 0;
}

void read_data(const string& FILE, vector<double>& data)
{
	ifstream infile_dat;
	infile_dat.open(FILE);
	if(infile_dat.is_open())
	{
	    double tmp;
	    while(!infile_dat.eof())
	    {
	        infile_dat >> tmp;
	        data.push_back(tmp);
	    }

	    data.pop_back();

	}
	else
	{
	    cout << "Unable to open file";
	    abort();
	}
	infile_dat.close();

}

double get_avg(const vector<double>& data)
{
	double avg = 0.0;
	for(auto it: data) avg += it;

	return avg/data.size();
}

//void write_data(ofstream& outfile_dat, const vector<double>& lambda_vals, const vector<double>& data)
//{
//	//Here we write the JN samples to file
//	ofstream outfile2;
//	outfile2.open(FOLDER + FILE_NAME_OUT2);
//
//	for(size_t i = 0; i < lambda_vals.size()-1; i++) outfile_dat << lambda_vals.at(i) << " ";
//	outfile_dat << outfile_dat.back() << endl;
//
//	for(size_t i = 0; i < data.size()-1; i++) outfile_dat << lambda_vals.at(i) << " ";
//	outfile_dat << f.back() << endl;
//
//	outfile2.close();
//}
