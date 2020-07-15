#include<iostream>
#include<fstream>
#include<cmath>
#include<algorithm>
#include<random>
#include<chrono>
#include<string>
#include<vector>
#include<iomanip>
//#include<eigen3/Eigen/Dense>
//#include<eigen3/Eigen/StdVector>
#include"../../eigen-3.3.7/Eigen/Dense"
#include"../../eigen-3.3.7/Eigen/StdVector"
#include<omp.h>
#include"io_data.hpp"
#include"initialize.hpp"
#include"observables.hpp"
#include"parameters.hpp"

using namespace std;
using namespace Eigen;

/*--------------------------------*/
/*-----------Prototypes---------- */
Matrix2cd get_SU2_rand();
void do_restricted_MC(const size_t& IntNum, const size_t& SubIntNum, const vector<double>& lambda_vals, const double& Qmin, const double& Qmax,
                      vector<int>* neib, vector<Matrix2cd,aligned_allocator<Matrix2cd> >* U_vec, double& Q);
void do_sweeps_pointer(const int& nsweeps, const double& Qmin, const double& Qmax, const double& lambda,
                       vector<int>* neib, vector<Matrix2cd,aligned_allocator<Matrix2cd> >* U_vec, double& Q);
void update_pointer(const double& Qmin, const double& Qmax, const double& lambda, vector<int>* neib,
                    vector<Matrix2cd,aligned_allocator<Matrix2cd> >* U_vec, double& Q);
/*--------------------------------*/
/*----Random Number Generator----*/

const unsigned seed(std::chrono::system_clock::now().time_since_epoch().count());
//std::default_random_engine generator(seed);
std::mt19937 generator(seed);

/*--------------------------------*/

string FOLDER_IN_PRC = "../prc/start_configs/Nt"+to_string(params::Nt)+"_Ns"+
                    to_string(params::Ns)+"_betastart"+to_string((int)(params::beta_start*10))+"/";
string FOLDER_IN_K = "../prc/stats/Nt"+to_string(params::Nt)+"_Ns"+
                    to_string(params::Ns)+"_beta"+to_string((int)(params::beta*10))+"/";
string FOLDER_IN = "../sim/start_configs/Nt"+to_string(params::Nt)+"_Ns"+
                    to_string(params::Ns)+"_betastart"+to_string((int)(params::beta_start*10))+"/";

string FOLDER_OUT = "../sim/data/Nt"+to_string(params::Nt)+"_Ns"+
                    to_string(params::Ns)+"_beta"+to_string((int)(params::beta*10))+"/";

int main(){

  vector<vector<Matrix2cd,aligned_allocator<Matrix2cd> > > U_vec (params::V, vector<Matrix2cd,aligned_allocator<Matrix2cd> >(params::dd, Matrix2cd::Identity()));
  vector<vector<int>> neib (params::V, vector<int>(2*params::dd, 0));
  D_dim_neib_init(params::dd, params::Ns, params::Nt, neib);
  vector<double> Q_steps = load_Qsteps_data(FOLDER_IN_PRC+"Qmin.dat");//load fom prc
  vector<double> lambda_mid = load_Qsteps_data(FOLDER_IN_K+"lambda_mid_vals.dat"); //load from prc

  //loop over Q_steps-1 intervals
  #pragma omp parallel for firstprivate(U_vec)
  for(size_t k = 0; k < Q_steps.size()-1; k++)
  {
    double Q = 0.0;

   // int n_delt = round((Q_steps.at(k+1) - Q_steps.at(k))/params::delta_sim);

    vector<double> Q_steps_i = load_Qsteps_data(FOLDER_IN+to_string(k)+"_Qmin.dat"); //load from sim
    vector<double> lambda_vals = lambda_vals_init_sim(lambda_mid.at(k), Q_steps.at(k+1) - Q_steps.at(k), params::delta_sim);//initialize lambda vals

    for(size_t i = 0; i < Q_steps_i.size()-1; i++)
    {
      //get configuration for this sub-interval
      string FILENAME = "start_config_Int"+to_string(k)+"_Subint"+to_string(i)+".dat";
      read_config(FOLDER_IN+FILENAME, U_vec);
      Q = meas_topo_clov(neib.data(),U_vec.data());

      cout << Q << endl;

      do_restricted_MC(k, i, lambda_vals, Q_steps_i.at(i), Q_steps_i.at(i+1), neib.data(), U_vec.data(), Q);

    }
  }

  return 0;
}

void do_restricted_MC(const size_t& IntNum, const size_t& SubIntNum, const vector<double>& lambda_vals, const double& Qmin, const double& Qmax, vector<int>* neib,
                      vector<Matrix2cd,aligned_allocator<Matrix2cd> >* U_vec,
                      double& Q)
{

  string FILENAME = to_string(IntNum)+"_"+to_string(SubIntNum)+"_lambda_vals.dat";
  ofstream outfile0(FOLDER_OUT+FILENAME, ios::out | ios::trunc);
  for(auto it: lambda_vals) outfile0 << it << endl;
  outfile0.close();


  //open file for restr. expectation values
  FILENAME = to_string(IntNum)+"_"+to_string(SubIntNum)+"_RE_Q_lambda.dat";
  ofstream outfile(FOLDER_OUT+FILENAME, ios::out | ios::trunc);

  //loop over lambda values
  for(size_t l = 0; l < lambda_vals.size(); l++)
  {
    //equilibrate the system in the interval with lambda
    auto start = chrono::high_resolution_clock::now();
    do_sweeps_pointer(params::nequi, Qmin, Qmax, lambda_vals.at(l), neib, U_vec, Q);
    auto end = chrono::high_resolution_clock::now();

    auto duration = chrono::duration<float, ratio<1>>(end-start);

    double rate = params::accepted/(1.0*params::tot_num_links);
    rate /= 1.0*params::nequi;
    rate /= 1.0*params::nmultihit;

    cout << endl << "Lambda is " << lambda_vals.at(l) << endl;
    cout << "Interval [" << Qmin << ", " << Qmax << "]" << endl;
    cout << "Equilibration took " << duration.count() << "s " << endl;
    cout << "Acceptance @ equilibration: " << rate << endl << endl;

    //check if topo is calculated correctly
    double Q_check = meas_topo_clov(neib, U_vec);
    if(fabs(Q/Q_check - 1.0) > 0.001)
    {
	    cout << "Something is wrong!" << endl;
	    cout << Q << " " << Q_check << endl;
	    abort();
    }

    vector<double> data;

    //loop over measurements
    for(int i = 0; i < params::nmeas; i++)
    {
      //decorrelation
      do_sweeps_pointer(params::nskip, Qmin, Qmax, lambda_vals.at(l), neib, U_vec, Q);

      //measure/write expectation vals
      data.push_back((Q - Qmin)/(Qmax-Qmin) - 0.5);

    }

    double mean = 0.0, err = 0.0;

    for(auto it: data)
    {
      mean += it;
      outfile << it << endl;
    }
    mean /= params::nmeas;

    for(auto it: data) err += (mean - it)*(mean - it);
    err = sqrt(err/(params::nmeas*(params::nmeas-1.0)));

    cout << Qmin << ", " << Qmax << ", " << lambda_vals.at(l) << " : " << mean << " +- " << err << endl;
    cout << "-------------------------------------------------------------------" << endl << endl;

  }
  outfile.close();
}

void do_sweeps_pointer(const int& nsweeps, const double& Qmin, const double& Qmax, const double& lambda,
                       vector<int>* neib, vector<Matrix2cd,aligned_allocator<Matrix2cd> >* U_vec, double& Q)
{
  params::accepted = 0;
	for(int i = 0; i < nsweeps; i++)
	{
		update_pointer(Qmin, Qmax, lambda, neib, U_vec, Q);
	}
}

void update_pointer(const double& Qmin, const double& Qmax, const double& lambda, vector<int>* neib,
                    vector<Matrix2cd,aligned_allocator<Matrix2cd> >* U_vec, double& Q)
{
	std::uniform_real_distribution<double> distribution(0,1);

	for(int site1 = 0; site1 < params::V; site1++)
	{
		for(int dir = 0; dir < params::dd; dir++)
		{
			//Check if we are even allowed to change the link
			if((*(neib+site1)).at(dir) != -1)
			{
				//get the staple for the link (site1, dir)
				Matrix2cd AA = Matrix2cd::Zero();
				for(int itdir = 0; itdir < params::dd; itdir++)
				{
					if(itdir != dir)
					{
						//check for boundary sites
						if((*(neib+site1)).at(itdir) != -1)
						{
							// AA += U_vec.at(neib.at(site1).at(dir)).at(itdir)* 																				//U_nu(n + \hat(mu))
							// 			U_vec.at(neib.at(site1).at(itdir)).at(dir).adjoint()* 															//U_-mu(n + \hat(mu) + \hat(nu))
							// 			U_vec.at(site1).at(itdir).adjoint(); 																								//U_-nu(n + \hat(nu))
							AA += (*(U_vec+(*(neib+site1)).at(dir))).at(itdir)*
										(*(U_vec+(*(neib+site1)).at(itdir))).at(dir).adjoint()*
										(*(U_vec+site1)).at(itdir).adjoint();
						}
						if((*(neib+site1)).at(itdir + params::dd) != -1)
						{
							// AA += U_vec.at(neib.at(neib.at(site1).at(dir)).at(itdir + dd)).at(itdir).adjoint()* 		//U_-nu(n+\hat(mu))
							// 			U_vec.at(neib.at(site1).at(itdir + dd)).at(dir).adjoint()* 												//U_-mu(n+\hat(mu) - \hat(nu))
							// 			U_vec.at(neib.at(site1).at(itdir + dd)).at(itdir); 																//U_nu(n - \hat(nu))
							AA += (*(U_vec+(*(neib+(*(neib+site1)).at(dir))).at(itdir + params::dd))).at(itdir).adjoint()*
										(*(U_vec+(*(neib+site1)).at(itdir + params::dd))).at(dir).adjoint()*
										(*(U_vec+(*(neib+site1)).at(itdir + params::dd))).at(itdir);
						}
					}
				}

        double Q_new = 0.0, Q_c_new = 0.0, Q_c_old = 0.0;

	//do multihit MC
	for(int k = 0; k < params::nmultihit; k++)
	{
		//get random matrix
		Matrix2cd X;
		{
			Matrix2cd tmp = get_SU2_rand();
			//We pick X and X ^ dagger with probability 1/2
			if(distribution(generator) < 0.5) X = tmp;
				else X = tmp.adjoint();
		}

 	  //get local Q contribution Q_c_old; we calculate this once then Q_c_old^(i+1) = Q_c_new^(i) if move is accepted
    if(k == 0) Q_c_old = get_local_topo_contribution_clov(site1, dir, neib, U_vec);

		//here we change the link variable - this needs to be undone if a move is not allowed!!!
		Matrix2cd U_store = (*(U_vec+site1)).at(dir);
		(*(U_vec+site1)).at(dir) = X*U_store;

		//get local Q contribution Q_c_new
		Q_c_new = get_local_topo_contribution_clov(site1, dir, neib, U_vec);
		Q_new = Q + Q_c_new - Q_c_old;

		if(Q_new >= Qmin && Q_new <= Qmax)
		{
      //build S_loc
  		Matrix2cd aux = ((*(U_vec+site1)).at(dir) - U_store)*AA;
  		// double rho = params::beta*real(aux.trace())/2.0 + lambda*(Q_c_new - Q_c_old);
      double rho = params::beta*real(aux.trace())/2.0 + lambda*(Q_new - Q);

			if(rho < 0.0)
			{
				if(distribution(generator) <= exp(rho))
				{
					params::accepted++;
				  Q = Q_new;
					Q_c_old = Q_c_new;
				}
				else (*(U_vec+site1)).at(dir) = U_store;//undo change
			}
			else //always accept
			{
				params::accepted++;
				Q = Q_new;
				Q_c_old = Q_c_new;
			}
		}
		else (*(U_vec+site1)).at(dir) = U_store;//undo change
    }
   }
  }
  // cout << site1 << " " << Q << " " << meas_topo_plaq(neib, U_vec) << endl;
 }
}

Matrix2cd get_SU2_rand()
{
  uniform_real_distribution<double> rand1 (-0.5, 0.5);
  double r0 = rand1(generator),
         r1 = rand1(generator),
         r2 = rand1(generator),
         r3 = rand1(generator);

  double norm = sqrt(r1*r1 + r2*r2 + r3*r3);

  r0 = (r0/fabs(r0))*sqrt(1.0 - params::eps*params::eps);
  r1 *= params::eps/norm;
  r2 *= params::eps/norm;
  r3 *= params::eps/norm;

  Matrix2cd mat;
  mat << complex<double> (r0, r3), complex<double> (r2, r1),
         complex<double> (-r2, r1), complex<double> (r0, -r3);

  return mat;
}
