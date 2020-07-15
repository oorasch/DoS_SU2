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
#include"../../eigen-3.3.7/Eigen/Dense"
//#include<eigen3/Eigen/StdVector>
#include"../../eigen-3.3.7/Eigen/StdVector"
#include<omp.h>
//#include "../LibEigen/Eigen/Dense"
#include "parameters.hpp"
#include "initialize.hpp"
#include "io_data.hpp"
#include "observables.hpp"

using namespace std;
using namespace Eigen;
using namespace std::chrono;

/*---------Local parameters--------*/

const int nmeas = 1e4;
const int ndecorr = 10;
const int nskip = 1;
const int nmultihit = 2;
const double eps = 0.01;

long long int accepted = 0;
const double delta_fine = 0.1;
const double delta_coarse = 0.2;
const int Int_refine = 4;
string F_params = "Nt"+to_string(params::Nt)+"_Ns"+to_string(params::Ns)+"_betastart"+to_string((int)(params::beta_start*10))+"/";
const string FOLDER_IN = "../prc/start_configs/"+F_params;
const string FOLDER_OUT = "../sim/start_configs/"+F_params;

/*--------------------------------*/
/*-----------Prototypes---------- */
Matrix2cd get_SU2_rand();
void do_sweeps_pointer(const double& Qmin, const double& Qmax, const int& nsweeps, vector<int>* neib, vector<Matrix2cd,aligned_allocator<Matrix2cd> >* U_vec, double& Q);
void update_pointer(const double& Qmin, const double& Qmax, vector<int>* neib, vector<Matrix2cd,aligned_allocator<Matrix2cd> >* U_vec, double& Q);
void gaugeT_test(vector<int>* neib, vector<vector<Matrix2cd,aligned_allocator<Matrix2cd> > > U_vec);
vector<double> fill_Qmin_vec(const int& i, const double& delta_prc, const double& delta_pr, const double& Qmin);
/*--------------------------------*/
/*----Random Number Generator----*/

const unsigned seed(std::chrono::system_clock::now().time_since_epoch().count());
//std::default_random_engine generator(seed);
std::mt19937 generator(seed);

/*--------------------------------*/

int main(){

  double Q = 0.0; // current value of top. charge
	vector<vector<int>> neib (params::V, vector<int>(2*params::dd, 0));
  string FILENAME = "Qmin.dat";
  vector<double> Qmin_vec_prc = load_Qsteps_data(FOLDER_IN+FILENAME);

	D_dim_neib_init(params::dd, params::Ns, params::Nt, neib);

	// #pragma omp parallel
	for(size_t i = 0; i < Qmin_vec_prc.size()-1; i++)
  {
    int IntsFull = -1;
    double Qmin = Qmin_vec_prc.at(i), Qmax = Qmin_vec_prc.at(i+1);
    int nSubInt = round((Qmax-Qmin)/params::delta_sim);
    vector<bool> Intervals_filled (nSubInt, false);

    cout << "---------------> " << Qmin << " " << Qmax << " " << nSubInt << endl;

    //read U_vec for the interval i
    vector<vector<Matrix2cd, aligned_allocator<Matrix2cd> > > U_vec (params::V, vector<Matrix2cd,aligned_allocator<Matrix2cd> >(params::dd, Matrix2cd::Identity()));
    FILENAME = "start_config_Int"+to_string(i)+".dat";
    read_config(FOLDER_IN+FILENAME, U_vec);
    Q = meas_topo_clov(neib.data(), U_vec.data());
    cout << Q << endl;

  	do_sweeps_pointer(Qmin, Qmax, ndecorr, neib.data(), U_vec.data(), Q);

    Q = meas_topo_clov(neib.data(), U_vec.data());
    cout << Q << endl;

    if(!(Qmin <= Q && Q <= Qmax))
    {
      cout << Q << " " << Qmin << " " << Qmax << endl;
      cout << "Q[U] not in interval ... abort!" << endl;
      abort();
    }

    vector<double> Qmin_vec_i = fill_Qmin_vec(i, Qmax - Qmin, params::delta_sim, Qmin);

  	// #pragma omp parallel for firstprivate(Q, U_vec)
  	for(int n = 0; n < nmeas; n++)
  	{
  		if(IntsFull >= 0) continue; //Once all Intervals are filled we skip the rest of the measurements.
  		//This is because OMP does not allow BREAK for OMP loops

  		do_sweeps_pointer(Qmin, Qmax, nskip, neib.data(), U_vec.data(), Q);

      cout << "Acceptance: " << accepted/((double)params::tot_num_links*nskip*nmultihit*4.0) << endl;

  		// #pragma omp critical
  		for(size_t k = 0; k < Qmin_vec_i.size()-1; k++)
  		{
  			if(Qmin_vec_i.at(k) <= Q && Q <= Qmin_vec_i.at(k+1) && !Intervals_filled.at(k))
  			{
  				cout << "Interval " << k << " full!" << endl;
          cout << Q << endl;

  				string FILENAME = "start_config_Int"+to_string(i)+"_Subint"+to_string(k)+".dat";
  				write_config(FOLDER_OUT+FILENAME, U_vec);

  				Intervals_filled.at(k) = true;
  			}
  		}

  		{
  			bool tmp = true;
  			for(auto it: Intervals_filled)
  			{
  				if(!it)
  				{
  					tmp = false;
  					break;
  				}
  			}

  			if(tmp) IntsFull = n;

  		}
  	}

    //-----------------------------------------------------------------------
    //--------------------Check if all bins have been filled-----------------
    {
      bool tmp = true;
      for(auto it: Intervals_filled) if(!it) tmp = false;

      if(tmp) cout << "All bins have been filled!" << endl;
      else{
        for(size_t k = 0; k < Intervals_filled.size(); k++)
        {
          if(!Intervals_filled.at(k)) cout << "Bin " << k << "has not been filled!" << endl;
        }
        cout << "Increase statistics!" << endl;
      }
    }
  }

	//-----------------------------------------------------------------------
	//------------Check if gauge invariance has not been violated------------
	// if(omp_get_thread_num() == 0)
	//#pragma omp critical
	// {
	// 	cout << meas_topo_clov(neib.data(), U_vec.data()) << endl;
	// 	gaugeT_test(neib.data(), U_vec);
	// }
  //-----------------------------------------------------------------------
  //--------------------Check local definition of Q------------------------

  // int kk = 0.5*params::Ns + 0.5*params::Ns*params::Ns + 0.5*params::Ns*params::Ns*params::Ns + 0.5*params::Nt*params::Ns*params::Ns*params::Ns;
  //
  // Q = meas_topo_clov(neib.data(), U_vec.data());
  // double Q_c_old = get_local_topo_contribution(kk, 0, neib.data(), U_vec.data());
  //
  // U_vec.at(kk).at(0) = get_SU2_rand()*U_vec.at(kk).at(0);
  //
  // double Q_c_new = get_local_topo_contribution(kk, 0, neib.data(), U_vec.data());
  //
  // cout << endl << Q << endl;
  // cout << Q + Q_c_new - Q_c_old << endl;
  // cout << meas_topo_clov(neib.data(), U_vec.data()) << endl;

  return 0;
}

Matrix2cd get_SU2_rand()
{
  uniform_real_distribution<double> rand1 (-0.5, 0.5);
  double r0 = rand1(generator),
         r1 = rand1(generator),
         r2 = rand1(generator),
         r3 = rand1(generator);

  double norm = sqrt(r1*r1 + r2*r2 + r3*r3);

  r0 = (r0/fabs(r0))*sqrt(1.0 - eps*eps);
  r1 *= eps/norm;
  r2 *= eps/norm;
  r3 *= eps/norm;

  Matrix2cd mat;
  mat << complex<double> (r0, r3), complex<double> (r2, r1),
         complex<double> (-r2, r1), complex<double> (r0, -r3);

  return mat;
}

void do_sweeps_pointer(const double& Qmin, const double& Qmax, const int& nsweeps, vector<int>* neib, vector<Matrix2cd,aligned_allocator<Matrix2cd> >* U_vec, double& Q)
{
	accepted = 0;
	for(int i = 0; i < nsweeps; i++)
	{
		update_pointer(Qmin, Qmax, neib, U_vec, Q);
	}
}

void update_pointer(const double& Qmin, const double& Qmax, vector<int>* neib, vector<Matrix2cd,aligned_allocator<Matrix2cd> >* U_vec, double& Q)
{
	std::uniform_real_distribution<double> distribution(0,1);
	std::uniform_int_distribution<int> site_dist(0,params::V-1);
	std::uniform_int_distribution<int> dir_dist(0,params::dd-1);

	for(int site1 = 0; site1 < params::V; site1++)
	{
		for(int dir = 0; dir < params::dd; dir++)
		{
			// int site1 = site_dist(generator), dir = dir_dist(generator);
	//		int site1 = 0, dir = 0;
			// cout << site1 << " " << dir << endl;

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
				for(int k = 0; k < nmultihit; k++)
				{
					//get random matrix
					Matrix2cd X;
					{
						Matrix2cd tmp = get_SU2_rand();
						//We pick X and X ^ dagger with probability 1/2
						if(distribution(generator) < 0.5) X = tmp;
						else X = tmp.adjoint();
					}

          Q_c_old = get_local_topo_contribution_clov(site1, dir, neib, U_vec);

          Matrix2cd U_store = (*(U_vec+site1)).at(dir);
          (*(U_vec+site1)).at(dir) = X*U_store;

          Q_c_new = get_local_topo_contribution_clov(site1, dir, neib, U_vec);

          Q_new = Q + Q_c_new - Q_c_old;

          if(Qmin <= Q_new && Q_new <= Qmax)
          {
    				//build S_loc
    				Matrix2cd aux = ((*(U_vec+site1)).at(dir) - U_store)*AA;
    				double rho = (-1.0)*params::beta_start*real(aux.trace())/2.0;

  					if(rho > 0.0)
  					{
  						if(distribution(generator) <= exp(-rho))
  						{
  							accepted++;
                Q = Q_new;
  						}
  						else (*(U_vec+site1)).at(dir) = U_store;
  					}
  					else //always accept
  					{
  						accepted++;
              Q = Q_new;
  					}
          }
          else (*(U_vec+site1)).at(dir) = U_store;
				}
			}
		}
	}
}

void gaugeT_test(vector<int>* neib, vector<vector<Matrix2cd,aligned_allocator<Matrix2cd> > > U_vec)
{
	vector<Matrix2cd,aligned_allocator<Matrix2cd> > Lambda;

	cout << endl << "Gauge transformation test:" << endl << endl;

	//Here we generate a set of gauge transformations
	for(int i = 0; i < params::V; i++)
	{
		Lambda.push_back(get_SU2_rand());
	}

	double qq = 0.0, qq_prime = 0.0;
	qq = meas_topo_clov(neib, U_vec.data());

	for(int i = 0; i < params::V; i++)
	{
		for(int dir = 0; dir < params::dd; dir++)
		{
			if((*(neib+i)).at(dir) != -1)
			{
				U_vec.at(i).at(dir) = Lambda.at(i)*(U_vec.at(i).at(dir)*Lambda.at((*(neib+i)).at(dir)).adjoint());
			}
		}
	}

	qq_prime = meas_topo_clov(neib, U_vec.data());
	cout << "Q[U]/Q[U'] = " << qq/qq_prime << "(" << qq << "," << qq_prime << ")" << endl;

}

vector<double> fill_Qmin_vec(const int& i, const double& delta_prc, const double& delta_sim, const double& Qmin)
{
	vector<double> Qmin_vec;
	double tmp = Qmin;
  int n_delt = round(delta_prc/delta_sim);

	for(int i = 0; i <= n_delt; i++)
	{
		Qmin_vec.push_back(tmp);
		tmp += delta_sim;
	}

	//write the vector to file
	string FILENAME = to_string(i)+"_Qmin.dat";
	// cout << FOLDER+FILENAME << endl;
	ofstream outfile(FOLDER_OUT+FILENAME, ios::out | ios::trunc);
	for(auto it: Qmin_vec) outfile << it << endl;
	outfile.close();

	return Qmin_vec;
}
