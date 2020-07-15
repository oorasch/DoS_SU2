#include "initialize.hpp"

void D_dim_neib_init(const int& d, const int &L, const int &Nt, vector<vector<int>>& neib)
{
	int V = pow(L, d-1)*Nt;

	for(int k = 0; k < V; k++)
	{
		vector<int> k_coords;
		int x = k%L;
		int k_slice = k;
		int tmp_k = 0, Ln = 0;
//		cout << k << " " << x << " ";
		k_coords.push_back(x);

		for(int i = 1; i < d; i++)
		{
			k_slice = (k_slice - x)/L;
			x = k_slice%L;
//			cout << x << " ";
			k_coords.push_back(x);

		}

//		cout << endl;

		//positive time direction
		int k_coords_save = k_coords.at(d-1);
		k_coords.at(d-1) += 1;

		if(k_coords.at(d-1) < Nt)
		{
			tmp_k = 0, Ln = 1;

			for(auto it: k_coords)
			{
				tmp_k += Ln*it;
				Ln *= L;
			}

			neib.at(k).at(0) = tmp_k;
		}
		else if(k_coords.at(d-1) == Nt) neib.at(k).at(0) = -1;

		k_coords.at(d-1) = k_coords_save;

		// if(k_coords.at(d-1) == Nt) k_coords.at(d-1) = 0;
		//
		// int tmp_k = 0, Ln = 1;
		//
		// for(auto it: k_coords)
		// {
		// 	tmp_k += Ln*it;
		// 	Ln *= L;
		// }
		//
		// neib.at(k).at(0) = tmp_k;
		// k_coords.at(d-1) = k_coords_save;

		//spatial neibs in positive direction
		for(int i = 0; i < d-1; i++)
		{
			int k_coords_save = k_coords.at(i);
			k_coords.at(i) += 1;

			if(k_coords.at(i) < L)
			{
				tmp_k = 0, Ln = 1;

				for(auto it: k_coords)
				{
					tmp_k += Ln*it;
					Ln *= L;
				}

				neib.at(k).at(1 + i) = tmp_k;
			}
			else if(k_coords.at(i) == L) neib.at(k).at(1 + i) = -1;

			k_coords.at(i) = k_coords_save;
		}

		//negative time direction
		k_coords_save = k_coords.at(d-1);
		k_coords.at(d-1) -= 1;

		if(k_coords.at(d-1) >= 0)
		{
			tmp_k = 0, Ln = 1;

			for(auto it: k_coords)
			{
				tmp_k += Ln*it;
				Ln *= L;
			}

			neib.at(k).at(d) = tmp_k;
		}
		else if(k_coords.at(d-1) < 0) neib.at(k).at(d) = -1;

		k_coords.at(d-1) = k_coords_save;

		//neibs in negative spatial direction
		for(int i = 0; i < d-1; i++)
		{
			int k_coords_save = k_coords.at(i);
			k_coords.at(i) -= 1;

			if(k_coords.at(i) >= 0)
			{
				int tmp_k = 0, Ln = 1;

				for(auto it: k_coords)
				{
					tmp_k += Ln*it;
					Ln *= L;
				}

				neib.at(k).at(d + 1 + i) = tmp_k;
			}
			else if(k_coords.at(i) < 0) neib.at(k).at(d + 1 + i) = -1;

			k_coords.at(i) = k_coords_save;
		}

	}
}

//we set the default value of lambda mid = 0
//in case we do preconditioning, we can adjust the mid point by the approximated value
vector<double> lambda_vals_init_prc()
{
	vector<double> lambda_vals;

 	double tmp = -params::lambda_start_prc;
 	for(int i = 0; i <= params::n_lambda_prc; i++)
 	{
 		lambda_vals.push_back(tmp);
 		tmp += params::lambda_inc_prc;
 	}

	return lambda_vals;
}

vector<double> lambda_vals_init_sim(const double& lambda_mid, const double& delta_i, const double delta_ij)
{
	vector<double> lambda_vals;

	double lambda_start = lambda_mid - delta_i*(params::lambda_inc_prc)/delta_ij;
	double lambda_end = lambda_mid + delta_i*(params::lambda_inc_prc)/delta_ij;
	double lambda_width = lambda_end - lambda_start;
	double inc = lambda_width/params::n_lambda_sim;

	for(int i = 0; i <= params::n_lambda_sim; i++)
	{
		lambda_vals.push_back(lambda_start + inc*i);
	}

	return lambda_vals;
}
