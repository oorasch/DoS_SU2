#ifndef INITIALIZE_HPP
#define INITIALIZE_HPP

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>
#include"parameters.hpp"

using namespace std;

void D_dim_neib_init(const int& d, const int &L, const int &Nt, vector<vector<int>>& neib);
vector<double> lambda_vals_init_prc();
vector<double> lambda_vals_init_sim(const double& lambda_mid, const double& delta_i, const double delta_ij);

#endif
