#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include<string>
#include<cmath>
#include<vector>

namespace params
{
  extern const std::string FOLDER;
  extern const int dd;
  extern const int Nt;
  extern const int Ns;
  extern const int V;
  extern const double beta;
  extern const double beta_start;
  extern const int nmeas;
  extern const int nequi;
  extern const int nskip;
  extern const int nmultihit;
  extern const double eps;
  extern long long int accepted;
  extern const long long int tot_num_links;
  extern const double pi;
  extern const std::vector<std::vector<int>> perms_red;
  extern const std::vector<int> signum_red;
  extern const int n_lambda_prc;
  extern const double lambda_start_prc;
  extern const double lambda_inc_prc;
  extern const double delta_sim;
  extern const int n_lambda_sim;
}

#endif
