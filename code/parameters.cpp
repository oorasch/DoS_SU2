#include "parameters.hpp"

namespace params
{
  extern const std::string FOLDER("start_configs/");

  /*Lattice parameters*/
  extern const int dd(4);
  extern const int Nt(8);
  extern const int Ns(8);
  extern const int V = Nt*pow(Ns, dd-1);

  /*Simulation parameters*/
  extern const double beta(2.0);
  extern const double beta_start(1.0);
  extern const int nmeas(4e3);
  extern const int nequi(1e3);
  extern const int nskip(1);
  extern const int nmultihit(4);
  extern const double eps(0.5);
  extern long long int accepted(0);
  extern const long long int tot_num_links = 4*(Nt-1)*(Ns-1)*(Ns-1)*(Ns-1) + 9*(Nt-1)*(Ns-1)*(Ns-1)+
                                             3*(Ns-1)*(Ns-1)*(Ns-1) + 6*(Nt-1)*(Ns-1) + 6*(Ns-1)*(Ns-1)+
                                             3*(Ns-1) +(Nt-1); //# of links of OBC lattice


  /*Observables*/
  extern const double pi(3.141592653589793);
  extern const std::vector<std::vector<int>> perms_red = {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2}};
  extern const std::vector<int> signum_red = {1, -1, 1};

  /*Parameters for prc*/
  extern const int n_lambda_prc = 20;
  extern const double lambda_start_prc = 100.0;
  extern const double lambda_inc_prc = lambda_start_prc/n_lambda_prc;

  /*Parameters for sim*/
  extern const double delta_sim = 0.025;
  extern const int n_lambda_sim = 40;
}
