#ifndef OBSERVABLES_HPP
#define OBSERVABLES_HPP

#include<iostream>
#include<vector>
//#include<eigen3/Eigen/Dense>
//#include<eigen3/Eigen/StdVector>
#include"../../eigen-3.3.7/Eigen/Dense"
#include"../../eigen-3.3.7/Eigen/StdVector"
#include"parameters.hpp"

using namespace std;
using namespace Eigen;

double meas_topo_clov(vector<int>* neib, vector<Matrix2cd,aligned_allocator<Matrix2cd> >* U_vec);
double get_local_topo_contribution_clov(const int& site, const int& dir, vector<int>* neib,
                                   vector<Matrix2cd,aligned_allocator<Matrix2cd> >* U_vec);
double meas_topo_plaq(vector<int>* neib, vector<Matrix2cd,aligned_allocator<Matrix2cd> >* U_vec);
double get_local_topo_contribution_plaq(const int& site, const int& dir, vector<int>* neib,
                                   vector<Matrix2cd,aligned_allocator<Matrix2cd> >* U_vec);


#endif
