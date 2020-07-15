#ifndef IO_DATA_HPP
#define IO_DATA_HPP

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<iomanip>
//#include<eigen3/Eigen/Dense>
//#include<eigen3/Eigen/StdVector>
#include"../../eigen-3.3.7/Eigen/Dense"
#include"../../eigen-3.3.7/Eigen/StdVector"
#include"parameters.hpp"

using namespace std;
using namespace Eigen;

vector<double> load_Qsteps_data(const string& FILENAME);
void read_config(const string& FILENAME, vector<vector<Matrix2cd,aligned_allocator<Matrix2cd> > >& U_vec);
void write_config(const string& FILENAME, vector<vector<Matrix2cd,aligned_allocator<Matrix2cd> > >& U_vec);

#endif
