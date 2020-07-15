#include "io_data.hpp"

vector<double> load_Qsteps_data(const string& FILENAME)
{
  vector<double> Q_steps;
  ifstream infile (FILENAME);
  if (infile.is_open())
  {
    double line;
    while ( infile >> line )
    {
      Q_steps.push_back(line);
    }
    infile.close();
  }
  else cout << "Unable to open file";

  return Q_steps;
}

void read_config(const string& FILENAME, vector<vector<Matrix2cd,aligned_allocator<Matrix2cd> > >& U_vec)
{
	ifstream infile(FILENAME);
	string line, cell;

	if(infile)
	{
		for(int sites = 0; sites < params::V; sites++)
		{
			for(int dir = 0; dir < params::dd; dir++)
			{
				//We read the lines from file, if len(File) != V*d then the prog throws an exception
				if(!getline(infile, line)) cout << "Something is wrong!! (" << sites << "," << dir << ")" << endl;
		    std::stringstream lineStream(line);

				int i = 0;
		    while (getline(lineStream, cell, ','))
		    {
					double elem = stod(cell, nullptr);
					//We read the (0,0) and (0,1) element of the SU2 matrices from the file
					//Then we build the FULL matrices according to
					//M(1,1) = MM(0,0)^*, M(1,0) = -M(0,1)^*
					switch(i)
					{
						case 0:
							U_vec.at(sites).at(dir)(0,0).real(elem);
							U_vec.at(sites).at(dir)(1,1).real(elem);
							break;
						case 1:
							U_vec.at(sites).at(dir)(0,0).imag(elem);
							U_vec.at(sites).at(dir)(1,1).imag(-elem);
							break;
						case 2:
							U_vec.at(sites).at(dir)(0,1).real(elem);
							U_vec.at(sites).at(dir)(1,0).real(-elem);
							break;
						case 3:
							U_vec.at(sites).at(dir)(0,1).imag(elem);
							U_vec.at(sites).at(dir)(1,0).imag(elem);
							break;
						case 4: cout << "Somethings wrong!!!" << endl; abort();
					}
					i++;
		    }
			}
		}
	}
	else
	{
	  cerr << "Could not open file!" << endl;
		abort();
	}
}

void write_config(const string& FILENAME, vector<vector<Matrix2cd,aligned_allocator<Matrix2cd> > >& U_vec)
{
	ofstream outfile(FILENAME, ios::out | ios::trunc);

	if(outfile)
	{
		for(auto it_site: U_vec)
		{
			for(auto it_link: it_site)
			{
				outfile << std::setprecision(15) << it_link(0,0).real() <<","<< it_link(0,0).imag()<<",";
				outfile << std::setprecision(15) << it_link(0,1).real() <<","<< it_link(0,1).imag()<<endl;
			}
		}
	  outfile.close();
	}
	else
	{
		cerr << "Could not open file!" << endl;
	}
}
