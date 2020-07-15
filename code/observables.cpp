#include "observables.hpp"

double meas_topo_clov(vector<int>* neib, vector<Matrix2cd,aligned_allocator<Matrix2cd> >* U_vec)
{
	double qq = 0.0;
	Matrix2cd C_munu = Matrix2cd::Zero(), C_rhosig = Matrix2cd::Zero();
	int mu, nu, rho, sig;
	int xpmu, xpnu, xpnummu, xmmu, xmnummu, xmnu, xmnupmu;
	int xprho, xpsig, xpsigmrho, xmrho, xmsigmrho, xmsig, xmsigprho;

	const vector<vector<int>> perms_red  = {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2}};
	const vector<int> signum_red = {1, -1, 1};

	for(int t = 1; t < params::Nt-1; t++)
	{
		for(int z = 1; z < params::Ns-1; z++)
		{
			for(int y = 1; y < params::Ns-1; y++)
			{
				for(int x = 1; x < params::Ns-1; x++)
				{
					int n = x + y*params::Ns + z*params::Ns*params::Ns + t*params::Ns*params::Ns*params::Ns;

					for(int k = 0; k < 3; k++)
					{
						mu = params::perms_red.at(k).at(0);
						nu = params::perms_red.at(k).at(1);

						xpmu = (*(neib+n)).at(mu);
						xmnu = (*(neib+n)).at(nu + params::dd);
						xpnu = (*(neib+n)).at(nu);
						xmmu = (*(neib+n)).at(mu + params::dd);

						xpnummu = (*(neib+xpnu)).at(mu + params::dd);
						xmnummu = (*(neib+xmnu)).at(mu + params::dd);
						xmnupmu = (*(neib+xmnu)).at(mu);

						rho = params::perms_red.at(k).at(2);
						sig = params::perms_red.at(k).at(3);

						xprho = (*(neib+n)).at(rho);
						xmrho = (*(neib+n)).at(rho + params::dd);
						xpsig = (*(neib+n)).at(sig);
						xmsig = (*(neib+n)).at(sig + params::dd);

						xpsigmrho = (*(neib+xpsig)).at(rho + params::dd);
						xmsigmrho = (*(neib+xmsig)).at(rho + params::dd);
						xmsigprho = (*(neib+xmsig)).at(rho);

						// 		C_munu  = U_vec.at(n).at(mu)*U_vec.at(xpmu).at(nu)*U_vec.at(xpnu).at(mu).adjoint()*U_vec.at(n).at(nu).adjoint();
						// 		C_munu -= U_vec.at(n).at(mu)*U_vec.at(xmnupmu).at(nu).adjoint()*U_vec.at(xmnu).at(mu).adjoint()*U_vec.at(xmnu).at(nu);
						// 		C_munu += U_vec.at(xmmu).at(mu).adjoint()*U_vec.at(xmnummu).at(nu).adjoint()*U_vec.at(xmnummu).at(mu)*U_vec.at(xmnu).at(nu);
						// 		C_munu -= U_vec.at(xmmu).at(mu).adjoint()*U_vec.at(xmmu).at(nu)*U_vec.at(xpnummu).at(mu)*U_vec.at(n).at(nu).adjoint();

						C_munu  = (*(U_vec+n)).at(mu)*(*(U_vec+xpmu)).at(nu)*(*(U_vec+xpnu)).at(mu).adjoint()*(*(U_vec+n)).at(nu).adjoint();
						C_munu -= (*(U_vec+n)).at(mu)*(*(U_vec+xmnupmu)).at(nu).adjoint()*(*(U_vec+xmnu)).at(mu).adjoint()*(*(U_vec+xmnu)).at(nu);
						C_munu += (*(U_vec+xmmu)).at(mu).adjoint()*(*(U_vec+xmnummu)).at(nu).adjoint()*(*(U_vec+xmnummu)).at(mu)*(*(U_vec+xmnu)).at(nu);
						C_munu -= (*(U_vec+xmmu)).at(mu).adjoint()*(*(U_vec+xmmu)).at(nu)*(*(U_vec+xpnummu)).at(mu)*(*(U_vec+n)).at(nu).adjoint();
				//
						{
							Matrix2cd tmp = C_munu.adjoint();
							C_munu -= tmp;
						}
				//
				// 		C_rhosig  = U_vec.at(n).at(rho)*U_vec.at(xprho).at(sig)*U_vec.at(xpsig).at(rho).adjoint()*U_vec.at(n).at(sig).adjoint();
				// 		C_rhosig -= U_vec.at(n).at(rho)*U_vec.at(xmsigprho).at(sig).adjoint()*U_vec.at(xmsig).at(rho).adjoint()*U_vec.at(xmsig).at(sig);
				// 		C_rhosig += U_vec.at(xmrho).at(rho).adjoint()*U_vec.at(xmsigmrho).at(sig).adjoint()*U_vec.at(xmsigmrho).at(rho)*U_vec.at(xmsig).at(sig);
				// 		C_rhosig -= U_vec.at(xmrho).at(rho).adjoint()*U_vec.at(xmrho).at(sig)*U_vec.at(xpsigmrho).at(rho)*U_vec.at(n).at(sig).adjoint();

						C_rhosig  = (*(U_vec+n)).at(rho)*(*(U_vec+xprho)).at(sig)*(*(U_vec+xpsig)).at(rho).adjoint()*(*(U_vec+n)).at(sig).adjoint();
						C_rhosig -= (*(U_vec+n)).at(rho)*(*(U_vec+xmsigprho)).at(sig).adjoint()*(*(U_vec+xmsig)).at(rho).adjoint()*(*(U_vec+xmsig)).at(sig);
						C_rhosig += (*(U_vec+xmrho)).at(rho).adjoint()*(*(U_vec+xmsigmrho)).at(sig).adjoint()*(*(U_vec+xmsigmrho)).at(rho)*(*(U_vec+xmsig)).at(sig);
						C_rhosig -= (*(U_vec+xmrho)).at(rho).adjoint()*(*(U_vec+xmrho)).at(sig)*(*(U_vec+xpsigmrho)).at(rho)*(*(U_vec+n)).at(sig).adjoint();

						{
							Matrix2cd tmp = C_rhosig.adjoint();
							C_rhosig -= tmp;
						}

						// cout << (C_munu*C_rhosig).trace() << endl;

						qq += params::signum_red.at(k)*real((C_munu*C_rhosig).trace());
					}
				}
			}
		}
	}

	return -qq/256.0/params::pi/params::pi;
}

double get_local_topo_contribution_clov(const int& site, const int& dir, vector<int>* neib, vector<Matrix2cd,aligned_allocator<Matrix2cd> >* U_vec)
{
	vector<int> dir_perp;
	double qq = 0.0;
	Matrix2cd C_munu = Matrix2cd::Zero(), C_rhosig = Matrix2cd::Zero();
	int mu, nu, rho, sig;
	int xpmu, xpnu, xpnummu, xmmu, xmnummu, xmnu, xmnupmu;
	int xprho, xpsig, xpsigmrho, xmrho, xmsigmrho, xmsig, xmsigprho;

	//get all direction (positive) perpendicular to dir
	switch(dir)
	{
		case 0:
			dir_perp = {1,2,3};
			break;
		case 1:
			dir_perp = {0,2,3};
			break;
		case 2:
			dir_perp = {0,1,3};
			break;
		case 3:
			dir_perp = {0,1,2};
			break;
		default:
			abort();
	}

	//get all the base points for the affected clover leaves
  int xneib = (*(neib+site)).at(dir);
	vector<int> site_vec = {site, xneib};
	for(auto it_dir: dir_perp)
	{
		site_vec.push_back((*(neib+site)).at(it_dir));
		site_vec.push_back((*(neib+site)).at(it_dir + params::dd));

		if(xneib >= 0)
	 	{
			site_vec.push_back((*(neib+xneib)).at(it_dir));
			site_vec.push_back((*(neib+xneib)).at(it_dir + params::dd));
		}
	}

	// for(auto it: site_vec) cout << it << endl;

	for(auto it_site: site_vec)
	{
		int n = it_site;

		if(n < 0) continue;

		for(int k = 0; k < 3; k++)
		{

			mu = params::perms_red.at(k).at(0);
			nu = params::perms_red.at(k).at(1);

			xpmu = (*(neib+n)).at(mu);
			xpnu = (*(neib+n)).at(nu);
			xmmu = (*(neib+n)).at(mu + params::dd);
			xmnu = (*(neib+n)).at(nu + params::dd);

			if(xpmu < 0 || xpnu < 0 || xmmu < 0 || xmnu < 0) continue;

			xpnummu = (*(neib+xpnu)).at(mu + params::dd);
			xmnummu = (*(neib+xmnu)).at(mu + params::dd);
			xmnupmu = (*(neib+xmnu)).at(mu);

			if(xpnummu < 0 || xmnummu < 0 || xmnupmu < 0) continue;

			rho = params::perms_red.at(k).at(2);
			sig = params::perms_red.at(k).at(3);

			xprho = (*(neib+n)).at(rho);
			xpsig = (*(neib+n)).at(sig);
			xmrho = (*(neib+n)).at(rho + params::dd);
			xmsig = (*(neib+n)).at(sig + params::dd);

			if(xprho < 0 || xpsig < 0 || xmrho < 0 || xmsig < 0) continue;

			xpsigmrho = (*(neib+xpsig)).at(rho + params::dd);
			xmsigmrho = (*(neib+xmsig)).at(rho + params::dd);
			xmsigprho = (*(neib+xmsig)).at(rho);

			if(xpsigmrho < 0 || xmsigmrho < 0 || xmsigprho < 0) continue;
	//
	// 		C_munu  = U_vec.at(n).at(mu)*U_vec.at(xpmu).at(nu)*U_vec.at(xpnu).at(mu).adjoint()*U_vec.at(n).at(nu).adjoint();
	// 		C_munu -= U_vec.at(n).at(mu)*U_vec.at(xmnupmu).at(nu).adjoint()*U_vec.at(xmnu).at(mu).adjoint()*U_vec.at(xmnu).at(nu);
	// 		C_munu += U_vec.at(xmmu).at(mu).adjoint()*U_vec.at(xmnummu).at(nu).adjoint()*U_vec.at(xmnummu).at(mu)*U_vec.at(xmnu).at(nu);
	// 		C_munu -= U_vec.at(xmmu).at(mu).adjoint()*U_vec.at(xmmu).at(nu)*U_vec.at(xpnummu).at(mu)*U_vec.at(n).at(nu).adjoint();

			C_munu  = (*(U_vec+n)).at(mu)*(*(U_vec+xpmu)).at(nu)*(*(U_vec+xpnu)).at(mu).adjoint()*(*(U_vec+n)).at(nu).adjoint();
			C_munu -= (*(U_vec+n)).at(mu)*(*(U_vec+xmnupmu)).at(nu).adjoint()*(*(U_vec+xmnu)).at(mu).adjoint()*(*(U_vec+xmnu)).at(nu);
			C_munu += (*(U_vec+xmmu)).at(mu).adjoint()*(*(U_vec+xmnummu)).at(nu).adjoint()*(*(U_vec+xmnummu)).at(mu)*(*(U_vec+xmnu)).at(nu);
			C_munu -= (*(U_vec+xmmu)).at(mu).adjoint()*(*(U_vec+xmmu)).at(nu)*(*(U_vec+xpnummu)).at(mu)*(*(U_vec+n)).at(nu).adjoint();
	//
			{
				Matrix2cd tmp = C_munu.adjoint();
				C_munu -= tmp;
			}
	//
	// 		C_rhosig  = U_vec.at(n).at(rho)*U_vec.at(xprho).at(sig)*U_vec.at(xpsig).at(rho).adjoint()*U_vec.at(n).at(sig).adjoint();
	// 		C_rhosig -= U_vec.at(n).at(rho)*U_vec.at(xmsigprho).at(sig).adjoint()*U_vec.at(xmsig).at(rho).adjoint()*U_vec.at(xmsig).at(sig);
	// 		C_rhosig += U_vec.at(xmrho).at(rho).adjoint()*U_vec.at(xmsigmrho).at(sig).adjoint()*U_vec.at(xmsigmrho).at(rho)*U_vec.at(xmsig).at(sig);
	// 		C_rhosig -= U_vec.at(xmrho).at(rho).adjoint()*U_vec.at(xmrho).at(sig)*U_vec.at(xpsigmrho).at(rho)*U_vec.at(n).at(sig).adjoint();

			C_rhosig  = (*(U_vec+n)).at(rho)*(*(U_vec+xprho)).at(sig)*(*(U_vec+xpsig)).at(rho).adjoint()*(*(U_vec+n)).at(sig).adjoint();
			C_rhosig -= (*(U_vec+n)).at(rho)*(*(U_vec+xmsigprho)).at(sig).adjoint()*(*(U_vec+xmsig)).at(rho).adjoint()*(*(U_vec+xmsig)).at(sig);
			C_rhosig += (*(U_vec+xmrho)).at(rho).adjoint()*(*(U_vec+xmsigmrho)).at(sig).adjoint()*(*(U_vec+xmsigmrho)).at(rho)*(*(U_vec+xmsig)).at(sig);
			C_rhosig -= (*(U_vec+xmrho)).at(rho).adjoint()*(*(U_vec+xmrho)).at(sig)*(*(U_vec+xpsigmrho)).at(rho)*(*(U_vec+n)).at(sig).adjoint();

			{
				Matrix2cd tmp = C_rhosig.adjoint();
				C_rhosig -= tmp;
			}
	//
	// 		// cout << (C_munu*C_rhosig).trace() << endl;
	//
			qq += params::signum_red.at(k)*real((C_munu*C_rhosig).trace());
		}
	}

	return -qq/256.0/params::pi/params::pi;

}

double meas_topo_plaq(vector<int>* neib, vector<Matrix2cd,aligned_allocator<Matrix2cd> >* U_vec)
{
	double qq = 0.0;
	Matrix2cd C_munu = Matrix2cd::Zero(), C_rhosig = Matrix2cd::Zero();
	int mu, nu, rho, sig;
	int ipmu, ipnu, iprho, ipsig;

	const vector<vector<int>> perms_red  = {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2}};
	const vector<int> signum_red = {1, -1, 1};

	for(int t = 0; t < params::Nt; t++)
	{
		for(int z = 0; z < params::Ns; z++)
		{
			for(int y = 0; y < params::Ns; y++)
			{
				for(int x = 0; x < params::Ns; x++)
				{
					int i = x + y*params::Ns + z*params::Ns*params::Ns + t*params::Ns*params::Ns*params::Ns;

					for(int k = 0; k < 3; k++)
					{
						mu = perms_red.at(k).at(0);
						nu = perms_red.at(k).at(1);
						rho = perms_red.at(k).at(2);
						sig = perms_red.at(k).at(3);

						ipmu = (*(neib+i)).at(mu);
						ipnu = (*(neib+i)).at(nu);
						iprho = (*(neib+i)).at(rho);
						ipsig = (*(neib+i)).at(sig);

						if(ipmu < 0 || ipnu < 0 || iprho < 0 || ipsig < 0) continue;

						C_munu = (*(U_vec+i)).at(mu)*(*(U_vec+ipmu)).at(nu)*(*(U_vec+ipnu)).at(mu).adjoint()*(*(U_vec+i)).at(nu).adjoint();
						C_rhosig = (*(U_vec+i)).at(rho)*(*(U_vec+iprho)).at(sig)*(*(U_vec+ipsig)).at(rho).adjoint()*(*(U_vec+i)).at(sig).adjoint();

						qq += signum_red.at(k)*real(((C_munu - C_munu.adjoint())*(C_rhosig - C_rhosig.adjoint())).trace());
					}
				}
			}
		}
	}
	return qq/32.0/params::pi/params::pi;
}

double get_local_topo_contribution_plaq(const int& site, const int& dir, vector<int>* neib, vector<Matrix2cd,aligned_allocator<Matrix2cd> >* U_vec)
{
	double qq = 0.0;
	Matrix2cd C_munu = Matrix2cd::Zero(), C_rhosig = Matrix2cd::Zero();
	int mu, nu, rho, sig;
	int ipmu, ipnu, iprho, ipsig;
	vector<int> dir_perp;

	const vector<vector<int>> perms_red  = {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2}};
	const vector<int> signum_red = {1, -1, 1};

	//get all direction (positive) perpendicular to dir
	switch(dir)
	{
		case 0:
			dir_perp = {1,2,3};
			break;
		case 1:
			dir_perp = {0,2,3};
			break;
		case 2:
			dir_perp = {0,1,3};
			break;
		case 3:
			dir_perp = {0,1,2};
			break;
		default:
			abort();
	}

	vector<int> site_vec = {site};
	for(auto it_dir: dir_perp)
	{
		site_vec.push_back((*(neib+site)).at(it_dir + params::dd));
	}

	for(auto i: site_vec)
	{
		if(i < 0) continue;

		for(int k = 0; k < 3; k++)
		{
			mu = perms_red.at(k).at(0);
			nu = perms_red.at(k).at(1);
			rho = perms_red.at(k).at(2);
			sig = perms_red.at(k).at(3);

			ipmu = (*(neib+i)).at(mu);
			ipnu = (*(neib+i)).at(nu);
			iprho = (*(neib+i)).at(rho);
			ipsig = (*(neib+i)).at(sig);

			if(ipmu < 0 || ipnu < 0 || iprho < 0 || ipsig < 0) continue;

			C_munu = (*(U_vec+i)).at(mu)*(*(U_vec+ipmu)).at(nu)*(*(U_vec+ipnu)).at(mu).adjoint()*(*(U_vec+i)).at(nu).adjoint();
			C_rhosig = (*(U_vec+i)).at(rho)*(*(U_vec+iprho)).at(sig)*(*(U_vec+ipsig)).at(rho).adjoint()*(*(U_vec+i)).at(sig).adjoint();

			qq += signum_red.at(k)*real(((C_munu - C_munu.adjoint())*(C_rhosig - C_rhosig.adjoint())).trace());
		}
	}

	return qq/32.0/params::pi/params::pi;
}
