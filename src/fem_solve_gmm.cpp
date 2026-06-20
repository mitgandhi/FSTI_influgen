#include "fem.h"

void fem::gmm_setup()
{
	gmm::copy(K_gmm,K_gmm_csc);
	K_gmm.clear_mat();
}
bool fem::gmm_solve(const loadmap_t& loadmap, std::vector<double>& b, std::vector<double>& x)
{
	//create the load vector
	loadmap2b(loadmap, b);

	//initialize the solution vector
	x.clear();
	x.resize(udof, 0);
	
	//Set solver params
	gmm::iteration iter(1e-12);
	iter.set_maxiter(2000);
	iter.set_noisy(0);
	
	gmm::ildlt_precond< gmm::csc_matrix<double> > PR(K_gmm_csc);
	//gmm::cg(K_gmm_csc, x, b, PR, iter);
	gmm::bicgstab(K_gmm_csc, x, b, PR, iter);

	if(iter.get_iteration() == iter.get_maxiter())
	{
		cout << "WARNING: Linear solver failed to converge in " << iter.get_maxiter() << " iterations!" << endl;
	} else {
		cout << "\tSolved linear system in " << iter.get_iteration() << " iterations." << endl;
	}

	if(options.inrel)
	{
		cout << "\tInertia Relief Trans. Acceleration = [ " << x[x.size()-3]*IR.scale << "\t"
															<< x[x.size()-2]*IR.scale << "\t"
															<< x[x.size()-1]*IR.scale << " ]" << endl;
		cout << "\tInertia Relief Rot. Acceleration = [ " << x[x.size()-6]*IR.scale << "\t"
															<< x[x.size()-5]*IR.scale << "\t"
															<< x[x.size()-4]*IR.scale << " ]" << endl;
	}

	if(iter.get_iteration() == iter.get_maxiter())
	{
		return false;
	}
	return true;
}
