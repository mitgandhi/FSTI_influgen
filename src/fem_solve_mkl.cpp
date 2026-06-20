#include "fem.h"
void fem::K_gmm2K_csr(bool only_triu)
{
	//this will load the K_gmm into the K_csr
	//note that K_gmm will be cleared in the process

	//first determint the nnz
	int nnz = 0;
	for(int i=0; i<udof; i++)
	{
		for(auto it = K_gmm[i].begin(); it != K_gmm[i].end(); ++it)
		{
			//only store the upper half if !only_triu
			if(it->first >= i || !only_triu)
			{
				++nnz;
			}
		}
	}
	mkl.K_csr = new double [nnz];
	mkl.ia = new _INTEGER_t [udof+1];
	mkl.ja = new _INTEGER_t [nnz];

	//fill the K_csr matrix
	int ctr = 0;
	for(int i=0; i<udof; i++)
	{
		mkl.ia[i] = ctr;
		for(auto it = K_gmm[i].begin(); it != K_gmm[i].end(); ++it)
		{
			//only store the upper half if !only_triu
			if(it->first >= i || !only_triu)
			{
				mkl.ja[ctr] = (_INTEGER_t) it->first;
				mkl.K_csr[ctr] = it->second;
				++ctr;
			}
		}

		//clear the current K_gmm row to minimize memory footprint
		K_gmm[i].clear();
	}
	mkl.ia[udof] = nnz;

};
void fem::mkl_phase0()
{
	//not needed during this phase
	_INTEGER_t nrhs = 0;
	double * b = NULL;
	double * x = NULL;

	//actually phase -1, but can't make a method that name
	mkl.phase = -1;
	pardiso(mkl.pt, &mkl.maxfct, &mkl.mnum, &mkl.mtype, &mkl.phase, &mkl.n, mkl.K_csr, mkl.ia, mkl.ja, mkl.perm, &nrhs, mkl.iparm, &mkl.msglvl, b, x, &mkl.error);

	if(mkl.error != 0)
	{
		std::cout << "ERROR: PHASE(" << mkl.phase << "). ERROR NO = " << mkl.error << "\n";
		exit(1);
	}
};
void fem::mkl_phase11()
{
	//this pass should ONLY be called once AFTER the K matrix has been formed and constrained.

	//ONLY used once to zero the internal MKL pointers
	for (int i = 0; i < 64; i++) {
		mkl.pt[i] = 0;
	}

	//only one matrix
	mkl.maxfct = 1;
	mkl.mnum = 1;

	//mtype -2 is real and symmetric indefinite
	//mtype 11 is real and unsymmetric - seems to be necessary for IR lamba constraint analysis w.o. correct b vector
	mkl.mtype = -2;

	//solver options
	for (int i = 0; i < 64; i++) {
		mkl.iparm[i] = 0;
	}

	//debugging message level (0: off, 1: on)
	mkl.msglvl = 0;

	//set the matrix size (n x n)
	mkl.n = udof;

	//SETUP SOLVER OPTIONS
	mkl.iparm[0] = 1; // No solver default 
	mkl.iparm[1] = 2; //Fill-in reordering from METIS
	mkl.iparm[3] = 0; // CGS 
	mkl.iparm[4] = 0; // No user fill-in reducing permutation 
	mkl.iparm[5] = 0; // Write solution into x 
	mkl.iparm[6] = 0; // Not in use 
	mkl.iparm[7] = 0; // Max numbers of iterative refinement steps 
	mkl.iparm[8] = 0; // Not in use 
	mkl.iparm[9] = 13; // Perturb the pivot elements with 1E-13 
	mkl.iparm[10] = 1; // Use nonsymmetric permutation and scaling MPS 
	mkl.iparm[11] = 0; // Not in use 
	mkl.iparm[12] = 1; // improved accuracy using (non-)symmetric weighted matchings.
	mkl.iparm[13] = 0; 
	mkl.iparm[14] = 0; // Not in use 
	mkl.iparm[15] = 0; // Not in use 
	mkl.iparm[16] = 0; // Not in use 
	mkl.iparm[17] = -1; // Output: Number of nonzeros in the factor LU 
	mkl.iparm[18] = -1; // Output: Mflops for LU factorization 
	mkl.iparm[20] = 1; // pivoting for symmetric indefinite matrices. 
	mkl.iparm[26] = 0; // check the data structure 
	mkl.iparm[27] = 0; // Double (0) or single (1) precision 
	mkl.iparm[30] = 0; //  partial solution for sparse right-hand sides and sparse solution.
	mkl.iparm[34] = 1; // zero indexing
	mkl.iparm[59] = 1; // OOC - only go out of core if over RAM

	//build the K_csr matrix
	K_gmm2K_csr();

	//not needed during this phase
	_INTEGER_t nrhs = 0;
	double * b = NULL;
	double * x = NULL;

	cout << "    Number of matrix non-zeros: " << mkl.ia[udof] << endl;

	mkl.phase = 11;
	pardiso(mkl.pt, &mkl.maxfct, &mkl.mnum, &mkl.mtype, &mkl.phase, &mkl.n, mkl.K_csr, mkl.ia, mkl.ja, mkl.perm, &nrhs, mkl.iparm, &mkl.msglvl, b, x, &mkl.error);
	if(mkl.error != 0)
	{
		std::cout << "ERROR: PHASE(" << mkl.phase << "). ERROR NO = " << mkl.error << "\n";
		exit(1);
	}

	cout << "    Number of matrix factors: " << mkl.iparm[17] << endl;
	cout << "    Peak RAM required for matrix analysis: " << mkl.iparm[14]/1024 << " MB" << endl;
	cout << "    Permanent RAM required for matrix analysis: " << mkl.iparm[15]/1024 << " MB" << endl;
	cout << "    Gflops required for matrix factorization: " << mkl.iparm[18]/1024 << endl;
	
	//need to use a double to prevent overflow in the case of large matrices
	double ram_Req = mkl.iparm[17];
	ram_Req *= 8;
	ram_Req /= 1024*1024;
	ram_Req += mkl.iparm[15]/1024;
	cout << "    Approximate RAM required for matrix factorization: " << floor(ram_Req) << " MB" << endl;


};
void fem::mkl_phase22()
{
	//not needed during this phase
	_INTEGER_t nrhs = 0;
	double * b = NULL;
	double * x = NULL;

	mkl.phase = 22;
	pardiso(mkl.pt, &mkl.maxfct, &mkl.mnum, &mkl.mtype, &mkl.phase, &mkl.n, mkl.K_csr, mkl.ia, mkl.ja, mkl.perm, &nrhs, mkl.iparm, &mkl.msglvl, b, x, &mkl.error);

	if(mkl.error != 0)
	{
		std::cout << "ERROR: PHASE(" << mkl.phase << "). ERROR NO = " << mkl.error << "\n";
		exit(1);
	}

	cout << "    Peak RAM required for matrix factorization: " << mkl.iparm[16]/1024+mkl.iparm[15]/1024 << " MB" << endl;
};
void fem::mkl_phase33(_INTEGER_t nrhs, double* b, double * x)
{
	mkl.phase = 33;
	pardiso(mkl.pt, &mkl.maxfct, &mkl.mnum, &mkl.mtype, &mkl.phase, &mkl.n, mkl.K_csr, mkl.ia, mkl.ja, mkl.perm, &nrhs, mkl.iparm, &mkl.msglvl, b, x, &mkl.error);

	if(mkl.error != 0)
	{
		std::cout << "ERROR: PHASE(" << mkl.phase << "). ERROR NO = " << mkl.error << "\n";
		exit(1);
	}
};
void fem::mkl_gap_ims(int startid, _INTEGER_t nrhs, std::vector<double>& b, std::vector<double> &x)
{
	//setup the vectors
	b.clear();
	b.resize(udof*nrhs, 0);
	x.clear();
	x.resize(udof*nrhs, 0);

	//get the load for each face
	for(int i=0; i<nrhs; ++i)
	{
		fem::loadmap_t loadmap;
		load_gap_face(startid+i, 100e5, loadmap);
		loadmap2b(loadmap, b.begin()+(udof*i));
	}

	mkl_phase33(nrhs, &b[0], &x[0]);
}
void fem::mkl_solve(const loadmap_t& loadmap, std::vector<double>& b, std::vector<double>& x)
{
	//setup the vectors
	b.clear();
	b.resize(udof, 0);
	x.clear();
	x.resize(udof, 0);
	
	loadmap2b(loadmap, b.begin());
	mkl_phase33(1, &b[0], &x[0]);
}
