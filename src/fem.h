#pragma once
#include <vector>
#include <map>
#include <set>
#include <string>
#include "mkl_pardiso.h"
#include "mkl_types.h"
#include <gmm/gmm.h>

#include "femoptions.h"
#include "node.h"
#include "element.h"
#include "elementset.h"
#include "nodeset.h"
#include "face.h"

class fem
{
public:
	typedef std::map<size_t, point> loadmap_t;

	enum im_direction_t
	{
		NORMAL,
		ANTI_NORMAL
	};

public:
	std::vector<node> nodes;
	size_t nodecnt;
	
	int udof;	//number of unconstrainted DOF

	std::vector<element> elements;
	size_t elementcnt;

	std::map<std::string, nodeset> nodesets;
	std::map<std::string, elementset> elementsets;

	//stiffness matrix
	//need to use row matrix type to aid in CSR conversion
	//This is the matrix in which the stiffness matrix is built
	gmm::row_matrix<gmm::wsvector<double> > K_gmm;
	
	//This a the compressed storage [column] format used by the gmm solvers
	gmm::csc_matrix<double>	K_gmm_csc;

	//node mass matrix (lumped form, so a diagonal vector)
	std::vector<double> nodes_mass;

	//
	std::map<face, size_t> externalelementfaces;

	//the gap set
	nodeset* gap_nodeset;

	//faces that makeup the gap set, used for searching
	std::set<const face*> gap_face_set;

	//faces that makeup the gap set, used to sequential loading
	std::vector<const face*> gap_face_vector;

	//gap node normals
	std::vector<point> gap_node_norms;

	//used for IR
	struct
	{
		double body_mass;
		point cog;
		double scale;
	} IR;

	//which way is defined as 'positive' for the im's
	im_direction_t im_direction;

	//MKL PARDISO solver data
	struct mkl_t
	{
		void* pt[64];
		_INTEGER_t maxfct;
		_INTEGER_t mnum;
		_INTEGER_t mtype;
		_INTEGER_t phase;
	
		_INTEGER_t iparm [64];
		_INTEGER_t msglvl;
		_INTEGER_t error;

		_INTEGER_t n;
		
		double * K_csr;
		_INTEGER_t * ia;
		_INTEGER_t * ja;

		_INTEGER_t * perm;	//not used

		mkl_t()
		{
			K_csr = NULL;
			ia = NULL;
			ja = NULL;
		}

		~mkl_t()
		{
			if(K_csr != NULL)
			{
				delete [] K_csr;
				K_csr = NULL;
			}
			if(ia != NULL)
			{
				delete [] ia;
				ia = NULL;
			}
			if(ja != NULL)
			{
				delete [] ja;
				ja = NULL;
			}
		}
		
	} mkl;

	//helper methods
	std::vector<const node*> nid2node(const std::vector<size_t> &nids);
	static std::vector<double> inverse_3x3(const std::vector<double> A);

public:
	//fem options
	femoptions options;
	
	//Member methods - IO
	void readoptions(std::string file);
	void loadinp();
	void writeVTK(std::string filename, const std::vector<double>& b, const std::vector<double>& x);
	void writeFacesetVTK(const std::string filename, const std::vector<const face*>& fs);
	void getim(const std::vector<double>::iterator x, std::vector<double>& im);

	//Member methods - loading
	void setupgap(std::string gap_set_name, im_direction_t direction);
	void nodeset2faceset(const nodeset& ns, std::vector<const face*>& fs);
	void load_gap_face(size_t faceid, double pressure, loadmap_t& loadmap);
	void load_set(std::string setname, double pressure, loadmap_t& loadmap);

	void loadmap2b(const loadmap_t& loadmap, std::vector<double>::iterator b);
	void loadmap2b(const loadmap_t& loadmap, std::vector<double> &b);
	
	//Member methods - helping
	void number_DOF();
	void stiffness_matrix(const bool only_triu = false);
	size_t get_gap_facecnt();
	size_t get_gap_nodecnt();
	int get_udof();

	//Solving - gmm
	void gmm_setup();
	bool gmm_solve(const loadmap_t& loadmap, std::vector<double>& b, std::vector<double>& x);

	//Solving - MKL
	void K_gmm2K_csr(bool only_triu = true);
	void mkl_phase0();
	void mkl_phase11();
	void mkl_phase22();
	void mkl_phase33(_INTEGER_t nrhs, double* b, double * x);
	void mkl_gap_ims(int startid, _INTEGER_t nrhs, std::vector<double>& b, std::vector<double> &x);
	void mkl_solve(const loadmap_t& loadmap, std::vector<double>& b, std::vector<double>& x);

};