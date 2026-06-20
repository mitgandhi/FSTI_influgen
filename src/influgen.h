#pragma once
#include <string>
#include <map>

#include "fem.h"


class influGen
{
public:
	enum part_type
	{
		PISTON,
		BUSHING,
		BLOCK,
		VALVEPLATE,
		SLIPPER,
		SWASHPLATE,
		NONE
	};

	enum solver_type
	{
		GMM,
		PARDISO
	};

	part_type part;
	const solver_type solver;
	fem model;
	
private:
	//these two functions are used to generate the correct file path for each interface
	std::string filename_gap(size_t id, bool is_vtk = false);
	std::string filename_gapim();
	std::string filename_set(std::string setname, bool is_vtk = false);

public:
	influGen(part_type Part, solver_type Solver);	
	void load(const std::string optionfile);
	void stiffness();
	void setup(const std::string gap_set_name);
	
	void make_SlipperSwash_mesh_files();
	void make_PistonBush_mesh_files();
	void make_BlockVP_mesh_files();
	
	void make_vizvtksets();

	void cleanup();
	void gmm_gap_im(size_t id, bool writevtk = false);
	void gmm_set_im(std::string setname, bool writevtk = false);
	void mkl_gap_ims();
	void mkl_set_im(std::string setname, bool writevtk);
	void mkl_gap_im(size_t id, bool writevtk = false);

};