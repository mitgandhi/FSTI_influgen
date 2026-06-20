#include "influgen.h"
#include "string_functions.h"
#include <iostream>
#include <stdlib.h>
#include "../lib/x64/fsti_netlic_client.h"

void help()
{
	cout << endl;
	cout << "Usage: " << endl;
	cout << "influgen_2.exe -INTERFACE -WORK [arg] [-WORK [arg]] ..." << endl;
	cout << endl;
	cout << "-INTERFACE must be only one of the following:" << endl;
	cout << "  -piston" << endl;
	cout << "  -bushing" << endl;
	cout << "  -block" << endl;
	cout << "  -valveplate" << endl;
	cout << "  -slipper" << endl;
	cout << "  -swashplate" << endl;
	cout << endl;
	cout << "-WORK may be any of the following:" << endl;
	cout << "  -imset SETNAME" << endl;
	cout << "  -vtkset SETNAME" << endl;
	cout << "  -vtkgap GAP_FACE_ID" << endl;
	cout << "  -imgap" << endl;
	cout << "  -vizvtksets (This will only create undeformed surface vtk's of the defined sets)" << endl;
	exit(2);
};

struct work
{
	enum types
	{
		imgap,
		imset,
		vtkgap,
		vtkset,
		makefiles,
		vizvtksets,
		none
	};
	types type;
	std::string arg;
};
void process_args(int argc, char *argv[], influGen::part_type& part, std::vector<work>& works, influGen::solver_type& solver)
{
	//set to none values
	part = influGen::NONE;
	works.clear();
	//always make the files
	{
		work w;
		w.type = work::makefiles;
		works.push_back(w);
	}

	//default the solver to PARDISO
	solver = influGen::PARDISO;

	for(int i=1; i<argc; ++i)
	{
		//first check the interface ty pes
		if(strcmp("-piston", argv[i]) == 0)
		{
			if(part != influGen::NONE)
			{
				//can't set two parts
				cout << "ERROR: Only set one interface!" << endl;
				help();
			}
			part = influGen::PISTON;
		} else if(strcmp("-bushing", argv[i]) == 0)
		{
			if(part != influGen::NONE)
			{
				//can't set two parts
				cout << "ERROR: Only set one interface!" << endl;
				help();
			}
			part = influGen::BUSHING;
		} else if(strcmp("-block", argv[i]) == 0)
		{
			if(part != influGen::NONE)
			{
				//can't set two parts
				cout << "ERROR: Only set one interface!" << endl;
				help();
			}
			part = influGen::BLOCK;
		} else if(strcmp("-valveplate", argv[i]) == 0)
		{
			if(part != influGen::NONE)
			{
				//can't set two parts
				cout << "ERROR: Only set one interface!" << endl;
				help();
			}
			part = influGen::VALVEPLATE;
		} else if(strcmp("-slipper", argv[i]) == 0)
		{
			if(part != influGen::NONE)
			{
				//can't set two parts
				cout << "ERROR: Only set one interface!" << endl;
				help();
			}
			part = influGen::SLIPPER;
		} else if(strcmp("-swashplate", argv[i]) == 0)
		{
			if(part != influGen::NONE)
			{
				//can't set two parts
				cout << "ERROR: Only set one interface!" << endl;
				help();
			}
			part = influGen::SWASHPLATE;
		}
		//now check the operations
		else if(strcmp("-imset", argv[i]) == 0)
		{
			work w;
			w.type = work::imset;
			if(i+1 < argc)
			{
				++i;
				w.arg = argv[i];
				works.push_back(w);
			} else {
				cout << "ERROR: -imset requires an argument!" << endl;
				help();
			}
		}
		else if(strcmp("-vtkset", argv[i]) == 0)
		{
			work w;
			w.type = work::vtkset;
			if(i+1 < argc)
			{
				++i;
				w.arg = argv[i];
				works.push_back(w);
			} else {
				cout << "ERROR: -vtkset requires an argument!" << endl;
				help();
			}
		}
		else if(strcmp("-vtkgap", argv[i]) == 0)
		{
			work w;
			w.type = work::vtkgap;
			if(i+1 < argc)
			{
				++i;
				w.arg = argv[i];
				works.push_back(w);
			} else {
				cout << "ERROR: -vtkgap requires an argument!" << endl;
				help();
			}
		}
		else if(strcmp("-imgap", argv[i]) == 0)
		{
			work w;
			w.type = work::imgap;
			works.push_back(w);

			//set the solver type to pardiso since we're generating the whole im gap
			solver = influGen::PARDISO;
		}
		else if(strcmp("-vizvtksets", argv[i]) == 0)
		{
			work w;
			w.type = work::vizvtksets;
			works.push_back(w);
		}
		else {
			if(string_functions::trim(argv[i]).length() > 0)
			{
				//argument not reconized
				cout << "ERROR: argument \"" << argv[i] << "\" not recognized!" << endl;
				help();
			}
		}
	}
	if(part == influGen::NONE)
	{
		cout << "ERROR: At least one interface must be specified!" << endl;
		help();
	}
	if(works.size() == 0)
	{
		cout << "ERROR: No work to do was defined!" << endl;
		help();
	}
}

int main(int argc, char *argv[])
{
	//the part type
	influGen::part_type part;
	//the solver type based on the operation to be preformed
	influGen::solver_type solver;
	//the work to be done	
	std::vector<work> works;
	//process the arguments
	process_args(argc, argv, part, works, solver);

	//check if influgen_mkl is authorized
	//convert the part enum to a string to pass to netlic just for logging on the netlic server
	std::string interface_name;
	switch(part)
	{
	case influGen::PISTON:
		interface_name = "PISTON";
		break;
	case influGen::BUSHING:
		interface_name = "BUSHING";
		break;
	case influGen::BLOCK:
		interface_name = "BLOCK";
		break;
	case influGen::VALVEPLATE:
		interface_name = "VALVEPLATE";
		break;
	case influGen::SLIPPER:
		interface_name = "SLIPPER";
		break;
	case influGen::SWASHPLATE:
		interface_name = "SWASHPLATE";
		break;
	case influGen::NONE:
		interface_name = "NONE";
		break;
	};
	void * license = fsti_netlic::checkout_license(fsti_netlic::NODE_THEN_NETWORK, interface_name);

	//create the two directories even if they really aren't needed
	system("IF NOT EXIST .\\vtk ( mkdir .\\vtk )");
	system("IF NOT EXIST .\\tmp ( mkdir .\\tmp )");

	//create the influGen instance
	influGen ig(part, solver);
	ig.load("./input/input.txt");

	//set the pardiso options even if just gmm will be used
	_putenv("MKL_PARDISO_OOC_PATH=.\\tmp\\influgen_tmp");
	if(ig.model.options.max_ram > 0)
	{
		_putenv(("MKL_PARDISO_OOC_MAX_CORE_SIZE=" + string_functions::n2s(ig.model.options.max_ram)).c_str());
	} else {
		_putenv("MKL_PARDISO_OOC_MAX_CORE_SIZE=2000");
	}
	_putenv("MKL_PARDISO_OOC_MAX_SWAP_SIZE=0");
	_putenv("MKL_PARDISO_OOC_KEEP_FILE=1");

	//generate the stiffness matrix
	ig.stiffness();

	//setup the remainder of the influgen process
	ig.setup("gap");

	//now do the work
	for(int w=0; w<works.size(); ++w)
	{
		switch(works[w].type)
		{
		case work::imset:
			if(solver == influGen::PARDISO)
			{
				ig.mkl_set_im(works[w].arg, false);
			} else {
				ig.gmm_set_im(works[w].arg, false);
			}
			break;
		case work::vtkset:
			if(solver == influGen::PARDISO)
			{
				ig.mkl_set_im(works[w].arg, true);
			} else {
				ig.gmm_set_im(works[w].arg, true);
			}
			break;
		case work::vtkgap:
			if(solver == influGen::PARDISO)
			{
				ig.mkl_gap_im(string_functions::s2i(works[w].arg), true);
			} else {
				ig.gmm_gap_im(string_functions::s2i(works[w].arg), true);
			}
			break;
		case work::imgap:
			//solver must always be PARDISO in this case
			ig.mkl_gap_ims();
			break;
		case work::makefiles:
			switch(part)
			{
			case influGen::PISTON:
			case influGen::BUSHING:
				ig.make_PistonBush_mesh_files();
				break;
			case influGen::BLOCK:
			case influGen::VALVEPLATE:
				ig.make_BlockVP_mesh_files();
				break;
			case influGen::SLIPPER:
			case influGen::SWASHPLATE:
				ig.make_SlipperSwash_mesh_files();
				break;
			}
			break;
		case work::vizvtksets:
			ig.make_vizvtksets();
			break;
		}
	}
	
	//clean
	ig.cleanup();

	//release the license
	fsti_netlic::release_license(license);

	return 0;
}
