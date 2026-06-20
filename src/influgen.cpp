#include "influgen.h"
#include "string_functions.h"
#include <iostream>
#include <iomanip>

using namespace std;
using namespace string_functions;

influGen::influGen(part_type Part, solver_type Solver) : solver(Solver)
{
	part = Part;
};

void influGen::load(const std::string optionfile)
{
	//load the options file
	cout << "Loading the option file \"" << optionfile << "\" ... ";
	model.readoptions(optionfile);
	cout << "done." << endl;
	
	cout << "Loading the mesh file \"" << model.options.meshfile << "\" ... " << endl;
	model.loadinp();
	cout << endl;
};
void influGen::setup(const std::string gap_set_name)
{
	cout << "Initializing the gap ..." << endl;
	
	fem::im_direction_t im_direction;
	switch(part)
	{
	case PISTON:
		im_direction = fem::NORMAL;
		break;
	case BUSHING:
		im_direction = fem::ANTI_NORMAL;
		break;
	case BLOCK:
		im_direction = fem::ANTI_NORMAL;
		break;
	case VALVEPLATE:
		im_direction = fem::NORMAL;
		break;
	case SLIPPER:
		im_direction = fem::ANTI_NORMAL;
		break;
	case SWASHPLATE:
		im_direction = fem::NORMAL;
		break;
	}

	model.setupgap(gap_set_name, im_direction);

	if(solver == PARDISO)
	{
		cout << "Initializing the linear solver ..." << endl;
		model.mkl_phase11();
		cout << "  Finished!" << endl;

		cout << "Factoring the stiffness matrix ..." << endl;
		model.mkl_phase22();
		cout << "  Finished!" << endl;
	} else {
		cout << "Initializing the linear solver ..." << endl;
		model.gmm_setup();
		cout << "  Finished!" << endl;
	}
};
void influGen::stiffness()
{
	//gmm requires the whole matrix, pardiso only half
	if(solver == PARDISO)
	{
		//create the stiffness matrix
		model.stiffness_matrix(true);
	} else {
		//create the stiffness matrix
		model.stiffness_matrix(false);
	}

}
void influGen::cleanup()
{
	if(solver == PARDISO)
	{
		model.mkl_phase0();
	}
}

void influGen::gmm_gap_im(size_t id, bool writevtk)
{
	if(solver != GMM)
	{
		cout << "ERROR: Wrong solver!" << endl;
		exit(1);
	}

	//the gap face cnt
	int gap_facecnt = (int) model.get_gap_facecnt();

	//the gap node cnt
	int gap_nodecnt = (int) model.get_gap_nodecnt();

	cout << endl;
	if(id < 0 || id >= gap_facecnt)
	{
		cout << "ERROR: The requested gap IM must be in the range of [0, " << gap_facecnt << "]" << endl;
		exit(1);
	}
	cout << "Solving for gap IM #: " << id << endl;
	fem::loadmap_t loadmap;
	model.load_gap_face(id, 100e5, loadmap);
	std::vector<double> b;	//the nodal load vector (used for vtk)
	std::vector<double> x;	//the solution
	model.gmm_solve(loadmap, b, x);
	if(writevtk)
	{
		model.writeVTK(filename_gap(id, writevtk), b, x);
	} else {			
		//this assumes this function is being called in gapset order
		fstream imfile(filename_gapim().c_str(), ios::out|ios::binary|ios::app);

		//get each IM from the global solution vector
		vector<double> im;
		model.getim(x.begin(), im);

		for(size_t n=0; n<gap_nodecnt; n++)
		{
			imfile.write((char*) &(im[n]), sizeof(double));
		}

		imfile.close();
	}
}
void influGen::gmm_set_im(std::string setname, bool writevtk)
{
	if(solver != GMM)
	{
		cout << "ERROR: Wrong solver!" << endl;
		exit(1);
	}

	cout << "Solving for set: " << setname << endl;
	fem::loadmap_t loadmap;
	model.load_set(setname, 100e5, loadmap);
	std::vector<double> b;	//the nodal load vector (used for vtk)
	std::vector<double> x;	//the solution
	model.gmm_solve(loadmap, b, x);
	if(writevtk)
	{
		model.writeVTK(filename_set(setname, writevtk), b, x);
	} else {		
		if(part == BLOCK || part == VALVEPLATE)
		{
			//write the block / valveplate as a txt
			fstream imfile(filename_set(setname, false).c_str(), ios::out);

			//get the IM from the global solution vector
			vector<double> im;
			model.getim(x.begin(), im);
								
			for(size_t n=0; n<im.size(); n++)
			{
				imfile << setprecision(16) << im[n] << endl;
			}		

			imfile.close();

		} else {
			//and everything else as binary

			fstream imfile(filename_set(setname, false).c_str(), ios::out|ios::binary);

			//get the IM from the global solution vector
			vector<double> im;
			model.getim(x.begin(), im);
								
			for(size_t n=0; n<im.size(); n++)
			{
				imfile.write((char*) &(im[n]), sizeof(double));
			}		

			imfile.close();
		}
	}
}
void influGen::mkl_gap_ims()
{
	if(solver != PARDISO)
	{
		cout << "ERROR: Wrong solver!" << endl;
		exit(1);
	}

	//the gap face cnt
	int gap_facecnt = (int) model.get_gap_facecnt();

	//the gap node cnt
	int gap_nodecnt = (int) model.get_gap_nodecnt();

	//how many rhs will we solve at once
	int nrhs = 100;	//should be tuneable, either external or based on DOF

	//how many times will we need to call the pardiso solve
	int nsolves = (gap_facecnt + nrhs - 1) / nrhs;

	fstream imfile(filename_gapim().c_str(), ios::out|ios::binary);
	for(int i=0; i<nsolves; ++i)
	{
		int startid = nrhs*i;
		int nextid = nrhs*(i+1);
		if(nextid > gap_facecnt)
		{
			nextid = gap_facecnt;
		}

		cout << "Solving for gap IM #: " << startid << " to " << nextid-1 << " ( " << i*100/nsolves << "% )" << endl;

		std::vector<double> b; //load vector
		std::vector<double> x; //solution vector
		model.mkl_gap_ims(startid, nextid-startid, b, x);

		//write the IM's
		{			
			for(int j=0; j<nextid-startid; ++j)
			{
				//get each IM from the global solution vector
				vector<double> im;
				model.getim(x.begin()+(model.get_udof()*j), im);

				for(size_t n=0; n<gap_nodecnt; n++)
				{
					imfile.write((char*) &(im[n]), sizeof(double));
				}
			}
		}
	}
	imfile.close();

}
void influGen::mkl_gap_im(size_t id, bool writevtk)
{
	if(solver != PARDISO)
	{
		cout << "ERROR: Wrong solver!" << endl;
		exit(1);
	}

	//the gap face cnt
	int gap_facecnt = (int) model.get_gap_facecnt();

	//the gap node cnt
	int gap_nodecnt = (int) model.get_gap_nodecnt();

	cout << endl;
	if(id < 0 || id >= gap_facecnt)
	{
		cout << "ERROR: The requested gap IM must be in the range of [0, " << gap_facecnt << "]" << endl;
		exit(1);
	}
	cout << "Solving for gap IM #: " << id << endl;
	fem::loadmap_t loadmap;
	model.load_gap_face(id, 100e5, loadmap);
	std::vector<double> b;	//the nodal load vector (used for vtk)
	std::vector<double> x;	//the solution
	model.mkl_solve(loadmap, b, x);
	if(writevtk)
	{
		model.writeVTK(filename_gap(id, writevtk), b, x);
	} else {			
		//this assumes this function is being called in gapset order
		fstream imfile(filename_gapim().c_str(), ios::out|ios::binary|ios::app);

		//get each IM from the global solution vector
		vector<double> im;
		model.getim(x.begin(), im);

		for(size_t n=0; n<gap_nodecnt; n++)
		{
			imfile.write((char*) &(im[n]), sizeof(double));
		}

		imfile.close();
	}
}
void influGen::mkl_set_im(std::string setname, bool writevtk)
{
	if(solver != PARDISO)
	{
		cout << "ERROR: Wrong solver!" << endl;
		exit(1);
	}

	cout << "Solving for set: " << setname << endl;
	fem::loadmap_t loadmap;
	model.load_set(setname, 100e5, loadmap);
	std::vector<double> b;	//the nodal load vector (used for vtk)
	std::vector<double> x;	//the solution
	model.mkl_solve(loadmap, b, x);
	if(writevtk)
	{
		model.writeVTK(filename_set(setname, writevtk), b, x);
	} else {			
		if(part == BLOCK || part == VALVEPLATE)
		{
			//write the block / valveplate as a txt
			fstream imfile(filename_set(setname, false).c_str(), ios::out);

			//get the IM from the global solution vector
			vector<double> im;
			model.getim(x.begin(), im);
								
			for(size_t n=0; n<im.size(); n++)
			{
				imfile << setprecision(16) << im[n] << endl;
			}		

			imfile.close();

		} else {
			//and everything else as binary

			fstream imfile(filename_set(setname, false).c_str(), ios::out|ios::binary);

			//get the IM from the global solution vector
			vector<double> im;
			model.getim(x.begin(), im);
								
			for(size_t n=0; n<im.size(); n++)
			{
				imfile.write((char*) &(im[n]), sizeof(double));
			}		

			imfile.close();
		}
	}
}
void influGen::make_SlipperSwash_mesh_files()
{
	cout << "Performing initial processing ... ";

	//write the gapface set to file
	ofstream file;
	file.precision(8);
	
	file.open("./faces.txt");
	
	//number the gap nodes
	std::map<size_t, int> nid2gap;
	{
		int renumber = 0;
		for(auto n=model.gap_nodeset->begin(); n!=model.gap_nodeset->end(); ++n)
		{
			nid2gap[*n] = renumber;
			++renumber;
		}
	}

	for(auto f=model.gap_face_vector.begin(); f!=model.gap_face_vector.end(); ++f)
	{
		for(size_t n=0; n<(*f)->nodes.size(); ++n)
		{
			file << nid2gap[(*f)->nodes[n]->nid] << "\t";
		}
		file << scientific;
		point centroid = (*f)->centroid();
		file << centroid.x() << "\t" << centroid.y() << endl;
	}

	file.close();

	file.open("./nodes.txt");

	size_t n_gap = 0;
	for(auto n=model.gap_nodeset->begin(); n!=model.gap_nodeset->end(); ++n)
	{
		file << nid2gap[*n] << "\t";
		file << scientific <<
			model.nodes[*n].x() << "\t" <<
			model.nodes[*n].y() << 
			endl;
		++n_gap;
	}

	file.close();

	cout << "done." << endl;

};
void influGen::make_PistonBush_mesh_files()
{
	cout << "Performing initial processing ... ";

	//write the gapface set to file
	ofstream file;
	file.precision(8);
	
	file.open("./xyz_faces.dat");
	
	file << model.gap_face_vector.size() << endl;
	for(auto f=model.gap_face_vector.begin(); f!=model.gap_face_vector.end(); ++f)
	{
		file << scientific;
		point centroid = (*f)->centroid();
		file << centroid.x() << "\t" << centroid.y() << "\t" << centroid.z() << endl;
	}

	file.close();

	file.open("./xyz_nodes.dat");

	file << model.gap_nodeset->size() << endl;
	for(auto n=model.gap_nodeset->begin(); n!=model.gap_nodeset->end(); ++n)
	{
		file << scientific <<
			model.nodes[*n].x() << "\t" <<
			model.nodes[*n].y() << "\t" <<
			model.nodes[*n].z() << 
			endl;
	}

	file.close();

	cout << "done." << endl;

};
void influGen::make_BlockVP_mesh_files()
{
	cout << "Performing initial processing ... ";

	//write the gapface set to file
	ofstream file;
	file.precision(8);
	
	file.open("./gap.txt");

	std::map<size_t, int> nid2gap;

	file << "NODES " << model.gap_nodeset->size() << endl;
	int gapnodecnt = 0;
	for(auto n=model.gap_nodeset->begin(); n!=model.gap_nodeset->end(); ++n)
	{
		nid2gap[model.nodes[*n].nid] = gapnodecnt;
		file << scientific <<
			model.nodes[*n].x() << "\t" <<
			model.nodes[*n].y() << "\t" <<
			model.nodes[*n].z() << 
			endl;
		++gapnodecnt;
	}

	file << "ELEMENTS " << model.gap_face_vector.size() << endl;
	for(auto f=model.gap_face_vector.begin(); f!=model.gap_face_vector.end(); ++f)
	{
		for(size_t n=0; n<(*f)->nodes.size(); ++n)
		{
			file << nid2gap[(*f)->nodes[n]->nid] << "\t";
		}
		file << endl;
	}

	file.close();

	cout << "done." << endl;

};

void influGen::make_vizvtksets()
{
	for(auto set=model.nodesets.begin(); set!=model.nodesets.end(); ++set)
	{
		std::vector<const face*> fs;
		model.nodeset2faceset(set->second, fs);
		model.writeFacesetVTK("./vtk/vizvtkset_"+set->first+".vtk", fs);
	}
}

std::string influGen::filename_gap(size_t id, bool is_vtk)
{
	string filename = "";


	switch(part)
	{
	case PISTON:
		if(is_vtk)
		{
			filename = "./vtk/im_piston_gap_" + n2s(id) +".vtk";
		} else {
			filename = "./im_piston_gap_" + n2s(id) +".bin";
		}
		break;
	case BUSHING:
		if(is_vtk)
		{
			filename = "./vtk/im_bushing_gap_" + n2s(id) +".vtk";
		} else {
			filename = "./im_bushing_gap_" + n2s(id) +".bin";
		}
		break;
	case BLOCK:
		if(is_vtk)
		{
			filename = "./vtk/IM.gap." + n2s(id) +".vtk";
		} else {
			filename = "./IM.gap." + n2s(id) +".bin";
		}
		break;
	case VALVEPLATE:
		if(is_vtk)
		{
			filename = "./vtk/IM.gap." + n2s(id) +".vtk";
		} else {
			filename = "./IM.gap." + n2s(id) +".bin";
		}
		break;
	case SLIPPER:
		if(is_vtk)
		{
			filename = "./vtk/im_slipper_gap_" + n2s(id) +".vtk";
		} else {
			filename = "./im_slipper_gap_" + n2s(id) +".bin";
		}
		break;
	case SWASHPLATE:
		if(is_vtk)
		{
			filename = "./vtk/im_swashplate_gap_" + n2s(id) +".vtk";
		} else {
			filename = "./im_swashplate_gap_" + n2s(id) +".bin";
		}
		break;
	}

	return filename;
}
std::string influGen::filename_set(std::string setname, bool is_vtk)
{
	string filename = "";

	switch(part)
	{
	case PISTON:
		if(is_vtk)
		{
			filename = "./vtk/" + setname + ".vtk";
		} else {
			filename = "./" + setname + ".bin";
		}
		break;
	case BUSHING:
		if(is_vtk)
		{
			filename = "./vtk/" + setname + ".vtk";
		} else {
			filename = "./" + setname + ".bin";
		}
		break;
	case BLOCK:
		if(is_vtk)
		{
			filename = "./vtk/IM." + setname + ".vtk";
		} else {
			filename = "./IM." + setname + ".txt";
		}
		break;
	case VALVEPLATE:
		if(is_vtk)
		{
			filename = "./vtk/IM." + setname + ".vtk";
		} else {
			filename = "./IM." + setname + ".txt";
		}
		break;
	case SLIPPER:
		if(is_vtk)
		{
			filename = "./vtk/im_slipper_" + setname + ".vtk";
		} else {
			filename = "./im_slipper_" + setname + ".bin";
		}
		break;
	case SWASHPLATE:
		if(is_vtk)
		{
			filename = "./vtk/im_swashplate_" + setname + ".vtk";
		} else {
			filename = "./im_swashplate_" + setname + ".bin";
		}
		break;
	}

	return filename;
}
std::string influGen::filename_gapim()
{
	string filename = "";

	switch(part)
	{
	case PISTON:
		filename = "./im_piston_gap.bin";
		break;
	case BUSHING:
		filename = "./im_bushing_gap.bin";
		break;
	case BLOCK:
		filename = "./IM.gap.bin";
		break;
	case VALVEPLATE:
		filename = "./IM.gap.bin";
		break;
	case SLIPPER:
		filename = "./im_slipper_gap.bin";
		break;
	case SWASHPLATE:
		filename = "./im_swashplate_gap.bin";
		break;
	}

	return filename;
}
