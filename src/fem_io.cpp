#include "fem.h"
#include "string_functions.h"

#include <map>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
using namespace string_functions;


void fem::readoptions(string file)
{
	//used in this type of option file reader
	struct section
	{
		string name;
		vector<istringstream*> fields;
	};


	ifstream in(file.c_str());
	if (!in.is_open()) 	
	{
		cout << "input::read(): Unable to open " << file << " input file!" << endl;
		exit(1);
	}

	string tmp;
	string line;
	double val = 0;
	vector<section*> sections;

	//read the file
	while(getline(in,line)) 
	{
		if(line.size() > 0)
		{
			istringstream iss (line,istringstream::in);
			iss >> tmp;
			size_t comment  = tmp.find("//"); // check if there is a comment
			// if is not a comment read
			if (comment == string::npos) 
			{
				// new section
				if(tmp.compare("section") == 0)
				{
					section* s = new section;
					sections.push_back(s);
					iss >> tmp;	
					if(tmp.compare("section") != 0)
					{
						s -> name = tmp;
						
						while(tmp.compare("endsection") != 0)
						{
							getline(in,line);
							if(line.size() > 0)
							{
								tmp.resize(0);
								istringstream iss(line,istringstream::in);
								iss >> tmp;
								if(tmp.size() > 0)
								{
									s -> fields.push_back(new istringstream (line,istringstream::in));
								}
							}
						}
					}
					else
					{
						cout << "Missing section name!" << endl;
						exit(1);
					}
				}
			}
		}
	}

	//process the sections
	string fieldname;
	string fieldvalue;

	bool read_general = false;
	bool read_FEM_solver = false;

	for(unsigned int i=0; i<sections.size(); i++)
	{
		// read section general
		if((sections[i] -> name).compare("general") == 0)
		{
			read_general = true;
			int validfields = 0;
			for(unsigned int j=0; j<sections[i] -> fields.size(); j++)
			{
				*(sections[i] -> fields[j]) >> fieldname;
				if(fieldname.compare("pump") == 0)
					*(sections[i] -> fields[j]) >> options.pump_name, validfields++;
				if(fieldname.compare("meshFile") == 0)
					*(sections[i] -> fields[j]) >> options.meshfile, validfields++;
				if(fieldname.compare("inrel") == 0)
					*(sections[i] -> fields[j]) >> options.inrel, validfields++;
				if(fieldname.compare("ram") == 0)
					*(sections[i] -> fields[j]) >> options.max_ram, validfields++;
				if(fieldname.compare("dimension") == 0)
				{
					string unit;
					*(sections[i] -> fields[j]) >> unit;
					if(unit.compare("m") == 0)
						options.scalefactor = 1.0,  validfields++;
					else if(unit.compare("mm") == 0)
						options.scalefactor = 1.0e-3,  validfields++;
					else
					{
						cout << "In file " << file 
								 << " wrong dimension specified: "
								 << unit << ". Please use either m or mm." << endl;
						exit(1);
					}
				}
			}
		}
		
		// read new material
		if((sections[i] -> name).compare("material") == 0)
		{
		
			material_type mat;

			int validfields = 0;
			for(unsigned int j=0; j<sections[i] -> fields.size(); j++)
			{
				*(sections[i] -> fields[j]) >> fieldname;
				if(fieldname.compare("name") == 0)
					*(sections[i] -> fields[j]) >> mat.name, validfields++;
				else if(fieldname.compare("E") == 0)
					*(sections[i] -> fields[j]) >> mat.E, validfields++;
				else if(fieldname.compare("nu") == 0)
					*(sections[i] -> fields[j]) >> mat.nu, validfields++;
				else if(fieldname.compare("rho") == 0)
					*(sections[i] -> fields[j]) >> mat.rho, validfields++;
				else if(fieldname.compare("sets") == 0)
				{
					validfields++;
					int nsets;
					*(sections[i] -> fields[j]) >> nsets;
					for(int h = 0; h < nsets; h++)
					{
						string thisset;
						*(sections[i] -> fields[j]) >> thisset;
						//force the set name to lower case
						thisset = str2lower(thisset);
						
						if(options.elementset_material.count(thisset) > 0)
						{
							cout << "Element set (" << thisset << ") has multiple materials assigned to it!" << endl;
							exit(1);
						}
						options.elementset_material[thisset] = options.materials.size();
					}
				}
				else if(fieldname.compare("//") != 0 && fieldname.compare("endsection") != 0)
				{
					cout << "Wrong keyword (" << fieldname << ") in material " 
						 << sections[i] -> name << " definition!" << endl;
					exit(1);
				}
			}
			if(validfields != 5)
			{
				cout << "Missing some field in material " << fieldname 
						 << " definition!" << endl;
				cout << validfields << endl;
				exit(1);
			}

			options.materials.push_back(mat);
		}

		// read new constraint
		if((sections[i] -> name).compare("constraint") == 0)
		{
			femoptions::elastic_constraint constraint;
			std::string setname = "";

			int validfields = 0;
			for(unsigned int j=0; j<sections[i] -> fields.size(); j++)
			{
				*(sections[i] -> fields[j]) >> fieldname;
				if(fieldname.compare("set") == 0)
				{	
					*(sections[i] -> fields[j]) >> setname, validfields++;
					//force the set name to lower case
					setname = str2lower(setname);
				}
				else if(fieldname.compare("x") == 0)
				{
					*(sections[i] -> fields[j]) >> constraint.x, validfields++;
				}
				else if(fieldname.compare("y") == 0)
				{
					*(sections[i] -> fields[j]) >> constraint.y, validfields++;
				}
				else if(fieldname.compare("z") == 0)
				{
					*(sections[i] -> fields[j]) >> constraint.z, validfields++;
				}
				else if(fieldname.compare("//") != 0 && fieldname.compare("endsection") != 0)
				{
					cout << "Wrong keyword (" << fieldname << ") in " 
						 << sections[i] -> name << " definition!" << endl;
					exit(1);
				}
			}
			
			if(setname.length() > 0)
			{
				//normal constraint
				if(options.elastic_constraints.count(setname) > 0)
				{
					cout << "ERROR: Constraint set ( " << setname << " ) is defined multiple times!" << endl;
					exit(1);
				}
				options.elastic_constraints[setname] = constraint;
			} 
			
			//either empty section or inertia relief so do nothing further

		}

	}	//done processing the sections

	//check for mandatory sections

	if(read_general == false)
	{
		cout << endl << "In input file " << file
				 << " the general section ismissing!" << endl;
		exit(1);
	}

};

void fem::loadinp()
{
	const string inpfile = options.meshfile;
	const double scale = options.scalefactor;

	nodecnt = 0;
	elementcnt = 0;

	ifstream inp(inpfile.c_str());
	if(!inp.is_open())
	{
		cout << "Error! Unable to open: " << inpfile << endl;
		exit(1);
	}

	//node id map used to reorder the nodes consecutively
	map<int, size_t> node_renumber;

	while(!inp.eof())
	{
		string l;
		getline(inp, l);;

		if(l.find("*NODE") == 0)
		{
			cout << "\tReading nodes... ";
			while(true)
			{
				streamoff pos = inp.tellg();
				getline(inp, l);
				if(l.find("*") != string::npos)
				{
					inp.seekg(pos);
					break;	//finished with the node section
				}
				vector<string> line = Tokenize(l, ",");
				nodes.push_back(node(scale*s2d(line[1]), scale*s2d(line[2]), scale*s2d(line[3])));
				nodes[nodes.size()-1].nid = nodes.size()-1;
				node_renumber[s2i(line[0])] = nodecnt;
				nodecnt++;
			}

			cout << "read " << nodes.size() << " nodes." << endl;

			//restart the main file loop
			continue;
		}

		/*
		if(l.find("*ELEMENT,TYPE=S3,ELSET=") == 0)
		{
			//get the face set name
			string name = l.substr(23);
			cout << "\tReading face set \"" << name << "\" ... ";

			while(true)
			{
				streamoff pos = inp.tellg();
				getline(inp, l);
				if(l.find("*") != string::npos)
				{
					inp.seekg(pos);
					break;	//finished with the node section
				}
				vector<string> line = Tokenize(l, ",");
				
				face Face;
				for(unsigned int i=1; i<line.size(); i++)
				{
					int nid = s2i(line[i])-1;
					if(nid < 0 || nid >= nodecnt)
					{
						cout << "Error! Problem with nid = " << nid << " of face = " << line[0] << " in faceset = " << name << endl;
						exit(1);
					}
					Face.push_back(&nodes[nid]);
				}

				facesets[name].push_back(Face);
			}
			cout << "read " << facesets[name].size() << " faces." << endl;
			
			//restart the main file loop
			continue;
		}
		*/

		if(l.find("*ELEMENT,TYPE=C3D4,ELSET=") == 0)
		{
			//get the element set name
			string name = l.substr(25);
			//force the setname to lower
			name = str2lower(name);
			cout << "\tReading element set \"" << name << "\" ... ";
			
			int set_elecnt = 0;

			while(true)
			{
				streamoff pos = inp.tellg();
				getline(inp, l);
				if(l.find("*") != string::npos)
				{
					inp.seekg(pos);
					break;	//finished with the node section
				}
				vector<string> line = Tokenize(l, ",");
				if(line.size() != 5)
				{
					cout << "Error! Problem with element = " << line[0] << " in element set = " << name << endl;
					exit(1);
				}
				
				
				element e;
				e.type = element::C3D4;

				for(int i=1; i<5; i++)
				{
					auto n = node_renumber.find(s2i(line[i]));
					if(n == node_renumber.end())
					{
						cout << "Error! Problem with node id = " << line[i] << " in element = " << line[0] << " of element set = " << name << endl;
						exit(1);
					}
					e.nodes.push_back(n->second);
				}

				elements.push_back(e);
				elementsets[name].push_back(elementcnt);
				
				elementcnt++;
				set_elecnt++;
			}

			cout << "read " << set_elecnt << " elements." << endl;

			//restart the main file loop
			continue;
		}
		
		if(l.find("*NSET, NSET=") == 0)
		{
			//get the face set name
			string name = l.substr(12);
			//force the node set to lower
			name = str2lower(name);
			cout << "\tReading node set \"" << name << "\" ... ";

			while(true)
			{
				streamoff pos = inp.tellg();
				getline(inp, l);
				if(l.find("*") != string::npos)
				{
					inp.seekg(pos);
					break;	//finished with the node section
				}
				vector<string> line = Tokenize(l, ",");
				
				for(unsigned int i=0; i<line.size(); i++)
				{
					auto n = node_renumber.find(s2i(line[i]));
					if(n == node_renumber.end())
					{
						cout << "Error! Problem with node id = " << line[i] << " in node set = " << name << endl;
						exit(1);
					}
					nodesets[name].push_back(n->second);
				}

			}
			cout << "read " << nodesets[name].size() << " nodes." << endl;
			

			//restart the main file loop
			continue;
		}
		

	}

	inp.close();

	/*
	cout << "\tBuilding element faces ... ";
	
	//update the elementfaces 
	for(int e=0; e<elecnt; e++)
	{
		vector<face> fs = elements[e]->getfaces();
		for(int f=0; f<fs.size(); f++)
		{
			elementfaces[fs[f]].push_back(e);
		}
	}

	//create the EXTERNAL faces in the mesh map
	for(map<face, vector<int> >::iterator f = elementfaces.begin(); f != elementfaces.end(); f++)
	{
		if(f->second.size() == 1)
		{
			externalelementfaces[f->first] = f->second;
		}
	}

	cout << "found " << externalelementfaces.size() << " external and " << elementfaces.size() << " total faces." << endl;

	//Build facesets from each node set
	for(map<string, nodeset>::iterator it=nodesets.begin(); it != nodesets.end(); it++)
	{
		//it->first - set name
		//it->second - the node set
		
		string name = it->first;

		//check so we don't clobber an existing faceset
		if(facesets.count(name) != 0)
		{
			error("face set = " + name + " is named the same as an existing nodeset and must be renamed!");
		}
		
		//create a c++ set that will allow for a fast lookup
		set<int> ns;
		for(int n=0; n<it->second.size(); n++)
		{
			ns.insert(it->second[n]);
		}

		//now we need to loop through all of the externalelementfaces and test for each face if
		//all the nodes are present in ns
		for(auto f=externalelementfaces.begin(); f != externalelementfaces.end(); f++)
		{
			bool faceINnset = true;	//assume all the nodes are present

			for(int n=0; n<f->first.nodes.size(); n++)
			{
				if(ns.count(f->first.nodes[n]->id) == 0)
				{
					faceINnset = false;	//this node isn't present so set to false
					break;	//break the inner for loop
				}
			}

			if(faceINnset)
			{
				//add the face to the faceset
				facesets[name].push_back(f->first);
			}
		}

		cout << "\tCreated face set \"" << name << "\" from the node set with " << facesets[name].size() << " faces." << endl;
	}

	//we can now clear the memory heavy set of all element faces
	elementfaces.clear();

	//Create a nodeset gap from the faceset if it doesn't exist
	if(nodesets.count("gap") == 0)
	{
		set<int> nodes;

		for(size_t f=0; f<facesets["gap"].size(); f++)
		{
			for(size_t n=0; n<facesets["gap"][f].nodes.size(); n++)
			{
				nodes.insert(facesets["gap"][f].nodes[n]->id);
			}
		}
		
		nodeset gapnodes;

		for(set<int>::iterator i=nodes.begin(); i != nodes.end(); i++)
		{
			gapnodes.push_back(*i);
		}
				
		nodesets["gap"] = gapnodes;

		cout << "\tCreated node set \"gap\" from the face set with " << gapnodes.size() << " nodes." << endl;
	}

	//the final step is to remove any newly created facesets that might 'clobber' faces in the gap faceset

	//first just check that gap exists
	if(!checkFaceset("gap"))
	{
		error("a node or face set named gap MUST exist!");
	}

	//let's create a c++ set for the gap faceset for quick lookup
	set<face> gapfaceset;
	for(int f=0; f<facesets["gap"].size(); f++)
	{
		gapfaceset.insert(facesets["gap"][f]);
	}

	//now loop through all NODE created facesets (except "gap") and remove any duplicate faces
	//if a faceset is actually defined in the input, we will not correct it
	for(map<string, nodeset>::iterator it=nodesets.begin(); it != nodesets.end(); it++)
	{
		string name = it->first;

		//obviously we need to skip the gap faceset
		if(name.compare("gap") == 0)
		{
			continue;
		}
		
		//let's just keep a count of how many faces we delete for logging purposes
		int del = 0;

		for(int f = 0; f<facesets[name].size(); f++)
		{
			if(gapfaceset.count(facesets[name][f]) > 0)
			{
				facesets[name].erase(facesets[name].begin() + f);
				f--;	//becuase we removed the current element, update the counter
				del++;
			}
		}

		if(del > 0)
		{
			cout << "\tRemoving " << del << " face(s) from the " << name << " face set because they duplicate the gap face set." << endl;
		}
	}
	*/

	//assign materials to element sets
	for(auto eset = elementsets.begin(); eset!=elementsets.end(); ++eset)
	{
		auto m = options.elementset_material.find(eset->first);
		if(m == options.elementset_material.end())
		{
			cout << "ERROR: Element set ( " << eset->first << " ) doesn't have a material assigned to it!" << endl;
			exit(1);
		}
		eset->second.material = options.materials[m->second];
	}

	//assign constraints to node sets
	for(auto constraint = options.elastic_constraints.begin(); constraint!=options.elastic_constraints.end(); ++constraint)
	{
		auto nset = nodesets.find(constraint->first);
		if(nset == nodesets.end())
		{
			cout << "ERROR: Constraint node set ( " << constraint->first << " ) is not defined in the mesh!" << endl;
			exit(1);
		}
		if(constraint->second.x != 0)
		{
			nset->second.const_x = true;
		}
		if(constraint->second.y != 0)
		{
			nset->second.const_y = true;
		}
		if(constraint->second.z != 0)
		{
			nset->second.const_z = true;
		}
	}
};
void fem::writeVTK(std::string filename, const std::vector<double>& b, const std::vector<double>& x)
{
	//this could be made flexible for more element types
	# define VTK_VOL_ELM_TYPE 10
	# define VTK_SURF_ELM_TYPE 5

	cout << "Writing VTK output ... ";

	ofstream vtk(filename.c_str());
	if (!vtk.is_open()) 
	{
		cout << "Error opening " << filename << endl;
		exit(1);
	}

	// write VTK header
	vtk << 
		"# vtk DataFile Version 2.0" << endl <<
		"vtk output" << endl <<
		"ASCII" << endl <<
		"DATASET UNSTRUCTURED_GRID" << endl << 
		"POINTS " << nodecnt << " double" << endl;
	
	// ----------------------------- mesh nodes -------------------------------//
	
	for(unsigned int i = 0; i < nodes.size(); i++) 
	{
		vtk << nodes[i][0] << "\t" 
				<< nodes[i][1] << "\t" 
				<< nodes[i][2] << endl; 
	}

	vtk << endl;

	// ------------------------ elements definition -------------------------- //
	
	//NOTE I HAVE FORCED ELEMENTE NODE CNT = 4 for now. -> need to be more general!
	const int elenodecnt = 4;
	vtk << "CELLS " << elementcnt << "\t" 
		<< (1 + elenodecnt)*elementcnt << endl;

	for(int i = 0; i < elementcnt; i++) {
		vtk << elenodecnt << "\t";
		for(unsigned int j = 0; j < elements[i].nodes.size(); j++)
			vtk << nodes[elements[i].nodes[j]].nid << "\t"; // index start from 0
		vtk << endl;
	}
	
	vtk << endl;

	// --------------------------- elements type ----------------------------- //

	vtk << "CELL_TYPES " << elementcnt << endl;
			
	for (int i = 0; i < elementcnt; i++) 
	{
		vtk << VTK_VOL_ELM_TYPE << endl;
	}


	vtk << "\nPOINT_DATA " << nodecnt << endl;

	vtk << "\nVECTORS displacement float" << endl;			
	for(int i=0; i<nodecnt; i++)
	{
		for(int j=0; j<3; j++)
		{
			const int Gdof = nodes[i].DOF[j];
			if(Gdof >= 0)
			{
				vtk << x[Gdof] << "\t";
			} else {
				//for now only 0 displacement BC's are supported
				vtk << 0 << "\t";
			}
		}
		vtk << endl;
	}
	
	//Applied Loads
	vtk << "\nVECTORS loads float" << endl;
	for(int i=0; i<nodecnt; i++)
	{
		for(int j=0; j<3; j++)
		{
			const int Gdof = nodes[i].DOF[j];
			if(Gdof >= 0)
			{
				vtk << b[Gdof] << "\t";
			} else {
				//constrained nodes can't be loaded
				vtk << 0 << "\t";
			}
		}
		vtk << endl;
	}

	vtk.close();
	vtk.clear();

	cout << "done." << endl;

	/*

	// --------------- write the boundary surfaces vtk files ----------------- //

	for(int i = 0; i < mesh.nfs; i++)
	{
		std::ostringstream oss;
		oss.str("");
		oss << in.pump << "." << i << ".vtk";
		vtk.open(oss.str().c_str());
		
		// determine the index of nodes to write
		std::vector<int> nodesidx(0);
		std::vector<std::vector<int>> newlabels(mesh.faceSets[i].size);
		// loop through all the faces in the face set
		for(int j = 0; j < mesh.faceSets[i].size; j++)
		{
			// loop to all the nodes in each faces
			for(unsigned int h = 0; h < mesh.faceSets[i].nodes[j].size(); h++)
			{
				// search is the node has already been considered
				std::vector<int>::iterator it;
				it = find(nodesidx.begin(),nodesidx.end(),mesh.faceSets[i].nodes[j][h]);
				// if not, add to the node list 
				if(it == nodesidx.end())
				{
					nodesidx.push_back(mesh.faceSets[i].nodes[j][h]);
					newlabels[j].push_back((int) nodesidx.size() - 1);
				}
				else
				{
					newlabels[j].push_back((int) (it - nodesidx.begin()));
				}
			}
		}
		
		// write vtk file for this boundary surface (defined in faceSets[i])
		vtk << 
		"# vtk DataFile Version 2.0" << endl <<
		"vtk output" << endl <<
		"ASCII" << endl <<
		"DATASET UNSTRUCTURED_GRID" << endl << 
		"POINTS " << nodesidx.size() << " double" << endl;
	
		// write nodes
		for(unsigned int j = 0; j < nodesidx.size(); j++)
		{
			vtk << mesh.nodes[nodesidx[j]].x() << "\t"
				  << mesh.nodes[nodesidx[j]].y() << "\t"
					<< mesh.nodes[nodesidx[j]].z() << endl;
		}

		// write faces
		int nelements = mesh.faceSets[i].size;

		vtk << "CELLS " << nelements << "\t" 
				<< (1 + FACE_NDS)*nelements << endl;

		for(int j = 0; j < nelements; j++) 
		{
			vtk << FACE_NDS << "\t";
			for(unsigned int h = 0; h < newlabels[j].size(); h++)
				vtk << newlabels[j][h] << "\t"; // index start from 0
			vtk << endl;
		}
			
		vtk << endl;

		// write face elements definition
		vtk << "CELL_TYPES " << nelements << endl;
			
		for (int j = 0; j < nelements; j++)
			vtk << VTK_SURF_ELM_TYPE << endl;

		vtk.close();
	}

	*/
}
void fem::writeFacesetVTK(const std::string filename, const std::vector<const face*>& fs)
{
	int vtkelecnt = (int) fs.size();

	//we need to build a map of nodes in this faceset and assign a vtk node id
	map<const node *, int> vtknodesmap;

	//we also need a vector of the nodes in this faceset ordered according to the local vtk node id
	vector<const node *> vtknodes;
	
	for(auto f=fs.begin(); f!=fs.end(); f++)
	{
		for(auto n=(*f)->nodes.begin(); n != (*f)->nodes.end(); n++)
		{
			if(vtknodesmap.count((*n)) == 0)
			{
				//we need to insert the node
				vtknodesmap[*n] = (int) vtknodes.size();
				vtknodes.push_back(*n);
			}
		}
	}

	int vtknodecnt = (int) vtknodes.size();

	//now we can write the vtk file
	ofstream vtk(filename.c_str());
	if (!vtk.is_open()) 
	{
		cout << "Error opening " + filename + "!" << endl;
		exit(1);
	}

	// write VTK header
	vtk << 
		"# vtk DataFile Version 2.0" << endl <<
		"vtk output" << endl <<
		"ASCII" << endl <<
		"DATASET UNSTRUCTURED_GRID" << endl << 
		"POINTS " << vtknodecnt << " double" << endl;
	
	// ----------------------------- mesh nodes -------------------------------//
	
	
	for(unsigned int i = 0; i < vtknodes.size(); i++) 
	{
		vtk << vtknodes[i]->x() << "\t" 
				<< vtknodes[i]->y() << "\t" 
				<< vtknodes[i]->z() << endl; 
	}
	
	vtk << endl;

	// ------------------------ elements definition -------------------------- //
	
	//NOTE I HAVE FORCED ELEMENTE NODE CNT = 3 for now. -> need to be more general!
	const int elenodecnt = 3;
	vtk << "CELLS " << vtkelecnt << "\t" 
		<< (1 + elenodecnt)*vtkelecnt << endl;

	for(auto f=fs.begin(); f != fs.end(); f++)
	{
		vtk << elenodecnt << "\t";
		for(auto n=(*f)->nodes.begin(); n != (*f)->nodes.end(); n++)
		{
			vtk << vtknodesmap[*n] << "\t";
		}
		vtk << endl;
	}

	vtk << endl;

	// --------------------------- elements type ----------------------------- //

	vtk << "CELL_TYPES " << vtkelecnt << endl;
			
	for (int i = 0; i < vtkelecnt; i++) 
	{
		vtk << VTK_SURF_ELM_TYPE << endl;
	}

	vtk.close();

}
void fem::getim(const std::vector<double>::iterator x, std::vector<double>& im)
{
	const size_t gap_nodecnt = gap_nodeset->size();

	im.clear();
	im.resize(gap_nodecnt, 0);

	for(size_t i=0; i<gap_nodecnt; ++i)
	{
		point deformation(0,0,0); //node deformation
		for(int j=0; j<3; ++j)
		{
			int gdof = nodes[(*gap_nodeset)[i]].DOF[j];
			if(gdof < 0)
			{
				//this really shouldn't be possible, to constrain a dof of the gap, but oh well
				//would normally need to grab the B.C. value if we didn't just assume 0 ebc
				deformation[j] = 0;
			} else {
				deformation[j] = *(x+gdof);
			}
		}

		//only get deformation normal to the gap direction
		im[i] = deformation.dot(gap_node_norms[i]);
		
		//by convention the gap norms are actually the 'anti-normals'
		//so reverse the magnitude in the case of a normal
		if(im_direction == NORMAL)
		{
			im[i] *= -1;
		}
	}
}