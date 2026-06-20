#pragma once
#include <string>
#include <vector>
#include <map> 

#include "material.h"

struct femoptions
{
	struct elastic_constraint
	{
		int x;
		int y;
		int z;
	};

	//general options
	std::string pump_name;
	std::string meshfile;
	double scalefactor;
	
	//Activate Inertial Relief 
	int inrel;

	//specify the RAM ammount for Pardiso solver
	int max_ram;
	
	//The list of materials / elastic constraints
	std::vector<material_type> materials;
	std::map<std::string, size_t> elementset_material;
	std::map<std::string, elastic_constraint> elastic_constraints;

	//constructor
	femoptions()
	{
		//default maxRAM to 2000 MB
		max_ram = 2000;

		//default inertia relief to 0
		inrel = 0;

		//default scale factor to 1e-3
		scalefactor = 1e-3;

		//ensure the lists are clear
		materials.clear();
		elastic_constraints.clear();
	}

	//Clear method
	void clear()
	{
		femoptions();
	}
};
