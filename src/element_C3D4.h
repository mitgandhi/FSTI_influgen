#pragma once
#include <vector>

#include "dmatrix.h"
#include "node.h"
#include "material.h"
#include "face.h"

class C3D4
{
private:
	std::vector<const node*> nodes;
	material_type material;

	dmatrix<double> Be();
	dmatrix<double> De();
	dmatrix<double> Ke();
	void calc_detJ();
	void calc_V();
	double x(const int n);
	double y(const int n);
	double z(const int n);

	double V;
	double detJ;

public:
	C3D4(const std::vector<const node*> &Nodes, material_type mat);
	void calc(dmatrix<double> &Kele, std::vector<double>& node_mass);
	std::vector<int> ldof2gdof();
	std::vector<face> getfaces();
};
