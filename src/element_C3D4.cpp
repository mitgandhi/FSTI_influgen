#include "element_C3D4.h"

C3D4::C3D4(const std::vector<const node*> &Nodes, material_type mat) : nodes(Nodes), material(mat)
{
}
dmatrix<double> C3D4::Be()
{
	dmatrix<double> B(6,12);

	const double a1 = y(2)*z(1) - y(1)*z(2) + y(1)*z(3) - y(3)*z(1) - y(2)*z(3) + y(3)*z(2);
	const double a2 = y(0)*z(2) - y(2)*z(0) - y(0)*z(3) + y(3)*z(0) + y(2)*z(3) - y(3)*z(2);
	const double a3 = y(1)*z(0) - y(0)*z(1) + y(0)*z(3) - y(3)*z(0) - y(1)*z(3) + y(3)*z(1);
	const double a4 = y(0)*z(1) - y(1)*z(0) - y(0)*z(2) + y(2)*z(0) + y(1)*z(2) - y(2)*z(1);
	const double b1 = x(1)*z(2) - x(2)*z(1) - x(1)*z(3) + x(3)*z(1) + x(2)*z(3) - x(3)*z(2);
	const double b2 = x(2)*z(0) - x(0)*z(2) + x(0)*z(3) - x(3)*z(0) - x(2)*z(3) + x(3)*z(2);
	const double b3 = x(0)*z(1) - x(1)*z(0) - x(0)*z(3) + x(3)*z(0) + x(1)*z(3) - x(3)*z(1);
	const double b4 = x(1)*z(0) - x(0)*z(1) + x(0)*z(2) - x(2)*z(0) - x(1)*z(2) + x(2)*z(1);
	const double c1 = x(2)*y(1) - x(1)*y(2) + x(1)*y(3) - x(3)*y(1) - x(2)*y(3) + x(3)*y(2);
	const double c2 = x(0)*y(2) - x(2)*y(0) - x(0)*y(3) + x(3)*y(0) + x(2)*y(3) - x(3)*y(2);
	const double c3 = x(1)*y(0) - x(0)*y(1) + x(0)*y(3) - x(3)*y(0) - x(1)*y(3) + x(3)*y(1);
	const double c4 = x(0)*y(1) - x(1)*y(0) - x(0)*y(2) + x(2)*y(0) + x(1)*y(2) - x(2)*y(1);

	B(0,0)  = a1;
	B(0,3)  = a2;
	B(0,6)  = a3;
	B(0,9)  = a4;

	B(1,1)  = b1;
	B(1,4)  = b2;
	B(1,7)  = b3;
	B(1,10) = b4;

	B(2,2)  = c1;
	B(2,5)  = c2;
	B(2,8)  = c3;
	B(2,11) = c4;

	B(3,0)  = b1;
	B(3,1)  = a1;
	B(3,3)  = b2;
	B(3,4)  = a2;
	B(3,6)  = b3;
	B(3,7)  = a3;
	B(3,9)  = b4;
	B(3,10) = a4;

	B(4,1)  = c1;
	B(4,2)  = b1;
	B(4,4)  = c2;
	B(4,5)  = b2;
	B(4,7)  = c3;
	B(4,8)  = b3;
	B(4,10) = c4;
	B(4,11) = b4;

	B(5,0)  = c1;
	B(5,2)  = a1;
	B(5,3)  = c2;
	B(5,5)  = a2;
	B(5,6)  = c3;
	B(5,8)  = a3;
	B(5,9)  = c4;
	B(5,11) = a4;

	B *= 1.0/(detJ);
	
	return B;
};
dmatrix<double> C3D4::De()
{
	dmatrix<double> D(6,6);

	D(0,0) = D(1,1) = D(2,2) = 1.0-material.nu;
	D(3,3) = D(4,4) = D(5,5) = 0.5-material.nu;
	D(0,1) = D(0,2) = material.nu;
	D(1,0) = D(1,2) = material.nu;
	D(2,0) = D(2,1) = material.nu;

	D *= material.E/((1.0+material.nu)*(1.0-2.0*material.nu));

	return D;
};
dmatrix<double> C3D4::Ke()
{
	dmatrix<double> b = Be();
	dmatrix<double> k = b.t()*De()*b;
	k *= V;

	return k;
};
void C3D4::calc_detJ()
{
	dmatrix<double> J(4,4);
	for(int n=0; n<4; n++)
	{
		J(0,n) = 1;
		for(int i=0; i<3; i++)
		{
			J(i+1,n) = nodes[n]->coord[i];
		}
	}

	detJ = J.det();
}
void C3D4::calc_V()
{
	V = 1.0/6.0*detJ;
};
//Utility functions used to compactly access the nodal x or y or z
inline double C3D4::x(const int n)
{
	return nodes[n]->coord[0];
};
inline double C3D4::y(const int n)
{
	return nodes[n]->coord[1];
};
inline double C3D4::z(const int n)
{
	return nodes[n]->coord[2];
};

void C3D4::calc(dmatrix<double> &Kele, std::vector<double>& node_mass)
{
	calc_detJ();
	calc_V();
	Kele = Ke();

	node_mass.clear();
	node_mass.resize(4, material.rho*V/4.0);
};
std::vector<int> C3D4::ldof2gdof()
{
	std::vector<int> dofmap;
	for(int n=0; n<4; ++n)
	{
		for(int dof=0; dof<nodes[n]->DOF.size(); ++dof)
		{
			dofmap.push_back(nodes[n]->DOF[dof]);
		}
	}

	return dofmap;
}
std::vector<face> C3D4::getfaces()
{
	std::vector<face> fs;
	
	//the four faces
	{
		face f;
		f.push_back(nodes[0]);
		f.push_back(nodes[1]);
		f.push_back(nodes[2]);
		fs.push_back(f);
	}

	{
		face f;
		f.push_back(nodes[0]);
		f.push_back(nodes[3]);
		f.push_back(nodes[1]);
		fs.push_back(f);
	}

	{
		face f;
		f.push_back(nodes[1]);
		f.push_back(nodes[3]);
		f.push_back(nodes[2]);
		fs.push_back(f);
	}

	{
		face f;
		f.push_back(nodes[2]);
		f.push_back(nodes[3]);
		f.push_back(nodes[0]);
		fs.push_back(f);
	}

	return fs;
}