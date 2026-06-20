#include "face.h"
#include <algorithm>

using namespace std;

double face::triarea(const int n1, const int n2, const int n3) const
{
	//Pythagorean sum of the areas of the respective projections on the three principal planes
	double area2 = 0;
	dmatrix<double> m (3,3);
	m(2,0) = 1; m(2,1) = 1; m(2,2) = 1;

	for(int j=0; j<3; j++)
	{
		m(0,0) = nodes[n1]->coord[j%3];
		m(0,1) = nodes[n2]->coord[j%3];
		m(0,2) = nodes[n3]->coord[j%3];
		m(1,0) = nodes[n1]->coord[(j+1)%3];
		m(1,1) = nodes[n2]->coord[(j+1)%3];
		m(1,2) = nodes[n3]->coord[(j+1)%3];
			
		area2 += pow(m.det(),2);
	}
	return 0.5*pow(area2,0.5);
};

void face::push_back(const node * n)
{
	nodes.push_back(n);
};

bool face::operator<(const face & f) const
{
	//first check node sizes
	if(nodes.size() < f.nodes.size())
	{
		return true;
	} else if (nodes.size() > f.nodes.size())
	{
		return false;
	}

	//we are going to do the rest of this comparison based on the pointer address of node *

	const int nn = int(nodes.size());

	//build a sortable vector
	vector<const node *> ns(nodes);
	vector<const node *> fns(f.nodes);

	//sort the nodal ids
	sort(ns.begin(), ns.end());
	sort(fns.begin(), fns.end());

	//now compare
	for(int n=0; n<nn; n++)
	{
		if(ns[n] < fns[n])
		{
			return true;
		} else if (ns[n] > fns[n])
		{
			return false;
		}
	}

	//must be equal, so not <
	return false;
};

double face::area() const
{
	//This should work for triangles and quads
	if(nodes.size() == 3)
	{
		return triarea(0,1,2);				
	} else if (nodes.size() == 4)
	{
		return triarea(0,1,2)+triarea(2,3,0);
	}

	return 0;
};

point face::normal() const
{
	//Returns a unit vector normal to the face
	if(nodes.size() == 3)
	{
		const point v1 = nodes[1]->coord - nodes[0]->coord;
		const point v2 = nodes[2]->coord - nodes[0]->coord;
		point c = point::cross(v1,v2);
		c *= 1.0/c.norm();
		return c;
	} else if (nodes.size() == 4)
	{
		//Size a quad can be skew, I am just going to take three of the four nodes.
		const point v1 = nodes[1]->coord - nodes[0]->coord;
		const point v2 = nodes[3]->coord - nodes[0]->coord;
		point c = point::cross(v1,v2);
		c *= 1.0/c.norm();
		return c;
	}

	return point();
};

point face::centroid() const
{
	point c(0,0,0);

	//Return a 'centroid' point
	for(int n=0; n<nodes.size(); n++)
	{
		c = c+nodes[n]->coord;
	}

	c *= (1.0/double(nodes.size()));

	return c;
}
