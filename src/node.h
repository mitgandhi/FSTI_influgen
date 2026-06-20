#pragma once
#include <array>
#include "point.h"

struct node
{
	//the node id in the global node vector
	size_t nid;

	//XYZ location
	point coord;

	//global DOF numbering
	//assume 3 DOF / node
	std::array<int, 3> DOF;

	node()
	{
		coord = point(0,0,0);
	};

	node(const point p)
	{
		coord = p;
	};

	node(const double x,const double y, const double z)
	{
		coord = point(x,y,z);
	};

	double operator[] (const int i) 
	{
		return coord[i];
	};

	double x() const
	{
		return coord.x();
	};

	double y() const
	{
		return coord.y();
	};

	double z() const
	{
		return coord.z();
	};
};