#pragma once
#include <vector>

#include "dmatrix.h"
#include "node.h"

//Face definition
class face
{
	double triarea(const int n1, const int n2, const int n3) const;

	public:
		std::vector<const node *> nodes;
		
		void push_back(const node * n);
		bool operator<(const face & f) const;
		double area() const;
		point normal() const;
		point centroid() const;
};