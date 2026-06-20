#pragma once
#include <string>
#include <vector>

#include "node.h"

struct element
{
	enum element_types
	{
		C3D4
	};

	element_types type;
	std::vector<size_t> nodes;
};
