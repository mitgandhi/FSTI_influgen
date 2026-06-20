#pragma once
#include <vector>

#include "material.h"

class elementset : public std::vector<size_t>
{
public:
	material_type material;
};
