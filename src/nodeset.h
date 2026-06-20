#pragma once
#include <vector>

//using a vector instead of a set to maintain read order in the event
//the read order was not ascending ordered
class nodeset : public std::vector<size_t>
{
public:
	bool const_x;
	bool const_y;
	bool const_z;

	nodeset()
	{
		const_x = false;
		const_y = false;
		const_z = false;
	}
};
