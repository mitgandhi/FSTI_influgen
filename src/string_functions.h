#pragma once 
#include <string>
#include <sstream>
#include <vector>

//Utility functions used for text reading
namespace string_functions
{
	std::vector<std::string> Tokenize(const std::string& str,const std::string& delimiters);
	double s2d(const std::string s);
	int s2i(const std::string s);
	std::string str2lower(std::string s);
	std::string trim(std::string s);

	template <class T> std::string n2s(const T number)
	{
		std::ostringstream oss (std::ostringstream::out);
		oss.str("");
		oss << number;
		return oss.str();
	};
}
