#include "string_functions.h"
#include <sstream>
#include <algorithm>
#include <cctype>
#include <locale>

using namespace std;

vector<string> string_functions::Tokenize(const string& str,const string& delimiters)
{
	vector<string> tokens;

    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    string::size_type pos = str.find_first_of(delimiters, lastPos);

	while (string::npos != pos || string::npos != lastPos)
    {
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiters, pos);
        pos = str.find_first_of(delimiters, lastPos);
    }

	return tokens;
};
double string_functions::s2d(const string s)
{
	istringstream iss (s);
	double n;
	iss >> n;
	return n;
};
int string_functions::s2i(const string s)
{
	istringstream iss (s);
	int n;
	iss >> n;
	return n;
};
std::string string_functions::str2lower(std::string s)
{
	std::transform(s.begin(), s.end(), s.begin(), ::tolower);
	return s;
};

static std::string &ltrim(std::string &s)
{
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) { return !std::isspace(ch); }));
	return s;
}

static std::string &rtrim(std::string &s)
{
	s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(), s.end());
	return s;
}

std::string string_functions::trim(std::string s)
{
	return ltrim(rtrim(s));
};


