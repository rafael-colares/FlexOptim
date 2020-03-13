#ifndef __CSVReader__h
#define __CSVReader__h
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <string>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include "Instance.h"

class CSVReader{
private:
	std::string fileName;
	std::string delimeter;
public:
	// Constructor
	CSVReader(std::string filename, std::string delm = ";") :
		fileName(filename), delimeter(delm)
	{ }

	// Function to fetch data from a CSV File
	std::vector<std::vector<std::string> > getData();
};


// Returns the substring of str between the first and last delimiters
std::string getInBetweenString(std::string str, std::string firstDelimiter, std::string lastDelimiter);

// Splits a given string into a vector by a given delimiter. Ex "1,2,3" becomes 1 2 3 if delimiter is ","
std::vector<std::string> splitBy(std::string str, std::string delimiter);


#endif