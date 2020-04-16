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

/************************************************
 * This class implements a reader of .csv files. 
 * It is used for reading the input files.
 ************************************************/
class CSVReader{
private:
	std::string fileName; 	/**< The file to be read. **/
	std::string delimeter;	/**< The delimiter used for separating data. **/
public:
	/** Constructor. **/
	CSVReader(std::string filename, std::string delm = ";") :
		fileName(filename), delimeter(delm)
	{ }

	/********************************************
	* Function to fetch data from a CSV File.
	* It goes through .csv file line by line and 
	* returns the data in a vector of vector 
	* of strings.
	*********************************************/
	std::vector<std::vector<std::string> > getData();
};


/** Returns the substring between the first and last delimiters. **/
std::string getInBetweenString(std::string str, std::string firstDelimiter, std::string lastDelimiter);

/** Splits a given string into a vector by a given delimiter. For instance, "1,2,3" becomes 1 2 3 if delimiter is ",". **/
std::vector<std::string> splitBy(std::string str, std::string delimiter);


#endif