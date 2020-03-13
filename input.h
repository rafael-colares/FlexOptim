#ifndef __input__h
#define __input__h

#include <iostream>
#include <fstream>
#include <string>

class Input {

private:
	const std::string PARAMETER_FILE;
	std::string linkFile;
	std::string demandFile;
	std::string assignmentFile;
	std::string onlineDemandFile;
    int nbDemandsAtOnce;

public:
	/************************************************/
	/*				   Constructors					*/
	/************************************************/
    Input(std::string parameterFile);
    Input(const Input &i);

	/************************************************/
	/*				     Getters					*/
	/************************************************/
    std::string getLinkFile() const { return linkFile; }
    std::string getDemandFile() const { return demandFile; }
    std::string getAssignmentFile() const { return assignmentFile; }
    std::string getOnlineDemandFile() const { return onlineDemandFile; }
    const std::string getParameterFile() const { return PARAMETER_FILE; }
    int getNbDemandsAtOnce() const {return nbDemandsAtOnce;}

	/************************************************/
	/*				     Methods					*/
	/************************************************/
    std::string getParameterValue(std::string pattern);
    void displayParameters();
};

#endif