#include "input.h"

Input::Input(std::string parameterFile) : PARAMETER_FILE(parameterFile){
    linkFile = getParameterValue("linkFile=");
    demandFile = getParameterValue("demandFile=");
    assignmentFile = getParameterValue("assignmentFile=");
    onlineDemandFile = getParameterValue("onlineDemandFile=");
    nbDemandsAtOnce = std::stoi(getParameterValue("nbDemandsAtOnce="));
    outputPath = getParameterValue("outputPath=");

    displayParameters();
}
Input::Input(const Input &i) : PARAMETER_FILE(i.getParameterFile()){
    linkFile = i.getLinkFile();
    demandFile = i.getDemandFile();
    assignmentFile = i.getAssignmentFile();
    onlineDemandFile = i.getOnlineDemandFile();
    nbDemandsAtOnce = i.getNbDemandsAtOnce();
    outputPath = i.getOutputPath();
}


// Searches for a pattern in the Parameter File (PF) and returns its associated value. I know it is not the most performant thing to do but PF has a reasonable size and the code becomes much clearer.
std::string Input::getParameterValue(std::string pattern){
    std::string line;
    std::string value = "";
    std::ifstream myfile (PARAMETER_FILE.c_str());
    if (myfile.is_open()) {
        while ( std::getline (myfile, line) ) {
            std::size_t pos = line.find(pattern);
            if (pos != std::string::npos){
                value = line.substr(pos + pattern.size());
                value.pop_back();
                return value;
            }
        }
        myfile.close();
    }
    else {
        std::cout << "ERROR: Unable to open parameters file" << std::endl; 
        abort();
    }
    return value;
}

void Input::displayParameters(){
    std::cout << "LINK FILE: " << linkFile << std::endl;
    std::cout << "DEMAND FILE: " << demandFile << std::endl;
    std::cout << "ASSIGNMENT FILE: " << assignmentFile << std::endl;
}

