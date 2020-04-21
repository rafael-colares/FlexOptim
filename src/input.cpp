#include "input.h"

/* Default constructor initializes the object with the information contained in the parameterFile. */
Input::Input(std::string parameterFile) : PARAMETER_FILE(parameterFile){
    linkFile = getParameterValue("linkFile=");
    demandFile = getParameterValue("demandFile=");
    assignmentFile = getParameterValue("assignmentFile=");
    onlineDemandFile = getParameterValue("onlineDemandFile=");
    nbDemandsAtOnce = std::stoi(getParameterValue("nbDemandsAtOnce="));
    outputPath = getParameterValue("outputPath=");
    nbSlicesInOutputFile = std::stoi(getParameterValue("nbSlicesInOutputFile="));
    
    chosenMethod = (Method) std::stoi(getParameterValue("method="));
    chosenPreprLvl = (PreprocessingLevel) std::stoi(getParameterValue("preprocessingLevel="));

    lagrangianMultiplier_zero = std::stod(getParameterValue("lagrangianMultiplier_zero="));
    lagrangianLambda_zero = std::stod(getParameterValue("lagrangianLambda_zero="));
    nbIterationsWithoutImprovement = std::stoi(getParameterValue("nbIterationsWithoutImprovement="));
    maxNbIterations = std::stoi(getParameterValue("maxNbIterations="));
    displayMainParameters();
}

/* Copy constructor. */
Input::Input(const Input &i) : PARAMETER_FILE(i.getParameterFile()){
    linkFile = i.getLinkFile();
    demandFile = i.getDemandFile();
    assignmentFile = i.getAssignmentFile();
    onlineDemandFile = i.getOnlineDemandFile();
    nbDemandsAtOnce = i.getNbDemandsAtOnce();
    outputPath = i.getOutputPath();
    nbSlicesInOutputFile = i.getnbSlicesInOutputFile();

    chosenMethod = i.getChosenMethod();
    chosenPreprLvl = i.getChosenPreprLvl();

    lagrangianMultiplier_zero = i.getInitialLagrangianMultiplier();
    lagrangianLambda_zero = i.getInitialLagrangianLambda();
    maxNbIterations = i.getMaxNbIterations();
}

/* Returns the path to the file containing all the parameters. */
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

/* Displays the main input file paths: link, demand and assignement. */
void Input::displayMainParameters(){
    std::cout << "LINK FILE: " << linkFile << std::endl;
    std::cout << "DEMAND FILE: " << demandFile << std::endl;
    std::cout << "ASSIGNMENT FILE: " << assignmentFile << std::endl;
}

