#include "input.h"

/* Default constructor initializes the object with the information contained in the parameterFile. */
Input::Input(std::string parameterFile) : PARAMETER_FILE(parameterFile){
    linkFile = getParameterValue("linkFile=");
    demandFile = getParameterValue("demandFile=");
    assignmentFile = getParameterValue("assignmentFile=");
    onlineDemandFolder = getParameterValue("onlineDemandFolder=");
    nbDemandsAtOnce = std::stoi(getParameterValue("nbDemandsAtOnce="));
    outputPath = getParameterValue("outputPath=");
    nbSlicesInOutputFile = std::stoi(getParameterValue("nbSlicesInOutputFile="));
    
    chosenMethod = (Method) std::stoi(getParameterValue("method="));
    chosenPreprLvl = (PreprocessingLevel) std::stoi(getParameterValue("preprocessingLevel="));
    chosenObj = to_ObjectiveMetric(getParameterValue("obj="));

    lagrangianMultiplier_zero = std::stod(getParameterValue("lagrangianMultiplier_zero="));
    lagrangianLambda_zero = std::stod(getParameterValue("lagrangianLambda_zero="));
    nbIterationsWithoutImprovement = std::stoi(getParameterValue("nbIterationsWithoutImprovement="));
    maxNbIterations = std::stoi(getParameterValue("maxNbIterations="));

    if (!onlineDemandFolder.empty()) {
        populateOnlineDemandFiles();
    }
    displayMainParameters();
}

/* Copy constructor. */
Input::Input(const Input &i) : PARAMETER_FILE(i.getParameterFile()){
    linkFile = i.getLinkFile();
    demandFile = i.getDemandFile();
    assignmentFile = i.getAssignmentFile();
    onlineDemandFolder = i.getOnlineDemandFolder();
    vecOnlineDemandFile = i.getOnlineDemandFiles();
    nbDemandsAtOnce = i.getNbDemandsAtOnce();
    outputPath = i.getOutputPath();
    nbSlicesInOutputFile = i.getnbSlicesInOutputFile();

    chosenMethod = i.getChosenMethod();
    chosenPreprLvl = i.getChosenPreprLvl();
    chosenObj = i.getChosenObj();

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

void Input::populateOnlineDemandFiles(){
    DIR *dir;
    dirent *pdir;
    dir = opendir(onlineDemandFolder.c_str());
    while ( (pdir = readdir(dir)) != NULL) {
        std::string file = onlineDemandFolder + "/" + pdir->d_name;    
        if (file.back() != '.'){
            vecOnlineDemandFile.push_back(file);
        }
    }
    closedir(dir);
}
/* Converts a string into an ObjectiveMetric. */
Input::ObjectiveMetric Input::to_ObjectiveMetric(std::string data){
    Input::ObjectiveMetric obj = (ObjectiveMetric) std::stoi(data);
    if (data == "1"){
        obj = OBJECTIVE_METRIC_1;
    }
    if (data == "1p"){
        obj = OBJECTIVE_METRIC_1p;
    }
    return obj;
}

/* Displays the main input file paths: link, demand and assignement. */
void Input::displayMainParameters(){
    std::cout << "LINK FILE: " << linkFile << std::endl;
    std::cout << "DEMAND FILE: " << demandFile << std::endl;
    std::cout << "ASSIGNMENT FILE: " << assignmentFile << std::endl;
}