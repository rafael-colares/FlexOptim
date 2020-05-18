#include "input.h"

/* Default constructor initializes the object with the information contained in the parameterFile. */
Input::Input(std::string parameterFile) : PARAMETER_FILE(parameterFile){
    linkFile = getParameterValue("linkFile=");
    std::cout << "Reading linkFile." << std::endl;
    demandFile = getParameterValue("demandFile=");
    std::cout << "Reading demandFile." << std::endl;
    assignmentFile = getParameterValue("assignmentFile=");
    std::cout << "Reading assignementFile." << std::endl;
    onlineDemandFolder = getParameterValue("onlineDemandFolder=");
    std::cout << "Reading onlineDemandFolder." << std::endl;
    outputPath = getParameterValue("outputPath=");
    std::cout << "Reading outputPath." << std::endl;

    nbDemandsAtOnce = std::stoi(getParameterValue("nbDemandsAtOnce="));
    std::cout << "Reading nbDemandsAtOnce." << std::endl;
    nbSlicesInOutputFile = std::stoi(getParameterValue("nbSlicesInOutputFile="));
    std::cout << "Reading nbSlicesInOutput." << std::endl;
    
    chosenMethod = (Method) std::stoi(getParameterValue("method="));
    std::cout << "Reading method." << std::endl;
    chosenPreprLvl = (PreprocessingLevel) std::stoi(getParameterValue("preprocessingLevel="));
    std::cout << "Reading preprocessingLevel." << std::endl;
    chosenObj = to_ObjectiveMetric(getParameterValue("obj="));
    std::cout << "Reading obj." << std::endl;
    chosenOutputLvl = (OutputLevel) std::stoi(getParameterValue("outputLevel="));
    std::cout << "Reading outputlvl." << std::endl;

    lagrangianMultiplier_zero = std::stod(getParameterValue("lagrangianMultiplier_zero="));
    std::cout << "Reading lagr1." << std::endl;
    lagrangianLambda_zero = std::stod(getParameterValue("lagrangianLambda_zero="));
    std::cout << "Reading lagr2." << std::endl;
    nbIterationsWithoutImprovement = std::stoi(getParameterValue("nbIterationsWithoutImprovement="));
    std::cout << "Reading lagr3." << std::endl;
    maxNbIterations = std::stoi(getParameterValue("maxNbIterations="));
    std::cout << "Reading lagr4." << std::endl;

    std::cout << "Populating onlineDemands." << std::endl;
    if (!onlineDemandFolder.empty()) {
        populateOnlineDemandFiles();
    }
    
    
    std::cout << "Finish input." << std::endl;
    displayMainParameters();
}

/* Copy constructor. */
Input::Input(const Input &i) : PARAMETER_FILE(i.getParameterFile()){
    linkFile = i.getLinkFile();
    demandFile = i.getDemandFile();
    assignmentFile = i.getAssignmentFile();
    onlineDemandFolder = i.getOnlineDemandFolder();
    vecOnlineDemandFile = i.getOnlineDemandFiles();
    outputPath = i.getOutputPath();

    nbDemandsAtOnce = i.getNbDemandsAtOnce();
    nbSlicesInOutputFile = i.getnbSlicesInOutputFile();

    chosenMethod = i.getChosenMethod();
    chosenPreprLvl = i.getChosenPreprLvl();
    chosenObj = i.getChosenObj();
    chosenOutputLvl = i.getChosenOutputLvl();

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

/****************************************************************************************/
/*										Destructor										*/
/****************************************************************************************/

/* Destructor. Clears the vector of strings. */
Input::~Input(){
    vecOnlineDemandFile.clear();
}