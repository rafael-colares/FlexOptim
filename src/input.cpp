#include "input.h"

/* Default constructor initializes the object with the information contained in the parameterFile. */
Input::Input(std::string parameterFile) : PARAMETER_FILE(parameterFile){
    std::cout << "Getting linkFile." << std::endl;
    linkFile = getParameterValue("linkFile=");

    std::cout << "Getting demandFile." << std::endl;
    demandFile = getParameterValue("demandFile=");

    std::cout << "Getting assignementFile." << std::endl;
    assignmentFile = getParameterValue("assignmentFile=");

    std::cout << "Getting onlineDemandFolder." << std::endl;
    onlineDemandFolder = getParameterValue("onlineDemandFolder=");

    std::cout << "Getting output path." << std::endl;
    outputPath = getParameterValue("outputPath=");



    std::cout << "Getting number of demands to be treated at once." << std::endl;
    nbDemandsAtOnce = std::stoi(getParameterValue("nbDemandsAtOnce="));

    std::cout << "Getting number of slices in output." << std::endl;
    nbSlicesInOutputFile = std::stoi(getParameterValue("nbSlicesInOutputFile="));
    
    std::cout << "Getting method." << std::endl;
    chosenMethod = (Method) std::stoi(getParameterValue("method="));

    std::cout << "Getting preprocessing level." << std::endl;
    chosenPreprLvl = (PreprocessingLevel) std::stoi(getParameterValue("preprocessingLevel="));

    std::cout << "Getting objective." << std::endl;
    chosenObj = to_ObjectiveMetric(getParameterValue("obj="));

    std::cout << "Getting output level." << std::endl;
    chosenOutputLvl = (OutputLevel) std::stoi(getParameterValue("outputLevel="));


    std::cout << "Getting subgradient parameters." << std::endl;
    lagrangianMultiplier_zero = std::stod(getParameterValue("lagrangianMultiplier_zero="));
    lagrangianLambda_zero = std::stod(getParameterValue("lagrangianLambda_zero="));
    nbIterationsWithoutImprovement = std::stoi(getParameterValue("nbIterationsWithoutImprovement="));
    maxNbIterations = std::stoi(getParameterValue("maxNbIterations="));

    std::cout << "Populating online demand files." << std::endl;
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
                if (value.empty()){
                    std::cout << "WARNING: Found " << pattern << " but field is empty." << std::endl; 
                }
                return value;
            }
        }
        myfile.close();
    }
    else {
        std::cerr << "ERROR: Unable to open parameters file" << std::endl; 
        exit(0);
    }
    
    std::cout << "WARNING: Did not found " << pattern << " inside parameters file." << std::endl; 
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