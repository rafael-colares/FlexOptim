#include "input.h"

/* Default constructor initializes the object with the information contained in the parameterFile. */
Input::Input(std::string parameterFile) : PARAMETER_FILE(parameterFile){
    std::cout << "Getting input file paths..." << std::endl;

    topologyFile = getParameterValue("topologyFile=");
    initialMappingDemandFile = getParameterValue("initialMappingDemandFile=");
    initialMappingAssignmentFile = getParameterValue("initialMappingAssignmentFile=");
    demandToBeRoutedFolder = getParameterValue("demandToBeRoutedFolder=");


    std::cout << "Getting GNPY parameters..." << std::endl;

    GNPY_activation = std::stoi(getParameterValue("GNPY_activation="));
    if (isGNPYEnabled()){
        GNPY_topologyFile = getParameterValue("GNPY_topologyFile=");
        GNPY_equipmentFile = getParameterValue("GNPY_equipmentFile=");
    }


    std::cout << "Getting formulation parameters..." << std::endl;

    nbDemandsAtOnce = std::stoi(getParameterValue("nbDemandsAtOnce="));
    chosenObj = to_ObjectiveMetric(getParameterValue("obj="));
    allowBlocking = std::stoi(getParameterValue("allowBlocking="));
    linearRelaxation = std::stoi(getParameterValue("linearRelaxation="));
    hopPenalty = std::stoi(getParameterValue("hopPenalty="));
    chosenFormulation = to_Formulation(getParameterValue("formulation="));
    chosenPartitionPolicy = to_PartitionPolicy(getParameterValue("partitionPolicy="));
    partitionLoad = std::stoi(getParameterValue("partitionLoad="));
    partitionSlice = std::stoi(getParameterValue("partitionSlice="));


    std::cout << "Getting optimization parameters..." << std::endl;

    chosenMipSolver = to_MIP_Solver(getParameterValue("solver="));
    chosenNodeMethod = to_NodeMethod(getParameterValue("method="));
    chosenPreprLvl = (PreprocessingLevel) std::stoi(getParameterValue("preprocessingLevel="));


    std::cout << "Getting execution parameters..." << std::endl;

    outputPath = getParameterValue("outputPath=");
    chosenOutputLvl = (OutputLevel) std::stoi(getParameterValue("outputLevel="));
    nbSlicesInOutputFile = std::stoi(getParameterValue("nbSlicesInOutputFile="));
    timeLimit = to_timeLimit(getParameterValue("timeLimit="));
    globalTimeLimit = to_timeLimit(getParameterValue("globalTimeLimit="));
    

    std::cout << "Getting subgradient parameters..." << std::endl;
    lagrangianMultiplier_zero = std::stod(getParameterValue("lagrangianMultiplier_zero="));
    lagrangianLambda_zero = std::stod(getParameterValue("lagrangianLambda_zero="));
    nbIterationsWithoutImprovement = std::stoi(getParameterValue("nbIterationsWithoutImprovement="));
    maxNbIterations = std::stoi(getParameterValue("maxNbIterations="));

    std::cout << "Populating online demand files..." << std::endl;
    populateOnlineDemandFiles();
    
    std::cout << "Finish reading input." << std::endl;
    displayMainParameters();
    checkConsistency();
}

/* Copy constructor. */
Input::Input(const Input &i) : PARAMETER_FILE(i.getParameterFile()){
    topologyFile = i.getTopologyFile();
    initialMappingDemandFile = i.getInitialMappingDemandFile();
    initialMappingAssignmentFile = i.getInitialMappingAssignmentFile();
    demandToBeRoutedFolder = i.getDemandToBeRoutedFolder();
    demandToBeRoutedFile = i.getDemandToBeRoutedFiles();

    GNPY_activation = i.isGNPYEnabled();
    if (isGNPYEnabled()){
        GNPY_topologyFile = i.getGNPYTopologyFile();
        GNPY_equipmentFile = i.getGNPYEquipmentFile();
    }

    nbDemandsAtOnce = i.getNbDemandsAtOnce();
    chosenObj = i.getChosenObj();
    allowBlocking = i.isBlockingAllowed();
    linearRelaxation = i.isRelaxed();
    hopPenalty = i.getHopPenalty();
    chosenFormulation = i.getChosenFormulation();
    chosenPartitionPolicy = i.getChosenPartitionPolicy();
    partitionSlice = i.getPartitionSlice();
    partitionLoad = i.getPartitionLoad();

    chosenMipSolver = i.getChosenMIPSolver();
    chosenNodeMethod = i.getChosenNodeMethod();
    chosenPreprLvl = i.getChosenPreprLvl();

    outputPath = i.getOutputPath();
    chosenOutputLvl = i.getChosenOutputLvl();
    nbSlicesInOutputFile = i.getnbSlicesInOutputFile();
    timeLimit = i.getIterationTimeLimit();
    globalTimeLimit = i.getOptimizationTimeLimit();

    lagrangianMultiplier_zero = i.getInitialLagrangianMultiplier();
    lagrangianLambda_zero = i.getInitialLagrangianLambda();
    maxNbIterations = i.getMaxNbIterations();
    nbIterationsWithoutImprovement = i.getNbIterationsWithoutImprovement();
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
                    std::cout << "WARNING: Field '" << pattern << "' is empty." << std::endl; 
                }
                return value;
            }
        }
        myfile.close();
    }
    else {
        std::cerr << "ERROR: Unable to open parameters file '" << PARAMETER_FILE << "'." << std::endl; 
        exit(0);
    }
    std::cout << "WARNING: Did not found '" << pattern << "' inside parameters file." << std::endl; 
    return value;
}

void Input::populateOnlineDemandFiles(){
    
    if (demandToBeRoutedFolder.empty()) {
        std::cout << "ERROR: No online demand folder was provided." << std::endl;
        exit(0);
    }
    DIR *dir;
    dirent *pdir;
    int numberOfFiles = 0;
    if((dir = opendir(demandToBeRoutedFolder.c_str())) == NULL) {
        std::cout << "ERROR: Could not open folder '" << demandToBeRoutedFolder << "'." << std::endl;
        exit(0);
    }
    while ( (pdir = readdir(dir)) != NULL) {
        std::string file = demandToBeRoutedFolder + "/" + pdir->d_name;    
        if (file.back() != '.'){
            demandToBeRoutedFile.push_back(file);
            numberOfFiles++;
        }
    }
    closedir(dir);
    if (numberOfFiles == 0){
        std::cout << "WARNING: The folder '" << demandToBeRoutedFolder << "' of online demands is empty!" << std::endl;
    }
}


/* Returns the name of the objective. */
std::string Input::getObjName(ObjectiveMetric obj) const{
    std::string name = "";
    switch (obj)
    {
    case OBJECTIVE_METRIC_0:
        name = "obj_0";
        break;
    case OBJECTIVE_METRIC_1:
        name = "obj_1";
        break;
    case OBJECTIVE_METRIC_1p:
        name = "obj_1p";
        break;
    case OBJECTIVE_METRIC_2:
        name = "obj_2";
        break;
    case OBJECTIVE_METRIC_2p:
        name = "obj_2p";
        break;
    case OBJECTIVE_METRIC_4:
        name = "obj_4";
        break;
    case OBJECTIVE_METRIC_8:
        name = "obj_8";
        break;
    default:
        std::cout << "ERROR: Unknown objective." << std::endl;
        exit(0);
        break;
    }
    return name;
}

/* Converts a string into an ObjectiveMetric. */
std::vector<Input::ObjectiveMetric> Input::to_ObjectiveMetric(std::string data){
    std::string delimeter = ",";
    std::vector<std::string> strVec;
	boost::algorithm::split(strVec, data, boost::is_any_of(delimeter));
    
    std::vector<Input::ObjectiveMetric> objVec;
    for (unsigned int i = 0; i < strVec.size(); i++){
        if (strVec[i] == "0"){
            objVec.push_back(OBJECTIVE_METRIC_0);
        }
        else if (strVec[i] == "1"){
            objVec.push_back(OBJECTIVE_METRIC_1);
        }
        else if (strVec[i] == "1p"){
            objVec.push_back(OBJECTIVE_METRIC_1p);
        }
        else if (strVec[i] == "2"){
            objVec.push_back(OBJECTIVE_METRIC_2);
        }
        else if (strVec[i] == "2p"){
            objVec.push_back(OBJECTIVE_METRIC_2p);
        }
        else if (strVec[i] == "4"){
            objVec.push_back(OBJECTIVE_METRIC_4);
        }
        else if (strVec[i] == "8"){
            objVec.push_back(OBJECTIVE_METRIC_8);
        }
        else{
            std::cout << "ERROR: Invalid objective metric." << std::endl;
            exit(0);
        }
    }
    return objVec;
}
/* Converts a string into a PartitionPolicy. */
Input::PartitionPolicy Input::to_PartitionPolicy(std::string data){
    Input::PartitionPolicy policy;
    if (!data.empty()){
        int policyId = std::stoi(data);
        if (policyId == 1){
            policy = PARTITION_POLICY_SOFT;
            return policy;
        }
        if (policyId == 2){
            policy = PARTITION_POLICY_HARD;
            return policy;
        }
    }
    policy = PARTITION_POLICY_NO;
    return policy;
}


/* Converts a string into a NodeMethod. */
Input::NodeMethod Input::to_NodeMethod(std::string data){
    Input::NodeMethod policy;
    if (!data.empty()){
        int policyId = std::stoi(data);
        if (policyId == 0){
            policy = NODE_METHOD_LINEAR_RELAX;
            return policy;
        }
        else if (policyId == 1){
            policy = NODE_METHOD_SUBGRADIENT;
            return policy;
        }
        else if (policyId == 2){
            policy = NODE_METHOD_VOLUME;
            return policy;
        }
        else{
            std::cout << "ERROR: Invalid node method." << std::endl;
            exit(0);
        }
    }
    else{
        std::cout << "ERROR: A node method must be specified." << std::endl;
        exit(0);
    }
    return policy;
}

/* Converts a string into a Formulation. */
Input::Formulation Input::to_Formulation(std::string data){
    Input::Formulation policy;
    if (!data.empty()){
        int policyId = std::stoi(data);
        if (policyId == 0){
            policy = FORMULATION_FLOW;
            return policy;
        }
        else if (policyId == 1){
            policy = FORMULATION_EDGE_NODE;
            return policy;
        }
        else{
            std::cout << "ERROR: Invalid formulation." << std::endl;
            exit(0);
        }
    }
    else{
        std::cout << "ERROR: A formulation must be specified." << std::endl;
        exit(0);
    }
    return policy;
}

/* Converts a string into a Formulation. */
Input::MIP_Solver Input::to_MIP_Solver(std::string data){
    Input::MIP_Solver policy;
    if (!data.empty()){
        int policyId = std::stoi(data);
        switch (policyId)
        {
        case 0: {
            policy = MIP_SOLVER_CPLEX;
            return policy;
            break;
        }
        case 1: {
            policy = MIP_SOLVER_CBC;
            std::cout << "ERROR: MIP_Solve=1 but CBC still needs to be implemented." << std::endl;
            exit(0);
            return policy;
            break;
        }
        case 2: {
            policy = MIP_SOLVER_GUROBI;
            std::cout << "ERROR: MIP_Solve=2 but Gurobi still needs to be implemented." << std::endl;
            exit(0);
            return policy;
            break;
        }

        default:
            std::cout << "ERROR: Invalid MIP_SOLVER." << std::endl;
            exit(0);
            break;
        }
    }
    else{
        std::cout << "ERROR: A solver must be specified." << std::endl;
        exit(0);
    }
}

void Input::checkConsistency(){
    if (getChosenMIPSolver() == MIP_SOLVER_GUROBI){
        std::cout << "ERROR: MIP_Solver Gurobi has been chosen but still needs to be implemented." << std::endl;
        exit(0);
    }
    if (getChosenMIPSolver() == MIP_SOLVER_CBC){
        std::cout << "ERROR: MIP_Solver CBC has been chosen but still needs to be implemented." << std::endl;
        exit(0);
    }
    if (getChosenNodeMethod() != NODE_METHOD_LINEAR_RELAX && getChosenMIPSolver() != MIP_SOLVER_CBC){
        std::cout << "ERROR: Subgradient methods should only be called with CBC." << std::endl;
        exit(0);
    }
    if (getChosenNodeMethod() != NODE_METHOD_LINEAR_RELAX){
        std::cout << "ERROR: Subgradient methods chosen but still needs to be implemented." << std::endl;
        exit(0);
    }
    std::cout << "All information from input is consistent." << std::endl;
}

/* Converts a string into a time limit. */
int Input::to_timeLimit(std::string data){
    if (data.empty()){
        return INT_MAX;
    }
    return std::stoi(data);
}

/* Displays the main input file paths: link, demand and assignement. */
void Input::displayMainParameters(){
    std::cout << "TOPOLOGY FILE: " << topologyFile << std::endl;
    std::cout << "ROUTED DEMANDS FILE: " << initialMappingDemandFile << std::endl;
    std::cout << "INITIAL ASSIGNMENT FILE: " << initialMappingAssignmentFile << std::endl;
}

/****************************************************************************************/
/*										Destructor										*/
/****************************************************************************************/

/* Destructor. Clears the vector of strings. */
Input::~Input(){
    demandToBeRoutedFile.clear();
}