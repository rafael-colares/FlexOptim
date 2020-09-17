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

    try
    {
        std::cout << "Getting number of demands to be treated at once." << std::endl;
        std::cout << "Value = " << getParameterValue("nbDemandsAtOnce=") << std::endl;
        nbDemandsAtOnce = std::stoi(getParameterValue("nbDemandsAtOnce="));
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << ' Tried to get nb demands at once.\n';
    }
    

    std::cout << "Getting number of slices in output." << std::endl;
    nbSlicesInOutputFile = std::stoi(getParameterValue("nbSlicesInOutputFile="));
    
    std::cout << "Getting time limits." << std::endl;
    timeLimit = to_timeLimit(getParameterValue("timeLimit="));
    globalTimeLimit = to_timeLimit(getParameterValue("globalTimeLimit="));

    std::cout << "Getting allow blocking property." << std::endl;
    allowBlocking = std::stoi(getParameterValue("allowBlocking="));

    std::cout << "Getting hop penalty." << std::endl;
    hopPenalty = std::stoi(getParameterValue("hopPenalty="));

    std::cout << "Getting node method." << std::endl;
    chosenNodeMethod = to_NodeMethod(getParameterValue("method="));
    
    std::cout << "Getting formulation." << std::endl;
    chosenFormulation = to_Formulation(getParameterValue("formulation="));

    std::cout << "Getting mip solver." << std::endl;
    chosenMipSolver = to_MIP_Solver(getParameterValue("solver="));

    std::cout << "Getting preprocessing level." << std::endl;
    chosenPreprLvl = (PreprocessingLevel) std::stoi(getParameterValue("preprocessingLevel="));

    std::cout << "Getting objective." << std::endl;
    chosenObj = to_ObjectiveMetric(getParameterValue("obj="));

    std::cout << "Getting output level." << std::endl;
    chosenOutputLvl = (OutputLevel) std::stoi(getParameterValue("outputLevel="));

    std::cout << "Getting partitioning policy parameters." << std::endl;
    chosenPartitionPolicy = to_PartitionPolicy(getParameterValue("partitionPolicy="));
    partitionLoad = std::stoi(getParameterValue("partitionLoad="));
    partitionSlice = std::stoi(getParameterValue("partitionSlice="));


    
    std::cout << "Getting number of slices in Left region." << std::endl;
    partitionSlice = std::stoi(getParameterValue("partitionSlice="));
    

    std::cout << "Getting subgradient parameters." << std::endl;
    lagrangianMultiplier_zero = std::stod(getParameterValue("lagrangianMultiplier_zero="));
    lagrangianLambda_zero = std::stod(getParameterValue("lagrangianLambda_zero="));
    nbIterationsWithoutImprovement = std::stoi(getParameterValue("nbIterationsWithoutImprovement="));
    maxNbIterations = std::stoi(getParameterValue("maxNbIterations="));

    std::cout << "Populating online demand files." << std::endl;
    populateOnlineDemandFiles();
    
    std::cout << "Finish reading input." << std::endl;
    displayMainParameters();
    checkConsistency();
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
    partitionSlice = i.getPartitionSlice();
    partitionLoad = i.getPartitionLoad();
    timeLimit = i.getIterationTimeLimit();
    globalTimeLimit = i.getOptimizationTimeLimit();
    allowBlocking = i.isBlockingAllowed();
    hopPenalty = i.getHopPenalty();

    chosenNodeMethod = i.getChosenNodeMethod();
    chosenMipSolver = i.getChosenMIPSolver();
    chosenFormulation = i.getChosenFormulation();
    chosenPreprLvl = i.getChosenPreprLvl();
    chosenObj = i.getChosenObj();
    chosenOutputLvl = i.getChosenOutputLvl();
    chosenPartitionPolicy = i.getChosenPartitionPolicy();

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
    
    if (onlineDemandFolder.empty()) {
        std::cout << "ERROR: No online demand folder was provided." << std::endl;
        exit(0);
    }
    DIR *dir;
    dirent *pdir;
    int numberOfFiles = 0;
    if((dir = opendir(onlineDemandFolder.c_str())) == NULL) {
        std::cout << "ERROR: Could not open folder '" << onlineDemandFolder << "'." << std::endl;
        exit(0);
    }
    while ( (pdir = readdir(dir)) != NULL) {
        std::string file = onlineDemandFolder + "/" + pdir->d_name;    
        if (file.back() != '.'){
            vecOnlineDemandFile.push_back(file);
            numberOfFiles++;
        }
    }
    closedir(dir);
    if (numberOfFiles == 0){
        std::cout << "WARNING: The folder '" << onlineDemandFolder << "' of online demands is empty!" << std::endl;
    }
}
/* Converts a string into an ObjectiveMetric. */
std::vector<Input::ObjectiveMetric> Input::to_ObjectiveMetric(std::string data){
    std::string delimeter = ",";
    std::vector<std::string> strVec;
	boost::algorithm::split(strVec, data, boost::is_any_of(delimeter));
    
    std::vector<Input::ObjectiveMetric> objVec;
    for (auto i = 0; i < strVec.size(); i++){
        Input::ObjectiveMetric obj;
        if (strVec[i] == "0"){
            objVec.push_back(OBJECTIVE_METRIC_0);
        }
        if (strVec[i] == "1"){
            objVec.push_back(OBJECTIVE_METRIC_1);
        }
        if (strVec[i] == "1p"){
            objVec.push_back(OBJECTIVE_METRIC_1p);
        }
        if (strVec[i] == "2"){
            objVec.push_back(OBJECTIVE_METRIC_2);
        }
        if (strVec[i] == "2p"){
            objVec.push_back(OBJECTIVE_METRIC_2p);
        }
        if (strVec[i] == "4"){
            objVec.push_back(OBJECTIVE_METRIC_4);
        }
        if (strVec[i] == "8"){
            objVec.push_back(OBJECTIVE_METRIC_8);
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
        if (policyId == 1){
            policy = NODE_METHOD_SUBGRADIENT;
            return policy;
        }
        if (policyId == 2){
            policy = NODE_METHOD_VOLUME;
            return policy;
        }
    }
    else{
        std::cout << "ERROR: A node method must be specified." << std::endl;
        exit(0);
    }
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
        if (policyId == 1){
            policy = FORMULATION_EDGE_NODE;
            return policy;
        }
    }
    else{
        std::cout << "ERROR: A formulation must be specified." << std::endl;
        exit(0);
    }
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