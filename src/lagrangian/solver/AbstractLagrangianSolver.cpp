#include "AbstractLagrangianSolver.h"


/****************************************************************************************/
/*										Constructor										*/
/****************************************************************************************/

AbstractLagSolver::AbstractLagSolver(const Instance &inst, const Status &s): MAX_NB_IT_WITHOUT_IMPROVEMENT(inst.getInput().getNbIterationsWithoutImprovement()), 
        MAX_NB_IT(inst.getInput().getMaxNbIterations()), INITIAL_STEPSIZE(inst.getInput().getInitialLagrangianLambda()),
        MIN_STEPSIZE(0.00001),currentStatus(s),time(ClockTime::getTimeNow()),generalTime(ClockTime::getTimeNow()){

    lagFormulationFactory factory;
    formulation = factory.createFormulation(inst);

    formulationConstTime = time.getTimeInSecFromStart();
    time.setStart(ClockTime::getTimeNow());

    heuristicFactory factoryHeuristic;
    heuristic = factoryHeuristic.createHeuristic(formulation,inst);

    heuristicConstTime = time.getTimeInSecFromStart();

    Input::LagFormulation formu = inst.getInput().getChosenLagFormulation();
    Input::LagMethod method = inst.getInput().getChosenLagMethod();

    int demands = formulation->getNbDemandsToBeRouted();

    std::string nom = "demands" + std::to_string(demands)+ "_method" + std::to_string(method) + "_formulation" + std::to_string(formu)+  "_sortie.txt";
    std::string nom2 =  "demands" + std::to_string(demands)+ "_method" + std::to_string(method) + "_formulation" + std::to_string(formu)+ "_sortie2.txt";
    std::cout << nom << std::endl;
    fichier.open(nom.c_str());
    fichier2.open(nom2.c_str());
}

/****************************************************************************************/
/*										METHODS  										*/
/****************************************************************************************/

/* Sets the initial lambda used for updating the step size. */
void AbstractLagSolver::initLambda(){
    setLambda(INITIAL_STEPSIZE);
    std::cout << "> Initial lambda was defined. " << std::endl;
}

bool AbstractLagSolver::isGradientMoving(){
    if (getStepSize() >= MIN_STEPSIZE){
        return true;
    }
    return false;
}

/****************************************************************************************/
/*										UPDATE  										*/
/****************************************************************************************/

/* Updates the known upper bound. */
void AbstractLagSolver::updateUB(double bound){
    if (bound < getUB()){
        setUB(bound);
        std::cout << "> Set UB: " << bound << std:: endl;
    }
}

/****************************************************************************************/
/*										DISPLAY  									*/
/****************************************************************************************/

void AbstractLagSolver::displayMainParameters(std::ostream & sortie){
    int k = getIteration();
    std::vector<int> sizeOfField;
    sizeOfField.resize(10);
    sizeOfField[0] = 5;
    sizeOfField[1] = 8;
    sizeOfField[2] = 10;
    sizeOfField[3] = 9;
    sizeOfField[4] = 9;
    sizeOfField[5] = 6;
    sizeOfField[6] = 6;
    sizeOfField[7] = 8;
    sizeOfField[8] = 11;
    sizeOfField[9] = 9;
    char space = ' ';

    std::vector<std::string> field;
    field.resize(10);
    if (k == 1){
        field[0] = "Iter";
        field[1] = "Wout Impr";
        field[2] = "LB";
        field[3] = "Lagr Cost";
        field[4] = "Real Cost";
        field[5] = "UB";
        field[6] = "Step";
        field[7] = "Lambda";
        field[8] = "Feasibility";
        field[9] = "Time";
        
        for (unsigned int i = 0; i < field.size(); i++){
            field[i].resize(sizeOfField[i], space);
            sortie << field[i] << " | ";
        }
        sortie << std::endl;;
    }
    field[0] = std::to_string(k);
    field[1] = std::to_string(getItWithoutImprovement());
    field[2] = std::to_string(getLB());
    field[3] = std::to_string(formulation->getLagrCurrentCost());
    field[4] = std::to_string(formulation->getRealCurrentCost());
    field[5] = std::to_string(getUB());
    field[6] = std::to_string(getStepSize());
    field[7] = std::to_string(getLambda());

    if (formulation->checkFeasibility()){
        field[8] = "YES";
    }
    else{
        field[8] = "NO";
    }

    field[9] = std::to_string(generalTime.getTimeInSecFromStart());
    
    for (unsigned int i = 0; i < field.size(); i++){
        field[i].resize(sizeOfField[i], space);
        sortie << field[i] << " | ";

    }
    sortie << std::endl;
}

void AbstractLagSolver::displayResults(std::ostream & sortie){
    
    std::string delimiter = ";";

    sortie << getUB() << delimiter;
    sortie << getLB() << delimiter;
    sortie << getIteration() << delimiter;
    sortie << getLambda() << delimiter;
    sortie << getStepSize() << delimiter;
    sortie << getStop() << delimiter;
    sortie << getTotalTime() << delimiter << delimiter;

    sortie << getFormulationConstTime() << delimiter;
    sortie << getHeuristicConstTime() << delimiter;
    sortie << getInitializationTime() << delimiter;
    sortie << getConstAuxGraphTime() << delimiter;
    sortie << getSolvingSubProblemTime() << delimiter;
    sortie << getUpdatingSlackTime() << delimiter;
    sortie << getUpdatingBoundsTime() << delimiter;
    sortie << getUpdatingHeuristicBoundTime() << delimiter;
    sortie << getUpdatingMultipliersTime() << delimiter;
    sortie << getUpdatingCostsTime() << delimiter;
    sortie << getStoppingCriterionTime() << delimiter;
    sortie << getUpdatingPrimalVariablesTime() << delimiter;
    sortie << std::endl;

}

/****************************************************************************************/
/*										DESTRUCTOR  									*/
/****************************************************************************************/

AbstractLagSolver::~AbstractLagSolver(){
    fichier.close();
    fichier2.close();
}
