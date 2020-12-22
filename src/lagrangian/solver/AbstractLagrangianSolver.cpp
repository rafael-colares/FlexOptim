#include "AbstractLagrangianSolver.h"


/****************************************************************************************/
/*										Constructor										*/
/****************************************************************************************/

AbstractLagSolver::AbstractLagSolver(const Instance &inst, const Status &s): MAX_NB_IT_WITHOUT_IMPROVEMENT(inst.getInput().getNbIterationsWithoutImprovement()), 
        MAX_NB_IT(inst.getInput().getMaxNbIterations()), INITIAL_STEPSIZE(inst.getInput().getInitialLagrangianLambda()),
        MIN_STEPSIZE(0.00001),currentStatus(s){

    lagFormulationFactory factory;
    formulation = factory.createFormulation(inst);

    heuristicFactory factoryHeuristic;
    heuristic = factoryHeuristic.createHeuristic(formulation,inst);

    std::string nom = "sortie.txt";
    std::string nom2 = "sortie2.txt";
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

void AbstractLagSolver::displayMainParameters(){
    int k = getIteration();
    std::vector<int> sizeOfField;
    sizeOfField.resize(9);
    sizeOfField[0] = 5;
    sizeOfField[1] = 8;
    sizeOfField[2] = 10;
    sizeOfField[3] = 9;
    sizeOfField[4] = 9;
    sizeOfField[5] = 6;
    sizeOfField[6] = 6;
    sizeOfField[7] = 8;
    sizeOfField[8] = 11;
    //sizeOfField[9] = 12;
    char space = ' ';

    std::vector<std::string> field;
    field.resize(9);
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
        
        for (unsigned int i = 0; i < field.size(); i++){
            field[i].resize(sizeOfField[i], space);
            std::cout << field[i] << " | ";
            fichier2 << field[i] << " | ";
        }
        fichier2 << std::endl;;
        std::cout << std::endl;
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
    
    for (unsigned int i = 0; i < field.size(); i++){
        field[i].resize(sizeOfField[i], space);
        std::cout << field[i] << " | ";
        fichier2 << field[i] << " | ";

    }
    fichier2 << std::endl;
    std::cout << std::endl;
   // displayMultiplier();
}

/****************************************************************************************/
/*										DESTRUCTOR  									*/
/****************************************************************************************/

AbstractLagSolver::~AbstractLagSolver(){
    fichier.close();
    fichier2.close();
}
