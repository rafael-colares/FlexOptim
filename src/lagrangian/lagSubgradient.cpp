#include "lagSubgradient.h"

lagSubgradient::lagSubgradient(const Instance &inst):MAX_NB_IT_WITHOUT_IMPROVEMENT(inst.getInput().getNbIterationsWithoutImprovement()), 
        MAX_NB_IT(inst.getInput().getMaxNbIterations()), INITIAL_STEPSIZE(inst.getInput().getInitialLagrangianLambda()),
        MIN_STEPSIZE(0.000001), formulation(inst){}

void lagSubgradient::run(){
    std::cout << "--- Subgradient was invoked ---" << std::endl;
    formulation.displayToBeRouted();

    initialization();
    std::cout << "> Subgradient was initialized. " << std::endl;
    bool STOP = false;
    //int index = 0;
    while (!STOP){
        //index++;
        runIteration();
        if (formulation.getStatus() != RSA::STATUS_INFEASIBLE){
            displayMainParameters();

            updateLambda();
            updateStepSize();
            formulation.updateMultiplier(getStepSize());
            formulation.updateCosts();

            //if (getLB() >= getUB() - DBL_EPSILON){
            if (getLB() >= getUB() - 0.0000001){   
                formulation.setStatus(RSA::STATUS_OPTIMAL);
                STOP = true;
                std::cout << "optimal" << std::endl;
            }
            if (getIteration() >= MAX_NB_IT){
                STOP = true;
                std::cout << "max_it" << std::endl;
            }
            if (isGradientMoving() == false){
                STOP = true;
                std::cout << "small_step" << std::endl;
            }
        }
        else{
            STOP = true;
            std::cout << "infeasible" << std::endl;
        }
    }

}

/* Sets the initial parameters for the subgradient to run. */
void lagSubgradient::initialization(){
    setIteration(0);
    setItWithoutImprovement(0);

    setLB(-__DBL_MAX__);
    //setUB(__DBL_MAX__/2);
    setUB(500000);

    formulation.setStatus(RSA::STATUS_UNKNOWN);

    initLambda();
    formulation.init();
    std::cout << "> Initialization is done. " << std::endl;
}

bool lagSubgradient::isGradientMoving(){
    if (getStepSize() >= MIN_STEPSIZE){
        return true;
    }
    return false;
}

void lagSubgradient::displayMainParameters(){
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
    sizeOfField[9] = 12;
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
        field[8] = "Length Feas";
        field[9] = "Flow Feas";
        
        for (unsigned int i = 0; i < field.size(); i++){
            field[i].resize(sizeOfField[i], space);
            std::cout << field[i] << " | ";
        }
        std::cout << std::endl;
    }
    field[0] = std::to_string(k);
    field[1] = std::to_string(getItWithoutImprovement());
    field[2] = std::to_string(getLB());
    field[3] = std::to_string(formulation.getLagrCurrentCost());
    field[4] = std::to_string(formulation.getRealCurrentCost());
    field[5] = std::to_string(getUB());
    field[6] = std::to_string(getStepSize());
    field[7] = std::to_string(getLambda());
    if (formulation.checkLengthFeasibility()){
        field[8] = "YES";
    }
    else{
        field[8] = "NO";
    }
    if (formulation.checkFlowFeasibility() && formulation.checkSourceTargetFeasibility()){
        field[9] = "YES";
    }
    else{
        field[9] = "NO";
    }


    for (unsigned int i = 0; i < field.size(); i++){
        field[i].resize(sizeOfField[i], space);
        std::cout << field[i] << " | ";
    }
    std::cout << std::endl;
   // displayMultiplier();
}



/* Sets the initial lambda used for updating the step size. */
void lagSubgradient::initLambda(){
    setLambda(INITIAL_STEPSIZE);
    std::cout << "> Initial lambda was defined. " << std::endl;
}


void lagSubgradient::runIteration(){ 

    incIteration();
    //std::cout << "> Running iteration " << getIteration()  << std::endl;
    formulation.run();
    updateLB(formulation.getLagrCurrentCost());

    if(formulation.getStatus() == RSA::STATUS_FEASIBLE){
        double feasibleSolutionCost = formulation.getRealCurrentCost();
        if (feasibleSolutionCost < getUB()){
            //updateOnPath();
            updateUB(feasibleSolutionCost);
        }
    }

    formulation.ShortestPathHeuristic();
    double feasibleSolutionCostHeur = formulation.getCurrentHeuristicCost();
    updateUB(feasibleSolutionCostHeur);
    
}

/* Updates the known lower bound. */
void lagSubgradient::updateLB(double bound){
    if (bound >= getLB() + __DBL_EPSILON__){
        setLB(bound);
        setItWithoutImprovement(0);
        //std::cout << "> Set LB: " << bound << std:: endl;
    }
    else{
        incItWithoutImprovement();
        //std::cout << "> Nb iterations without improvement: " << getItWithoutImprovement() << std:: endl;
    } 
}

/* Updates the known upper bound. */
void lagSubgradient::updateUB(double bound){
    if (bound < getUB()){
        setUB(bound);
        std::cout << "> Set UB: " << bound << std:: endl;
    }
}


/* Updates the lambda used in the update of step size. Lambda is halved if LB has failed to increade in some fixed number of iterations */
void lagSubgradient::updateLambda(){
    double last_lambda = getLambda();
    if (getItWithoutImprovement() >= MAX_NB_IT_WITHOUT_IMPROVEMENT){
        setItWithoutImprovement(0);
        setLambda(last_lambda / 2.0);
    }
}

/* Updates the step size with the rule: lambda*(UB - Z[u])/|slack| */
void lagSubgradient::updateStepSize(){
    double numerator = getLambda()*(getUB() - formulation.getLagrCurrentCost());
    double denominator = formulation.getSlackModule();
    //std::cout << getUB() << std::endl;
    //std::cout << numerator << std::endl;
    //std::cout << denominator << std::endl;
    //std::cout << (numerator/denominator) << std::endl;
    double new_stepSize = (numerator/denominator); 
    setStepSize(new_stepSize);
}