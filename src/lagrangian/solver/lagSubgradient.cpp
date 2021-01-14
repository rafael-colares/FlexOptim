#include "lagSubgradient.h"


void lagSubgradient::run(){
    std::cout << "--- Subgradient was invoked ---" << std::endl;
    //formulation->displayToBeRouted();

    initialization();
    std::cout << "> Subgradient was initialized. " << std::endl;

    bool STOP = false;
    while (!STOP){
        runIteration();
        if (formulation->getStatus() != RSA::STATUS_INFEASIBLE){
            //displayMainParameters();

            updateLambda();

            updateStepSize();
            
            formulation->updateMultiplier(getStepSize());
            formulation->updateCosts();

            //fichier << "\n\n> ******************** Iteration: " << getIteration() << " ********************"<< std::endl;
            //formulation->displaySlack(fichier);
            //formulation->displayMultiplier(fichier);

            //formulation->displaySlack();
            //formulation->displayMultiplier();

            if (getLB() >= getUB() - DBL_EPSILON){ 
                formulation->setStatus(RSA::STATUS_OPTIMAL);
                STOP = true;
                setStop("Optimal");
                std::cout << "Optimal" << std::endl;
            }
            else if (getIteration() >= MAX_NB_IT){
                STOP = true;
                setStop("Max It");
                std::cout << "Max It" << std::endl;
            }
            else if (isGradientMoving() == false){
                std::cout << getStepSize() << std::endl;
                STOP = true;
                setStop("Small Step");
                std::cout << "Small Step" << std::endl;
            }
            bool alternativeStop = formulation->getInstance().getInput().getAlternativeStop();
            if(alternativeStop){
                if(getGlobalItWithoutImprovement() >= 3*MAX_NB_IT_WITHOUT_IMPROVEMENT){
                    STOP = true;
                    setStop("Alternative stop");
                    std::cout << "Alternative stop" << std::endl;
                }
            }
        }
        else{
            STOP = true;
            setStop("Infeasible");
            std::cout << "Infeasible" << std::endl;
        }
    }

}

/* Sets the initial parameters for the subgradient to run. */
void lagSubgradient::initialization(){
    setIteration(0);
    setItWithoutImprovement(0);
    setGlobalItWithoutImprovement(0);
    setStepSize(0.000);

    setLB(-__DBL_MAX__);
    //setUB(__DBL_MAX__/2);
    setUB(50000);

    formulation->setStatus(RSA::STATUS_UNKNOWN);

    initLambda();
    formulation->init();
    //std::cout << "> Initialization is done. " << std::endl;

    //fichier << "> Initialization is done. " << std::endl;
    //fichier << "> ******************** Iteration: " << getIteration() << " ********************"<< std::endl;
    
    //formulation.createGraphFile(getIteration());
    //formulation->displaySlack(fichier);
    //formulation->displayMultiplier(fichier);
}

void lagSubgradient::runIteration(){ 

    incIteration();
    //std::cout << "> Running iteration " << getIteration()  << std::endl;
    formulation->run();

    updateLB(formulation->getLagrCurrentCost());

    if(formulation->getStatus() == RSA::STATUS_FEASIBLE){
        double feasibleSolutionCost = formulation->getRealCurrentCost();
        if (feasibleSolutionCost < getUB()){
            updateUB(feasibleSolutionCost);
        }
    }
    if(getIteration()==1 || getIteration()%30 ==0){
        heuristic->run();
        double feasibleSolutionCostHeur = heuristic->getCurrentHeuristicCost();
        updateUB(feasibleSolutionCostHeur);
    }
    
}

/* Updates the lambda used in the update of step size. Lambda is halved if LB has failed to increade in some fixed number of iterations */
void lagSubgradient::updateLambda(){
    double last_lambda = getLambda();
    if (last_lambda >= 10e-3){
        if (getItWithoutImprovement() >= MAX_NB_IT_WITHOUT_IMPROVEMENT){
            setItWithoutImprovement(0);
            setLambda(last_lambda / 2.0);
        }
    }
}

/* Updates the known lower bound. */
void lagSubgradient::updateLB(double bound){
    if (bound >= getLB() + __DBL_EPSILON__){
        setLB(bound);
        setItWithoutImprovement(0);
        setGlobalItWithoutImprovement(0);
        //std::cout << "> Set LB: " << bound << std:: endl;
    }
    else{
        incItWithoutImprovement();
        incGlobalItWithoutImprovement();
        //std::cout << "> Nb iterations without improvement: " << getItWithoutImprovement() << std:: endl;
    } 
}

/* Updates the step size with the rule: lambda*(UB - Z[u])/|slack| */
void lagSubgradient::updateStepSize(){
    double numerator = getLambda()*(getUB() - formulation->getLagrCurrentCost());
    double denominator = formulation->getSlackDirectionProd();
    double new_stepSize = (numerator/denominator); 
    setStepSize(new_stepSize);
}


