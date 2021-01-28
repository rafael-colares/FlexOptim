#include "lagSubgradient.h"

/*************************************************************************/
/*					            INITIALIZATION 		    		         */
/*************************************************************************/
/* Sets the initial parameters for the subgradient to run. */
void lagSubgradient::initialization(){
    /** Setting initial time **/
    setInitializationTime(0.0);
    setConstAuxGraphTime(0.0);
    setSolvingSubProblemTime(0.0);
    setUpdatingSlackTime(0.0);
    setUpdatingBoundsTime(0.0);
    setHeuristicBoundTime(0.0);
    setUpdatingMultipliersTime(0.0);
    setUpdatingCostsTime(0.0);
    setStoppingCriterionTime(0.0);
    setUpdatingPrimalVariablesTime(0.0);

    time.setStart(ClockTime::getTimeNow());
    
    setIteration(0);
    setItWithoutImprovement(0);
    setGlobalItWithoutImprovement(0);
    setStepSize(0.000);

    setLB(-__DBL_MAX__);
    setUB(50000); //setUB(__DBL_MAX__/2);

    Input::ObjectiveMetric chosenMetric = formulation->getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        setUB(formulation->getInstance().getMaxSlice());
        //setUB(13);
    }

    formulation->setStatus(RSA::STATUS_UNKNOWN);

    initLambda();
    formulation->init();

    setInitializationTime(time.getTimeInSecFromStart());

    setConstAuxGraphTime(formulation->getConstAuxGraphTime()); 

    std::cout << "> Initialization is done. " << std::endl;

    fichier2 << "> Initialization is done. " << std::endl;
    fichier2 << " ****** Iteration: " << getIteration() << " ****** "<< std::endl;
    //formulation->createGraphFile(getIteration());
    formulation->displaySlack(fichier2);
    formulation->displayMultiplier(fichier2);
}

/*************************************************************************/
/*					            RUNNING METHODS 		    		     */
/*************************************************************************/

void lagSubgradient::run(){
    std::cout << "--- Subgradient was invoked ---" << std::endl;
    //formulation->displayToBeRouted();

    initialization();
    std::cout << "> Subgradient was initialized. " << std::endl;

    bool STOP = false;
    while (!STOP){
        runIteration();
        if (formulation->getStatus() != RSA::STATUS_INFEASIBLE){
            displayMainParameters(fichier);

            updateLambda();
            updateStepSize();
            
            time.setStart(ClockTime::getTimeNow());
            formulation->updateMultiplier(getStepSize());
            incUpdatingMultipliersTime(time.getTimeInSecFromStart());

            time.setStart(ClockTime::getTimeNow());
            formulation->updateCosts();
            incUpdatingCostsTime(time.getTimeInSecFromStart());

            time.setStart(ClockTime::getTimeNow());
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
            incStoppingCriterionTime(time.getTimeInSecFromStart());
        }
        else{
            STOP = true;
            setStop("Infeasible");
            std::cout << "Infeasible" << std::endl;
        }
        if(STOP==true){
            setTotalTime(getGeneralTime().getTimeInSecFromStart());
        }
    }

}

void lagSubgradient::runIteration(){ 
    incIteration();

    /** Solving sub problem **/
    time.setStart(ClockTime::getTimeNow());
    formulation->run();
    incSolvingSubProblemTime(time.getTimeInSecFromStart());

    /** Updating slack **/
    time.setStart(ClockTime::getTimeNow());
    formulation->updateSlack();
    formulation->updateSlack_v2();
    formulation->updateDirection();
    incUpdatingSlackTime(time.getTimeInSecFromStart());
    
    time.setStart(ClockTime::getTimeNow());
    /** Updating feasibility **/
    if (formulation->checkFeasibility() == true){
        formulation->setStatus(RSA::STATUS_FEASIBLE);
    }

    /** Updating lower bound **/
    updateLB(formulation->getLagrCurrentCost());
    
    /** Updating upper bound **/
    if(formulation->getStatus() == RSA::STATUS_FEASIBLE){
        double feasibleSolutionCost = formulation->getRealCurrentCost();
        if (feasibleSolutionCost < getUB()){
            updateUB(feasibleSolutionCost);
        }
    }   
    incUpdatingBoundsTime(time.getTimeInSecFromStart());
    
    time.setStart(ClockTime::getTimeNow());
    Input::ObjectiveMetric chosenMetric = formulation->getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric != Input::OBJECTIVE_METRIC_8){
        if(getIteration()==1 || getIteration()%30 ==0){
            heuristic->run();
            double feasibleSolutionCostHeur = heuristic->getCurrentHeuristicCost();
            updateUB(feasibleSolutionCostHeur);
        }
    }
    incHeuristicBoundTime(time.getTimeInSecFromStart());
    
    
    fichier2 << "\n\n> ****** Iteration: " << getIteration() << " ******"<< std::endl;
    formulation->displaySlack(fichier2);
    formulation->displayMultiplier(fichier2);
    //formulation->displaySlack();
    //formulation->displayMultiplier();
}

/*************************************************************************/
/*					            UPDATE 		    		                 */
/*************************************************************************/

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


