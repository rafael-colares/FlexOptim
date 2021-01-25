#include "lagVolume.h"

/*************************************************************************/
/*					            INITIALIZATION 		    		         */
/*************************************************************************/

/* Sets the initial parameters for the subgradient to run. */
void lagVolume::initialization(){
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
    setGreenIt(false);
    setNbRedIt(0);
    setNbYellowIt(0);
    setItWithoutImprovement(0);

    setStepSize(0.000);
    initLambda();

    setLB(-__DBL_MAX__);
    setUB(50000);

    formulation->init();
    formulation->setStatus(RSA::STATUS_UNKNOWN);

    setInitializationTime(time.getTimeInSecFromStart());

    setConstAuxGraphTime(formulation->getConstAuxGraphTime()); 

    std::cout << "> Initialization is done. " << std::endl;
    //fichier << "> Initialization is done. " << std::endl;
    //fichier << "> ******* Iteration: " << getIteration() << " ********" << std::endl;
    //formulation->displaySlack(fichier);
    //formulation->displayMultiplier(fichier);
}

/*************************************************************************/
/*					            RUNNING METHODS 		    		     */
/*************************************************************************/

void lagVolume::run(){
    std::cout << "--- Volume was invoked ---" << std::endl;
    formulation->displayToBeRouted();

    initialization();
    std::cout << "> Volume was initialized. " << std::endl;

    bool STOP = false;
    while (!STOP){
        runIteration();
        if (formulation->getStatus() != RSA::STATUS_INFEASIBLE){
            displayMainParameters(fichier);

            updateLambda();
            updateStepSize();
            
            time.setStart(ClockTime::getTimeNow());
            formulation->updateMultiplier_v2(getStepSize());
            incUpdatingMultipliersTime(time.getTimeInSecFromStart());

            time.setStart(ClockTime::getTimeNow());
            formulation->updateCosts();
            incUpdatingCostsTime(time.getTimeInSecFromStart());

            time.setStart(ClockTime::getTimeNow());
            if(getLB() >= getUB() - DBL_EPSILON){ 
                formulation->setStatus(RSA::STATUS_OPTIMAL);
                STOP = true;
                setStop("Optimal");
                std::cout << "Optimal" << std::endl;
            }
            if(getIteration() >= MAX_NB_IT){
                STOP = true;
                setStop("Max It");
                std::cout << "Max It" << std::endl;
            }
            if(formulation->getMeanSlackModule_v2()  < MIN_MODULE_VALUE){
                STOP = true;
                setStop("Small Module");
                std::cout << "Small Module" << std::endl;
            }
            //if((getIteration()>2)&&(std::abs(formulation->getPrimalObjective()-getLB())/getLB()) < MIN_DIF_OBJ_VALUE){
            //    std::cout<< "Primal Obj " << formulation->getPrimalObjective() << std::endl;
            //    std::cout<< "LB " << getLB() << std::endl;
            //    STOP = true;
            //    setStop("Small Objective Difference");
            //    std::cout << "Small Objective Difference" << std::endl;
            //}
            incStoppingCriterionTime(time.getTimeInSecFromStart());
        }
        else{
            STOP = true;
            std::cout << "Infeasible" << std::endl;
        }
        if(STOP==true){
            setTotalTime(getGeneralTime().getTimeInSecFromStart());
        }
    }

}

void lagVolume::runIteration(){ 

    incIteration();

    /* Running sub problem */
    time.setStart(ClockTime::getTimeNow());
    formulation->run();
    incSolvingSubProblemTime(time.getTimeInSecFromStart());

    /* Updating primal solution */
    time.setStart(ClockTime::getTimeNow());
    if(getIteration()==1){
        formulation->initPrimalSolution();
    }else{
        double alpha = getAlphaValue();         
        formulation->updatePrimalSolution(alpha);
    }
    incUpdatingPrimalVariablesTime(time.getTimeInSecFromStart());

    /* Updating slacks */
    time.setStart(ClockTime::getTimeNow());
    formulation->updateSlack();
    formulation->updateSlack_v2();
    formulation->updateDirection();
    incUpdatingSlackTime(time.getTimeInSecFromStart());
    
    time.setStart(ClockTime::getTimeNow());
    /* Updating feasibility */
    if (formulation->checkFeasibility() == true){
        formulation->setStatus(RSA::STATUS_FEASIBLE);
    }

    /* Updating LOWER BOUND*/
    updateLB(formulation->getLagrCurrentCost());

    /* Updating UPPER BOUND*/
    if(formulation->getStatus() == RSA::STATUS_FEASIBLE){
        double feasibleSolutionCost = formulation->getRealCurrentCost();
        if (feasibleSolutionCost < getUB()){
            updateUB(feasibleSolutionCost);
        }
    }
    incUpdatingBoundsTime(time.getTimeInSecFromStart());

    time.setStart(ClockTime::getTimeNow());
    if(getIteration()==1 || getIteration()%30 ==0){
        heuristic->run();
        double feasibleSolutionCostHeur = heuristic->getCurrentHeuristicCost();
        updateUB(feasibleSolutionCostHeur);
    }
    incHeuristicBoundTime(time.getTimeInSecFromStart());

    if(getIteration()%100==0){
        updateParamAlpha();
    }
    if(getIteration()==1){
        previous_lb = getLB();
    }

    //fichier << "\n\n> ****** Iteration: " << getIteration() << " ******"<< std::endl;
    //formulation->displaySlack(fichier);
    //formulation->displayMultiplier(fichier);
    //formulation->displaySlack();
    //formulation->displayMultiplier();
}

/*************************************************************************/
/*					            UPDATE 		    		                 */
/*************************************************************************/

/* Updates the known lower bound. */
void lagVolume::updateLB(double bound){
    if (bound >= getLB() + __DBL_EPSILON__){
        setLB(bound);
        formulation->updateStabilityCenter();
        double prod = formulation->get_prod_slack();
        if(prod < 0){ //yellow iteration
            setGreenIt(false);
            setNbRedIt(0);
            incNbYellowIt();
        }else{  // green iteration
            setGreenIt(true);
            setNbRedIt(0);
            setNbYellowIt(0);
        }   
    }
    else{  // red iteration
        setGreenIt(false);
        setNbYellowIt(0);
        incNbRedIt();
    } 
}

/* Updates the known step size. */
void lagVolume::updateStepSize(){
   double numerator = getLambda()*(getUB() - getLB());
    double denominator = formulation->getSlackModule_v2();
    double new_stepSize = (numerator/denominator); 
    setStepSize(new_stepSize);
}

/* Updates the known lambda - according to red, green, yellow iteration. */
void lagVolume::updateLambda(){
    double last_lambda = getLambda();
    if(getGreenIt()){
        setLambda(std::min(2.0,last_lambda*UPD_STEPSIZE_GREEN_YELLOW));
    }
    if(getNbYellowIt() >= MAX_NB_IT_YELLOW){
        setLambda(std::min(2.0,last_lambda*UPD_STEPSIZE_GREEN_YELLOW));
        setNbYellowIt(0);
    }
    if(last_lambda >= 0.0005){
        if(getNbRedIt() >= MAX_NB_IT_RED){
            setLambda(last_lambda*UPD_STEPSIZE_RED);
            setNbRedIt(0);
        }
    }
}

/* Updates the alpha for the convex combinational of the primal vector. */
double lagVolume::getAlphaValue(){

    double alpha;

    double prod_hv = formulation->get_prod_slack();
    double prod_vv = formulation->getSlackModule();
    double prod_hh = formulation->getSlackModule_v2();

    double aux = (MAX_ALPHA*prod_vv - prod_hv)/ (prod_vv- prod_hv);

    if(prod_vv - 2*prod_hv + prod_hh > 0){
        alpha = (prod_hh- prod_hv)/(prod_vv - 2*prod_hv + prod_hh);
    }else{
        alpha = MAX_ALPHA;
    }

    if(alpha > MAX_ALPHA)
        alpha = MAX_ALPHA;
    if(alpha < aux)
        alpha = aux;
    if(alpha > 1.0)
        alpha = MAX_ALPHA;
    if(alpha < 0.0)
        alpha = MAX_ALPHA/10;
    
    return alpha;
}

/* Update the bound for alpha */
void lagVolume::updateParamAlpha(){
    if(MAX_ALPHA >= 10e-5){
        if((getLB()-previous_lb)/previous_lb < 0.01){
            MAX_ALPHA = MAX_ALPHA/2;
        }
    }
    previous_lb = getLB();
}