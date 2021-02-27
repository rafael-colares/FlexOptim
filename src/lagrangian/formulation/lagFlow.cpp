#include "lagFlow.h"

void lagFlow::getDualSolution(double *rowprice){
    int notComputedMultipliers = (getNbDemandsToBeRouted()*instance.getNbNodes()) + getNbDemandsToBeRouted();
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){  
        notComputedMultipliers+=(countNodes(*vecGraph[d])-2);
    }
    int nbDemands = getNbDemandsToBeRouted();
    std::copy(lagrangianMultiplierLength.begin(),lagrangianMultiplierLength.end(),rowprice+notComputedMultipliers);
    notComputedMultipliers = notComputedMultipliers+nbDemands;

    int nbSlicesEdges;
    for (int i = 0; i < instance.getNbEdges(); i++){
        nbSlicesEdges = getNbSlicesLimitFromEdge(i);
        std::copy(lagrangianMultiplierOverlap[i].begin(),lagrangianMultiplierOverlap[i].end(),rowprice+notComputedMultipliers);
        notComputedMultipliers = notComputedMultipliers + nbSlicesEdges;
    }

    if(instance.getInput().isObj8(0)){
        std::copy(lagrangianMultiplierMaxUsedSliceOverall.begin(),lagrangianMultiplierMaxUsedSliceOverall.end(),rowprice+notComputedMultipliers);
    }
}

/* ******************************************************************************************************************
*                                               INITIALIZATION METHODS
******************************************************************************************************************* */

void lagFlow::init(bool initMult){
    /** Lagrangian Values **/
    if(initMult){
        initMultipliers();
    }
    initSlacks();
    initDirection();
    initCoeff();
    initAssignmentMatrix();

    /** Time **/
    setConstAuxGraphTime(0.0);
    setUpdateVariablesTime(0.0);
    setShorstestPathTime(0.0);
    setSubstractMultipliersTime(0.0);
    setCostTime(0.0);
}

/************************************************** MULTIPLIERS ****************************************************/

/* The dual multipliers considering the source constraints, the flow constraints and the target constraints
are not defined. They are always feasible (not relaxed in the lagrangian model), so the dual is always zero. */
void lagFlow::startMultipliers(double *row,int size,int objSignal){
    int notComputedMultipliers = (getNbDemandsToBeRouted()*instance.getNbNodes()) + getNbDemandsToBeRouted();
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){  
        notComputedMultipliers+=(countNodes(*vecGraph[d])-2);
    }
    int nbDemands = getNbDemandsToBeRouted();
    std::copy(row+notComputedMultipliers,row+(notComputedMultipliers+nbDemands) ,std::back_inserter(lagrangianMultiplierLength));
    notComputedMultipliers = notComputedMultipliers+nbDemands;

    int nbSlicesEdges = 0;
    for (int i = 0; i < instance.getNbEdges(); i++){
        nbSlicesEdges = getNbSlicesLimitFromEdge(i);
        std::copy(row+notComputedMultipliers,row+(notComputedMultipliers+nbSlicesEdges) ,std::back_inserter(lagrangianMultiplierOverlap[i]));
        notComputedMultipliers = notComputedMultipliers + nbSlicesEdges;
    }

    if(instance.getInput().isObj8(0)){
        std::copy(row+notComputedMultipliers,row+(notComputedMultipliers+nbDemands) ,std::back_inserter(lagrangianMultiplierMaxUsedSliceOverall));
    }
}

/* Sets the initial lagrangian multipliers for the subgradient to run. */
void lagFlow::initMultipliers(){
    double initialMultiplier = instance.getInput().getInitialLagrangianMultiplier();
    initializeLengthMultipliers(initialMultiplier);
    initializeOverlapMultipliers(initialMultiplier);

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        initializeMaxUsedSliceOverallMultipliers(initialMultiplier);
        //initializeMaxUsedSliceOverallAuxMultipliers(initialMultiplier);
        //initializeMaxUsedSliceOverall2Multipliers(initialMultiplier);
        //initializeMaxUsedSliceOverall3Multipliers(initialMultiplier);
    }
    std::cout << "> Initial Lagrangian multipliers were defined. " << std::endl;
}

/************************************************ STABILITY CENTER **************************************************/

void lagFlow::initStabilityCenter(){
    initializeLengthSC();
    initializeOverlapSC();

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        initializeMaxUsedSliceOverallSC();
        //initializeMaxUsedSliceOverallAuxSC();
        //initializeMaxUsedSliceOverall2SC();
        //initializeMaxUsedSliceOverall3SC();
    }
    std::cout << "> Initial Stability Centers were defined. " << std::endl;
}

/****************************************************** SLACK ********************************************************/

/* Initializes the slack of relaxed constraints. */
void lagFlow::initSlacks(){
    initializeLengthSlacks();
    initializeOverlapSlacks();

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        initializeMaxUsedSliceOverallSlacks();
        //initializeMaxUsedSliceOverallAuxSlacks();
        //initializeMaxUsedSliceOverall2Slacks();
        //initializeMaxUsedSliceOverall3Slacks();
    }
    std::cout << "> Initial Slacks were defined. " << std::endl;
}

/* Resets the slack of relaxed constraints. */
void lagFlow::resetSlacks(){
    resetLengthSlacks();
    resetOverlapSlacks();

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        resetMaxUsedSliceOverallSlacks();
        //resetMaxUsedSliceOverallAuxSlacks();
        //resetMaxUsedSliceOverall2Slacks();
        //resetMaxUsedSliceOverall3Slacks();
    }
}

/************************************************ SLACK PRIMAL APPROXIMATION *********************************************/

/* Initializes the slack of relaxed constraints. */
void lagFlow::initPrimalSlacks(){
    initializeLengthPrimalSlacks();
    initializeOverlapPrimalSlacks();

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        initializeMaxUsedSliceOverallPrimalSlacks();
        //initializeMaxUsedSliceOverallAuxPrimalSlacks();
        //initializeMaxUsedSliceOverall2Slacks();
        //initializeMaxUsedSliceOverall3Slacks();
    }
    std::cout << "> Initial Slacks were defined. " << std::endl;
}

/**************************************************** DIRECTION ******************************************************/

void lagFlow::initDirection(){
    initializeLengthDirection();
    initializeOverlapDirection();

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        initializeMaxUsedSliceOverallDirection();
        //initializeMaxUsedSliceOverallAuxDirection();
        //initializeMaxUsedSliceOverall2Direction();
        //initializeMaxUsedSliceOverall3Direction();
    }
    std::cout << "> Initial Direction were defined. " << std::endl;
}

/* *********************************************************************************************************************
*                                                   UPDATE METHODS
********************************************************************************************************************** */

/*************************************************** MULTIPLIERS *******************************************************/

/* Updates lagrangian multiplier with the rule: u[k+1] = u[k] + t[k]*violation */
void lagFlow::updateMultiplier(double step){
    updateLengthMultiplier(step);
    updateOverlapMultiplier(step);

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        updateMaxUsedSliceOverallMultiplier(step);
        //updateMaxUsedSliceOverallAuxMultiplier(step);
        //updateMaxUsedSliceOverall2Multiplier(step);
        //updateMaxUsedSliceOverall3Multiplier(step);
    }
}

/************************************** MULTIPLIER CONSIDERING THE STABILITY CENTER *****************************************/

void lagFlow::updateMultiplier_v2(double step){
    updateLengthMultiplier_v2(step);
    updateOverlapMultiplier_v2(step);

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        updateMaxUsedSliceOverallMultiplier_v2(step);
        //updateMaxUsedSliceOverallAuxMultiplier_v2(step);
        //updateMaxUsedSliceOverall2Multiplier_v2(step);
        //updateMaxUsedSliceOverall3Multiplier_v2(step);
    }
}

/*************************************************** STABILITY CENTER *******************************************************/

void lagFlow::updateStabilityCenter(){
    updateLengthSC();
    updateOverlapSC();

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        updateMaxUsedSliceOverallSC();
        //updateMaxUsedSliceOverallAuxSC();
        //updateMaxUsedSliceOverall2SC();
        //updateMaxUsedSliceOverall3SC();
    }
}  

/************************************************* SLACK PRIMAL APPROXIMATION ***********************************************/

void lagFlow::updatePrimalSlack(double alpha){
    updateLengthPrimalSlack(alpha);
    updateOverlapPrimalSlack(alpha);

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        updateMaxUsedSliceOverallPrimalSlack(alpha);
        //updateMaxUsedSliceOverallAuxPrimalSlack(alpha);
        //updateMaxUsedSliceOverall2PrimalSlack(alpha);
        //updateMaxUsedSliceOverall3PrimalSlack(alpha);
    }
}
        
/****************************************************** DIRECTION ***********************************************************/

void lagFlow::updateDirection(){
    double theta = getDirectionMult();
    updateLengthDirection(theta);
    updateOverlapDirection(theta);

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        updateMaxUsedSliceOverallDirection(theta);
        //updateMaxUsedSliceOverallAuxDirection(theta);
        //updateMaxUsedSliceOverall2Direction(theta);
        //updateMaxUsedSliceOverall3Direction(theta);
    }
}

/**************************************************** ASSIGNMENT MATRIX *****************************************************/

/* Updates the assignment of a demand based on the a given path. For general objective. */
void lagFlow::updateAssignment_k(int d, DijkstraCost &path, const ListDigraph::Node &SOURCE, const ListDigraph::Node &TARGET){
    std::fill(assignmentMatrix_d[d].begin(), assignmentMatrix_d[d].end(), false);
    
    ListDigraph::Node currentNode = TARGET;
    while (currentNode != SOURCE){
        ListDigraph::Arc arc = path.predArc(currentNode);
        int index = getArcIndex(arc, d);

        /* Update Assignment */
        assignmentMatrix_d[d][index] = true;
        
        /* Update Slack */
        updateLengthSlack(d,arc);
        updateOverlapSlack(d,arc);

        currentNode = path.predNode(currentNode);
    }
}

/* Updates the assignment of a demand based on the a given path. For objective 8. */
void lagFlow::updateAssignment_k(int d, DijkstraCostObj8 &path, const ListDigraph::Node &SOURCE, const ListDigraph::Node &TARGET){
    std::fill(assignmentMatrix_d[d].begin(), assignmentMatrix_d[d].end(), false);
    
    ListDigraph::Node currentNode = TARGET; 
    while (currentNode != SOURCE){
        ListDigraph::Arc arc = path.predArc(currentNode);
        int index = getArcIndex(arc, d);

        /* Update Assignment */
        assignmentMatrix_d[d][index] = true;

        /* Update Slack */
        updateLengthSlack(d,arc);
        updateOverlapSlack(d,arc);
        updateMaxUsedSliceOverallSlack(d,arc);
        //updateMaxUsedSliceOverallAuxSlack(d,arc);
        //updateMaxUsedSliceOverall2Slack(d,arc);
        //updateMaxUsedSliceOverall3Slack(d,arc);

        currentNode = path.predNode(currentNode);
    }
}

/* Updates the assignment of a demand based on the a given path. For general objective. */
void lagFlow::updateAssignment_k(int d, CostScaling<ListDigraph,int,double> &costScale, const ListDigraph::Node &SOURCE, const ListDigraph::Node &TARGET){
    std::fill(assignmentMatrix_d[d].begin(), assignmentMatrix_d[d].end(), false);

    IterableValueMap<ListDigraph,ListDigraph::Arc,double> auxiliary(*vecGraph[d]);
    costScale.flowMap(auxiliary);
    int flow = 1;

    for(IterableValueMap<ListDigraph,ListDigraph::Arc,double>::ItemIt arc(auxiliary,flow); arc != INVALID; ++arc){
        int index = getArcIndex(arc, d);

        /* Update Assignment */
        assignmentMatrix_d[d][index] = true;
        
        /* Update Slack */
        updateLengthSlack(d,arc);
        updateOverlapSlack(d,arc);
    }
}

/**************************************************** CHECK FEASIBILITY *****************************************************/

/* Checks if sub problem solution is feasible. */
bool lagFlow::checkFeasibility(){
    if (checkLengthFeasibility() == false){
        return false;
    }
    if (checkOverlapFeasibility() == false){
        return false;
    }
    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        if(checkMaxUsedSliceOverallFeasibility() == false){
            return false;
        }
    }
    return true;
}

/* Checks if primal approximation solution is feasible. */
bool lagFlow::checkFeasibility_v2(){
    if (checkLengthFeasibility_v2() == false){
        return false;
    }
    if (checkOverlapFeasibility_v2() == false){
        return false;
    }
    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        if(checkMaxUsedSliceOverallFeasibility_v2() == false){
            return false;
        }
    }
    return true;
}

bool lagFlow::checkSlacknessCondition(){
     if (checkLengthSlacknessCondition() == false){
        return false;
    }
    if (checkOverlapSlacknessCondition() == false){
        return false;
    }
    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        if(checkMaxUsedSliceOverallSlacknessCondition() == false){
            return false;
        }
    }
    return true;
}

/* *********************************************************************************************************************
*                                                         GET METHODS
********************************************************************************************************************* */

double lagFlow::getRealCostFromPath(int d, DijkstraCost &path, const ListDigraph::Node &SOURCE, const ListDigraph::Node &TARGET){
    double total = 0.0;
    ListDigraph::Node currentNode = TARGET;
    while (currentNode != SOURCE){
        ListDigraph::Arc arc = path.predArc(currentNode);
        total += (*coeff[d])[arc];
        currentNode = path.predNode(currentNode);
    }
    return total;
}

double lagFlow::getRealCostFromPath(int d, CostScaling<ListDigraph,int,double> &costScale, const ListDigraph::Node &SOURCE, const ListDigraph::Node &TARGET){
    double total = 0.0;
    IterableValueMap<ListDigraph,ListDigraph::Arc,double> auxiliary(*vecGraph[d]);
    costScale.flowMap(auxiliary);
    int flow = 1;

    for(IterableValueMap<ListDigraph,ListDigraph::Arc,double>::ItemIt arc(auxiliary,flow); arc != INVALID; ++arc){
        total += (*coeff[d])[arc];
    }
    return total;
}


/* ******************************************************* MODULES *****************************************************/

/* Returns |slack|^2 */
double lagFlow::getSlackModule(double alpha) {
    double denominator = 0.0;
    /* Length */
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double slack = alpha*(-getLengthSlack_k(d)) + (1.0 - alpha)*(-getLengthSlack_v2_k(d));
        double mult = getLengthMultiplier_k(d);
        if((alpha == -1.0) || !((slack < - DBL_EPSILON) && (mult > -DBL_EPSILON && mult < DBL_EPSILON))){ 
            denominator += std::pow(getLengthSlack_k(d),2);
        }
    }
    /* Overlap */
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            double slack = alpha*(-getOverlapSlack_k(e,s)) + (1.0 - alpha)*(-getOverlapSlack_v2_k(e,s));
            double mult = getOverlapMultiplier_k(e,s);
            if((alpha == -1.0) || !((slack < - DBL_EPSILON) && (mult > -DBL_EPSILON && mult < DBL_EPSILON))){ 
                denominator += std::pow(getOverlapSlack_k(e,s),2);
            }
        }
    }

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            double slack = alpha*(-getMaxUsedSliceOverallSlack_k(d)) + (1.0 - alpha)*(-getMaxUsedSliceOverallSlack_v2_k(d));
            double mult = getMaxUsedSliceOverallMultiplier_k(d);
            if((alpha == -1.0) || !((slack < - DBL_EPSILON) && (mult > -DBL_EPSILON && mult < DBL_EPSILON))){ 
                denominator += std::pow(getMaxUsedSliceOverallSlack_k(d),2);
            }       
        }
        /*
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            double slack = alpha*(-getMaxUsedSliceOverallAuxSlack_k(d)) + (1.0 - alpha)*(-getMaxUsedSliceOverallAuxSlack_v2_k(d));
            double mult = getMaxUsedSliceOverallAuxMultiplier_k(d);
            if((alpha == -1.0) || !((slack < - DBL_EPSILON) && (mult > -DBL_EPSILON && mult < DBL_EPSILON))){ 
                denominator += std::pow(getMaxUsedSliceOverallAuxSlack_k(d),2);
            }       
        }
        */
        /*
        for (int e = 0; e < instance.getNbEdges(); e++){
            for (int s = 0; s < getNbSlicesLimitFromEdge(e); s++){
                double slack = alpha*(-getMaxUsedSliceOverall2Slack_k(e,s)) + (1.0 - alpha)*(-getMaxUsedSliceOverall2Slack_v2_k(e,s));
                double mult = getMaxUsedSliceOverall2Multiplier_k(e,s);
                if((alpha == -1.0) || !((slack < - DBL_EPSILON) && (mult > -DBL_EPSILON && mult < DBL_EPSILON))){ 
                    denominator += std::pow(getMaxUsedSliceOverall2Slack_k(e,s),2);
                }    
            }
        }
        */
        /*
        for (int v = 0; v < instance.getNbNodes(); v++){
            for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
                double slack = alpha*(-getMaxUsedSliceOverall3Slack_k(v,s)) + (1.0 - alpha)*(-getMaxUsedSliceOverall3Slack_v2_k(v,s));
                double mult = getMaxUsedSliceOverall3Multiplier_k(v,s);
                if((alpha == -1.0) || !((slack < - DBL_EPSILON) && (mult > -DBL_EPSILON && mult < DBL_EPSILON))){ 
                    denominator += std::pow(getMaxUsedSliceOverall3Slack_k(v,s),2);
                }     
            }
        }
        */
    }
    return denominator;
}

/* Returns |slack|^2, with the slack considering the primal variables*/
double lagFlow::getSlackModule_v2(double alpha) {
    double denominator = 0.0;
    /* Length */
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double slack = alpha*(-getLengthSlack_k(d)) + (1.0 - alpha)*(-getLengthSlack_v2_k(d));
        double mult = getLengthMultiplier_k(d);
        if((alpha == -1.0) || !((slack < - DBL_EPSILON) && (mult > -DBL_EPSILON && mult < DBL_EPSILON))){ 
            denominator += std::pow(getLengthSlack_v2_k(d),2);
        }
    }
    /* Overlap */
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            double slack = alpha*(-getOverlapSlack_k(e,s)) + (1.0 - alpha)*(-getOverlapSlack_v2_k(e,s));
            double mult = getOverlapMultiplier_k(e,s);
            if((alpha == -1.0) || !((slack < - DBL_EPSILON) && (mult > -DBL_EPSILON && mult < DBL_EPSILON))){ 
                denominator += std::pow(getOverlapSlack_v2_k(e,s),2);
            }
        }
    }

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            double slack = alpha*(-getMaxUsedSliceOverallSlack_k(d)) + (1.0 - alpha)*(-getMaxUsedSliceOverallSlack_v2_k(d));
            double mult = getMaxUsedSliceOverallMultiplier_k(d);
            if((alpha == -1.0) || !((slack < - DBL_EPSILON) && (mult > -DBL_EPSILON && mult < DBL_EPSILON))){ 
                denominator += std::pow(getMaxUsedSliceOverallSlack_v2_k(d),2);
            }       
        }
        /*
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            double slack = alpha*(-getMaxUsedSliceOverallAuxSlack_k(d)) + (1.0 - alpha)*(-getMaxUsedSliceOverallAuxSlack_v2_k(d));
            double mult = getMaxUsedSliceOverallAuxMultiplier_k(d);
            if((alpha == -1.0) || !((slack < - DBL_EPSILON) && (mult > -DBL_EPSILON && mult < DBL_EPSILON))){ 
                denominator += std::pow(getMaxUsedSliceOverallAuxSlack_v2_k(d),2);
            }       
        }
        */
        /*
        for (int e = 0; e < instance.getNbEdges(); e++){
            for (int s = 0; s < getNbSlicesLimitFromEdge(e); s++){
                double slack = alpha*(-getMaxUsedSliceOverall2Slack_k(e,s)) + (1.0 - alpha)*(-getMaxUsedSliceOverall2Slack_v2_k(e,s));
                double mult = getMaxUsedSliceOverall2Multiplier_k(e,s);
                if((alpha == -1.0) || !((slack < - DBL_EPSILON) && (mult > -DBL_EPSILON && mult < DBL_EPSILON))){ 
                    denominator += std::pow(getMaxUsedSliceOverall2Slack_v2_k(e,s),2);
                }    
            }
        }
        */
        /*
        for (int v = 0; v < instance.getNbNodes(); v++){
            for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
                double slack = alpha*(-getMaxUsedSliceOverall3Slack_k(v,s)) + (1.0 - alpha)*(-getMaxUsedSliceOverall3Slack_v2_k(v,s));
                double mult = getMaxUsedSliceOverall3Multiplier_k(v,s);
                if((alpha == -1.0) || !((slack < - DBL_EPSILON) && (mult > -DBL_EPSILON && mult < DBL_EPSILON))){ 
                    denominator += std::pow(getMaxUsedSliceOverall3Slack_v2_k(v,s),2);
                }     
            }
        }
        */
    }
    return denominator;
}

/* Returns |Direction|^2 */
double lagFlow::getDirectionModule(){

    double denominator = 0.0;

    /* Length */
    denominator += std::accumulate(lengthDirection.begin(),lengthDirection.end(),0.0,[](double sum, double direction){return sum +=std::pow(direction,2);});
    
    /* Overlap */
    for (int e = 0; e < instance.getNbEdges(); e++){
        denominator += std::accumulate(overlapDirection[e].begin(),overlapDirection[e].end(),0.0,[](double sum, double direction){return sum +=std::pow(direction,2);});
    }

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        denominator += std::accumulate(maxUsedSliceOverallDirection.begin(),maxUsedSliceOverallDirection.end(),0.0,[](double sum, double direction){return sum +=std::pow(direction,2);});
        //denominator += std::accumulate(maxUsedSliceOverallAuxDirection.begin(),maxUsedSliceOverallAuxDirection.end(),0.0,[](double sum, double direction){return sum +=std::pow(direction,2);});
     
        /*
        for (int e = 0; e < instance.getNbEdges(); e++){
            denominator += std::accumulate(maxUsedSliceOverall2Direction[e].begin(),maxUsedSliceOverall2Direction[e].end(),0.0,[](double sum, double direction){return sum +=std::pow(direction,2);});
        }
        for (int v = 0; v < instance.getNbNodes(); v++){
            denominator += std::accumulate(maxUsedSliceOverall3Direction[v].begin(),maxUsedSliceOverall3Direction[v].end(),0.0,[](double sum, double direction){return sum +=std::pow(direction,2);});
        }
        */
    }
    return denominator;
}

/* Returns (Slack*Direction) */
double lagFlow::getSlackDirectionProd(){

    Input::ProjectionType projection = getInstance().getInput().getChosenProjection();
    if((projection == Input::IMPROVED) || (projection ==Input::PROJECTED)){
        return getSlackDirectionProdProjected(projection);
    }   
    return getSlackDirectionProdNormal(); 
}

/* Returns (Slack*Direction), normal direction */
double lagFlow::getSlackDirectionProdNormal(){
    double denominator = 0.0;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        denominator += getLengthSlack_k(d)*getLengthDirection_k(d);
    }
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            denominator += getOverlapSlack_k( e, s)*getOverlapDirection_k(e,s);
        }
    }

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            denominator += getMaxUsedSliceOverallSlack_k(d)*getMaxUsedSliceOverallDirection_k(d);
            //denominator += getMaxUsedSliceOverallAuxSlack_k(d)*getMaxUsedSliceOverallAuxDirection_k(d);
        }
        
        /*
        for (int e = 0; e < instance.getNbEdges(); e++){
            int sliceLimit = getNbSlicesLimitFromEdge(e);
            for (int s = 0; s < sliceLimit; s++){
                denominator += getMaxUsedSliceOverall2Slack_k(e,s)*getMaxUsedSliceOverall2Direction_k(e,s);
            }
        }

        for (int v = 0; v < instance.getNbNodes(); v++){
            for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
                denominator += getMaxUsedSliceOverall3Slack_k(v,s)*getMaxUsedSliceOverall3Direction_k(v,s);
            }
        }
        */
        
    }
    return denominator;
}

/* Returns (Slack*Direction), projected or improved direction */
double lagFlow::getSlackDirectionProdProjected(Input::ProjectionType projection){
    double denominator = 0.0;
    if(projection == Input::IMPROVED){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            if(!((-getLengthDirection_k(d) < - DBL_EPSILON) && (getLengthMultiplier_k(d) > -DBL_EPSILON && getLengthMultiplier_k(d) < DBL_EPSILON))){ // if non negative or multiplier different from zero
                denominator += getLengthSlack_k(d)*getLengthDirection_k(d);
            }
        }
        for (int e = 0; e < instance.getNbEdges(); e++){
            for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
                if(!((-getOverlapDirection_k(e,s) < - DBL_EPSILON)&&(getOverlapMultiplier_k(e,s)> -DBL_EPSILON && getOverlapMultiplier_k(e,s) < DBL_EPSILON))){ // if non negative or multiplier different from zero
                    denominator += getOverlapSlack_k( e, s)*getOverlapDirection_k(e,s);
                }
            }
        }

        Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
        if(chosenMetric == Input::OBJECTIVE_METRIC_8){
            for (int d = 0; d < getNbDemandsToBeRouted(); d++){
                if(!((-getMaxUsedSliceOverallDirection_k(d) < - DBL_EPSILON)&&(getMaxUsedSliceOverallMultiplier_k(d) > - DBL_EPSILON && getMaxUsedSliceOverallMultiplier_k(d) < DBL_EPSILON))){
                    denominator += getMaxUsedSliceOverallSlack_k(d)*getMaxUsedSliceOverallDirection_k(d);
                }
                //if(!((-getMaxUsedSliceOverallAuxDirection_k(d) < - DBL_EPSILON)&&(getMaxUsedSliceOverallAuxMultiplier_k(d) > - DBL_EPSILON && getMaxUsedSliceOverallAuxMultiplier_k(d) < DBL_EPSILON))){
                //    denominator += getMaxUsedSliceOverallAuxSlack_k(d)*getMaxUsedSliceOverallAuxDirection_k(d);
                //}
            }

            /*
            for (int e = 0; e < instance.getNbEdges(); e++){
                int sliceLimit = getNbSlicesLimitFromEdge(e);
                for (int s = 0; s < sliceLimit; s++){
                    if(!((-getMaxUsedSliceOverall2Direction_k(e,s) < - DBL_EPSILON)&&(getMaxUsedSliceOverall2Multiplier_k(e,s) > - DBL_EPSILON && getMaxUsedSliceOverall2Multiplier_k(e,s) < DBL_EPSILON))){
                        denominator += getMaxUsedSliceOverall2Slack_k(e,s)*getMaxUsedSliceOverall2Direction_k(e,s);
                    }
                }
            }

            for (int v = 0; v < instance.getNbNodes(); v++){
                for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
                    if(!((-getMaxUsedSliceOverall3Direction_k(v,s) < - DBL_EPSILON)&&(getMaxUsedSliceOverall3Multiplier_k(v,s) > - DBL_EPSILON && getMaxUsedSliceOverall3Multiplier_k(v,s) < DBL_EPSILON))){
                        denominator += getMaxUsedSliceOverall3Slack_k(v,s)*getMaxUsedSliceOverall3Direction_k(v,s);
                    }
                }
            }
            */
        }
    }else if(projection == Input::PROJECTED){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            if(!(-getLengthDirection_k(d) < - DBL_EPSILON)){ // if non negative
                denominator += getLengthSlack_k(d)*getLengthDirection_k(d);
            }
        }
        for (int e = 0; e < instance.getNbEdges(); e++){
            for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
                if(!(-getOverlapDirection_k(e,s) < - DBL_EPSILON)){ // if non negative 
                    denominator += getOverlapSlack_k( e, s)*getOverlapDirection_k(e,s);
                }
            }
        }

        Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
        if(chosenMetric == Input::OBJECTIVE_METRIC_8){
            for (int d = 0; d < getNbDemandsToBeRouted(); d++){
                if(!(-getMaxUsedSliceOverallDirection_k(d) < - DBL_EPSILON )){
                    denominator += getMaxUsedSliceOverallSlack_k(d)*getMaxUsedSliceOverallDirection_k(d);
                }
                //if(!(-getMaxUsedSliceOverallAuxDirection_k(d) < - DBL_EPSILON )){
                //    denominator += getMaxUsedSliceOverallAuxSlack_k(d)*getMaxUsedSliceOverallAuxDirection_k(d);
                //}
            }
            /*
            for (int e = 0; e < instance.getNbEdges(); e++){
                int sliceLimit = getNbSlicesLimitFromEdge(e);
                for (int s = 0; s < sliceLimit; s++){
                    if(!(-getMaxUsedSliceOverall2Direction_k(e,s) < - DBL_EPSILON)){
                        denominator += getMaxUsedSliceOverall2Slack_k(e,s)*getMaxUsedSliceOverall2Direction_k(e,s);
                    }
                }
            }

            for (int v = 0; v < instance.getNbNodes(); v++){
                for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
                    if(!(-getMaxUsedSliceOverall3Direction_k(v,s) < - DBL_EPSILON)){
                        denominator += getMaxUsedSliceOverall3Slack_k(v,s)*getMaxUsedSliceOverall3Direction_k(v,s);
                    }
                }
            }
            */
        }

    }
    return denominator;
}

/* Returns mean of |slack|, slack considering the primal variables */
double lagFlow::getMeanSlackModule_v2(){
    double module = std::accumulate(lengthSlack_v2.begin(),lengthSlack_v2.end(),0.0,[](double sum, double slack){return sum +=std::abs(slack);});
    for (int e = 0; e < instance.getNbEdges(); e++){
        module += std::accumulate(overlapSlack_v2[e].begin(),overlapSlack_v2[e].end(),0.0,[](double sum, double slack){return sum +=std::abs(slack);});
    }
    
    double numRest = getNbDemandsToBeRouted();
    for (int e = 0; e < instance.getNbEdges(); e++){
        numRest += instance.getPhysicalLinkFromIndex(e).getNbSlices();
    }

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        module += std::accumulate(maxUsedSliceOverallSlack_v2.begin(),maxUsedSliceOverallSlack_v2.end(),0.0,[](double sum, double slack){return sum +=std::abs(slack);});
        //module += std::accumulate(maxUsedSliceOverallAuxSlack_v2.begin(),maxUsedSliceOverallAuxSlack_v2.end(),0.0,[](double sum, double slack){return sum +=std::abs(slack);});

        /*
        for (int e = 0; e < instance.getNbEdges(); e++){
            module += std::accumulate(maxUsedSliceOverall2Slack_v2[e].begin(),maxUsedSliceOverall2Slack_v2[e].end(),0.0,[](double sum, double slack){return sum +=std::abs(slack);});
        }

        for (int v = 0; v < instance.getNbNodes(); v++){
            module += std::accumulate(maxUsedSliceOverall3Slack_v2[v].begin(),maxUsedSliceOverall3Slack_v2[v].end(),0.0,[](double sum, double slack){return sum +=std::abs(slack);});
        }
        */

        numRest += getNbDemandsToBeRouted();
        //numRest += getNbDemandsToBeRouted();
        /*
        for (int e = 0; e < instance.getNbEdges(); e++){
            numRest += getNbSlicesLimitFromEdge(e);
        }
        for (int e = 0; e < instance.getNbNodes(); e++){
            numRest += getNbSlicesGlobalLimit();
        }
        */
    }
    return module/numRest;
}

/* Returns (slack*slack_v2)*/ 
double lagFlow::getSlackPrimalSlackProd(double alpha){
    double denominator = 0.0;
    /* Length */
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double slack = alpha*(-getLengthSlack_k(d)) + (1.0 - alpha)*(-getLengthSlack_v2_k(d));
        double mult = getLengthMultiplier_k(d);
        if((alpha == -1.0) || !((slack < - DBL_EPSILON) && (mult > -DBL_EPSILON && mult < DBL_EPSILON))){ 
            denominator += getLengthSlack_k(d)*getLengthSlack_v2_k(d);
        }
    }
    /* Overlap */
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            double slack = alpha*(-getOverlapSlack_k(e,s)) + (1.0 - alpha)*(-getOverlapSlack_v2_k(e,s));
            double mult = getOverlapMultiplier_k(e,s);
            if((alpha == -1.0) || !((slack < - DBL_EPSILON) && (mult > -DBL_EPSILON && mult < DBL_EPSILON))){ 
                denominator += getOverlapSlack_k(e,s)*getOverlapSlack_v2_k(e,s);
            }
        }
    }

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            double slack = alpha*(-getMaxUsedSliceOverallSlack_k(d)) + (1.0 - alpha)*(-getMaxUsedSliceOverallSlack_v2_k(d));
            double mult = getMaxUsedSliceOverallMultiplier_k(d);
            if((alpha == -1.0) || !((slack < - DBL_EPSILON) && (mult > -DBL_EPSILON && mult < DBL_EPSILON))){ 
                denominator += getMaxUsedSliceOverallSlack_k(d)*getMaxUsedSliceOverallSlack_v2_k(d);
            }       
        }
        /*
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            double slack = alpha*(-getMaxUsedSliceOverallAuxSlack_k(d)) + (1.0 - alpha)*(-getMaxUsedSliceOverallAuxSlack_v2_k(d));
            double mult = getMaxUsedSliceOverallAuxMultiplier_k(d);
            if((alpha == -1.0) || !((slack < - DBL_EPSILON) && (mult > -DBL_EPSILON && mult < DBL_EPSILON))){ 
                denominator += getMaxUsedSliceOverallAuxSlack_k(d)*getMaxUsedSliceOverallAuxSlack_v2_k(d);
            }       
        }
        */
        /*
        for (int e = 0; e < instance.getNbEdges(); e++){
            for (int s = 0; s < getNbSlicesLimitFromEdge(e); s++){
                double slack = alpha*(-getMaxUsedSliceOverall2Slack_k(e,s)) + (1.0 - alpha)*(-getMaxUsedSliceOverall2Slack_v2_k(e,s));
                double mult = getMaxUsedSliceOverall2Multiplier_k(e,s);
                if((alpha == -1.0) || !((slack < - DBL_EPSILON) && (mult > -DBL_EPSILON && mult < DBL_EPSILON))){ 
                    denominator += getMaxUsedSliceOverall2Slack_k(e,s)*getMaxUsedSliceOverall2Slack_v2_k(e,s);
                }    
            }
        }
        */
        /*
        for (int v = 0; v < instance.getNbNodes(); v++){
            for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
                double slack = alpha*(-getMaxUsedSliceOverall3Slack_k(v,s)) + (1.0 - alpha)*(-getMaxUsedSliceOverall3Slack_v2_k(v,s));
                double mult = getMaxUsedSliceOverall3Multiplier_k(v,s);
                if((alpha == -1.0) || !((slack < - DBL_EPSILON) && (mult > -DBL_EPSILON && mult < DBL_EPSILON))){ 
                    denominator += getMaxUsedSliceOverall3Slack_k(v,s)*getMaxUsedSliceOverall3Slack_v2_k(v,s);
                }     
            }
        }
        */
    }
    return denominator;
}

/* Returns the physical length of the path. */
double lagFlow::getPathLength(int d, DijkstraCost &path, const ListDigraph::Node &s, const ListDigraph::Node &t){
    double pathLength = 0.0;
    ListDigraph::Node n = t;
    while (n != s){
        ListDigraph::Arc arc = path.predArc(n);
        n = path.predNode(n);
        pathLength += getArcLength(arc, d);
    }
    return pathLength;
}

/* Returns the actual cost of the path according to the metric used. */
double lagFlow::getPathCost(int d, DijkstraCost &path, const ListDigraph::Node &s, const ListDigraph::Node &t){
    double pathCost = 0.0;
    ListDigraph::Node n = t;
    while (n != s){
        ListDigraph::Arc arc = path.predArc(n);
        n = path.predNode(n);
        pathCost += (*coeff[d])[arc];
    }
    return pathCost;
}


/* ***********************************************************************************************************************
*                                                      RUNING METHODS
*********************************************************************************************************************** */

/* To solve an iteration of the Lagrangian - one sub problem */
void lagFlow::run(bool adaptedSubproblem){
    setCurrentLagrCost(0.0);
    setCurrentRealCost(0.0);
    
    if(getStatus() == STATUS_FEASIBLE){
        setStatus(STATUS_UNKNOWN);
    }

    resetSlacks();
    if(adaptedSubproblem){
        runAdaptedGeneralObj();
    }else{
        if(instance.getInput().isObj8(0)){
            runObj8();
        }else{
            runGeneralObj();
        }
    }
}

void lagFlow::runGeneralObj(){
    operatorCost oper(lagrangianMultiplierOverlap); 
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){

        /* Time to compute costs */
        time.setStart(ClockTime::getTimeNow());  

        double scale = getLengthMultiplier_k(d)/getToBeRouted_k(d).getMaxLength();
        ScaleMapCost scaleMap((*vecArcLength[d]),scale);  // Multiply length Map Arc by length constraint multiplier
        AddMapCost addMap((*coeff[d]),scaleMap);  // Add to the coefficients
        oper.setDemandLoad(getToBeRouted_k(d).getLoad());
        CombineMapCost combine((*vecArcLabel[d]),(*vecArcSlice[d]),oper);  // Add the overlap part
        AddMapFinalCost addMapFinal(combine,addMap);

        incCostTime(time.getTimeInSecFromStart());

        time.setStart(ClockTime::getTimeNow());
        
        const ListDigraph::Node SOURCE = getFirstNodeFromLabel(d, getToBeRouted_k(d).getSource());
        const ListDigraph::Node TARGET = getFirstNodeFromLabel(d, getToBeRouted_k(d).getTarget());

        DijkstraCost shortestPath((*vecGraph[d]), addMapFinal);
        shortestPath.run(SOURCE, TARGET);

        //BellmanFord<ListDigraph,ListDigraph::ArcMap<double>>  shortestPath((*vecGraph[d]), (*cost[d]));
        //shortestPath.run(SOURCE);
        //displayPath(shortestPath, s, t);

        if (shortestPath.reached(TARGET) == false){
            setStatus(STATUS_INFEASIBLE);
            std::cout << "> RSA is infeasible because there is no path from " << getToBeRouted_k(d).getSource()+1 << " to " << getToBeRouted_k(d).getTarget()+1 << " required for routing demand " << getToBeRouted_k(d).getId()+1 << "." << std::endl;
            return;
        }
        incShorstestPathTime(time.getTimeInSecFromStart());

        time.setStart(ClockTime::getTimeNow());
        updateAssignment_k(d, shortestPath, SOURCE, TARGET);
        incCurrentLagrCost(shortestPath.dist(TARGET));
        incCurrentRealCost(getRealCostFromPath(d, shortestPath, SOURCE, TARGET));
        incUpdateVariablesTime(time.getTimeInSecFromStart());
    }
    time.setStart(ClockTime::getTimeNow());
    subtractConstantValuesFromLagrCost();
    incSubstractMultipliersTime(time.getTimeInSecFromStart());
}

void lagFlow::runAdaptedGeneralObj(){
    operatorCost oper(lagrangianMultiplierOverlap); 
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){

        /* Time to compute costs */
        time.setStart(ClockTime::getTimeNow());  

        double scale = getLengthMultiplier_k(d)/getToBeRouted_k(d).getMaxLength();
        ScaleMapCost scaleMap((*vecArcLength[d]),scale);  // Multiply length Map Arc by length constraint multiplier
        AddMapCost addMap((*coeff[d]),scaleMap);  // Add to the coefficients
        oper.setDemandLoad(getToBeRouted_k(d).getLoad());
        CombineMapCost combine((*vecArcLabel[d]),(*vecArcSlice[d]),oper);  // Add the overlap part
        AddMapFinalCost addMapFinal(combine,addMap);

        const ListDigraph::Node SOURCE = getFirstNodeFromLabel(d, getToBeRouted_k(d).getSource());
        const ListDigraph::Node TARGET = getFirstNodeFromLabel(d, getToBeRouted_k(d).getTarget());

        CostScaling<ListDigraph,int,double> costScale((*vecGraph[d]));

        if(instance.getInput().isObj8(0)){
            operatorCostObj8 operCostObj8(getFirstNodeFromLabel(d, getToBeRouted_k(d).getSource()),getMaxUsedSliceOverallMultiplier_k(d),0.0);
            //operatorCostObj8 operCostObj8(getFirstNodeFromLabel(d, getToBeRouted_k(d).getSource()),getMaxUsedSliceOverallMultiplier_k(d),getMaxUsedSliceOverallAuxMultiplier_k(d));
            SourceMap<ListDigraph> sourceMap((*vecGraph[d])); 
            CombineMapCostObj8 combineObj8((*vecArcSlice[d]),sourceMap,operCostObj8); 
            AddMapFinalCostObj8 addMapFinalObj8(addMapFinal,combineObj8);
            costScale.costMap(addMapFinalObj8);
        }else{
            std::cout<<"yoi" <<std::endl;
            costScale.costMap(addMapFinal);
        }

        incCostTime(time.getTimeInSecFromStart());

        time.setStart(ClockTime::getTimeNow());

        costScale.upperMap(*upperBound[d]);
        costScale.lowerMap(*lowerBound[d]);
        costScale.stSupply(SOURCE,TARGET,1);
        std::cout<<"oi" <<std::endl;
        CostScaling<ListDigraph,int,double>::ProblemType problemType = costScale.run();
        std::cout<< d <<std::endl;

        if(problemType == CostScaling<ListDigraph,int,double>::INFEASIBLE){
            setStatus(STATUS_INFEASIBLE);
            std::cout << "> RSA is infeasible because there is no path from " << getToBeRouted_k(d).getSource()+1 << " to " << getToBeRouted_k(d).getTarget()+1 << " required for routing demand " << getToBeRouted_k(d).getId()+1 << "." << std::endl;
            return;
        }
        if(problemType == CostScaling<ListDigraph,int,double>::UNBOUNDED){
            std::cout << "> The problem should not be unbounded because all the variables considered have upper bound at most equal to 1." << std::endl;
        }
        incShorstestPathTime(time.getTimeInSecFromStart());

        time.setStart(ClockTime::getTimeNow());
        updateAssignment_k(d, costScale, SOURCE, TARGET);
        incCurrentLagrCost(costScale.totalCost());
        if(!instance.getInput().isObj8(0)){
            incCurrentRealCost(getRealCostFromPath(d, costScale, SOURCE, TARGET));
        }
        incUpdateVariablesTime(time.getTimeInSecFromStart());
    }
    time.setStart(ClockTime::getTimeNow());
    subtractConstantValuesFromLagrCost();
    incSubstractMultipliersTime(time.getTimeInSecFromStart());

    if(instance.getInput().isObj8(0)){
        solveProblemMaxUsedSliceOverall();
    }
}

void lagFlow::runObj8(){
    operatorCost oper(lagrangianMultiplierOverlap); 
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){

        time.setStart(ClockTime::getTimeNow());  

        double scale = getLengthMultiplier_k(d)/getToBeRouted_k(d).getMaxLength();
        ScaleMapCost scaleMap((*vecArcLength[d]),scale);
        AddMapCost addMap((*coeff[d]),scaleMap);
        oper.setDemandLoad(getToBeRouted_k(d).getLoad());
        CombineMapCost combine((*vecArcLabel[d]),(*vecArcSlice[d]),oper);
        AddMapFinalCost addMapFinal(combine,addMap);

        operatorCostObj8 operCostObj8(getFirstNodeFromLabel(d, getToBeRouted_k(d).getSource()),getMaxUsedSliceOverallMultiplier_k(d),0.0);
        //operatorCostObj8 operCostObj8(getFirstNodeFromLabel(d, getToBeRouted_k(d).getSource()),getMaxUsedSliceOverallMultiplier_k(d),getMaxUsedSliceOverallAuxMultiplier_k(d));
        SourceMap<ListDigraph> sourceMap((*vecGraph[d])); 
        CombineMapCostObj8 combineObj8((*vecArcSlice[d]),sourceMap,operCostObj8); 
        AddMapFinalCostObj8 addMapFinalObj8(addMapFinal,combineObj8);

        incCostTime(time.getTimeInSecFromStart());

        time.setStart(ClockTime::getTimeNow());

        const ListDigraph::Node SOURCE = getFirstNodeFromLabel(d, getToBeRouted_k(d).getSource());
        const ListDigraph::Node TARGET = getFirstNodeFromLabel(d, getToBeRouted_k(d).getTarget());

        DijkstraCostObj8 shortestPath((*vecGraph[d]), addMapFinalObj8);   
        shortestPath.run(SOURCE, TARGET);

        if (shortestPath.reached(TARGET) == false){
            setStatus(STATUS_INFEASIBLE);
            std::cout << "> RSA is infeasible because there is no path from " << getToBeRouted_k(d).getSource()+1 << " to " << getToBeRouted_k(d).getTarget()+1 << " required for routing demand " << getToBeRouted_k(d).getId()+1 << "." << std::endl;
            return;
        }
        incShorstestPathTime(time.getTimeInSecFromStart());

        time.setStart(ClockTime::getTimeNow());
        updateAssignment_k(d, shortestPath, SOURCE, TARGET);
        incCurrentLagrCost(shortestPath.dist(TARGET));
        incUpdateVariablesTime(time.getTimeInSecFromStart());
    }

    time.setStart(ClockTime::getTimeNow());
    subtractConstantValuesFromLagrCost();
    incSubstractMultipliersTime(time.getTimeInSecFromStart());

    solveProblemMaxUsedSliceOverall();
}

void lagFlow::subtractConstantValuesFromLagrCost(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double val = - getLengthMultiplier_k(d);
        incCurrentLagrCost(val);  
    }
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            double val = - getOverlapMultiplier_k(e, s);
            incCurrentLagrCost(val);
        }
    }
}

void lagFlow::solveProblemMaxUsedSliceOverall(){
    double exp = 0.0;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        exp += getMaxUsedSliceOverallMultiplier_k(d);
    }
    if(exp > 1.000){
        maxUsedSliceOverall = (getNbSlicesGlobalLimit()-1);
    }else{
        maxUsedSliceOverall = 0.0;
    }
    updateMaxUsedSliceOverallSlack_aux();

    incCurrentLagrCost(maxUsedSliceOverall*(1.0-exp));
    incCurrentRealCost(maxUsedSliceOverall);
}

/*void lagFlow::solveProblemMaxUsedSliceOverall(){
    std::fill(varAuxZ.begin(), varAuxZ.end(), false);
    double min = __DBL_MAX__; int ind = 0;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        if(getMaxUsedSliceOverallAuxMultiplier_k(d) < min){
            min = getMaxUsedSliceOverallAuxMultiplier_k(d);
            ind = d;
        }
    }
    varAuxZ[ind] = true; 
    
    double exp = 0.0;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        exp += getMaxUsedSliceOverallMultiplier_k(d);
        exp -= getMaxUsedSliceOverallAuxMultiplier_k(d);
    }
    if(exp > 1.000){
        maxUsedSliceOverall = (getInstance().getMaxSlice()-1);
    }else{
        maxUsedSliceOverall = 0.0;
    }

    updateMaxUsedSliceOverallSlack_aux();
    updateMaxUsedSliceOverallAuxSlack_aux();

    double sum = 0.0;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        sum += getMaxUsedSliceOverallAuxMultiplier_k(d)*(getInstance().getMaxSlice()-1)*(1-varAuxZ[d]);
    }
    incCurrentLagrCost(maxUsedSliceOverall*(1-exp)-sum);
    incCurrentRealCost(maxUsedSliceOverall);
}*/

/* Assigns length as the main cost of the arcs. */
void lagFlow::setLengthCost(){
    for(int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            (*coeff[d])[a] = getArcLength(a, d) ;
        }
    }
}

/* Stores the path found in the arcMap onPath. */
void lagFlow::updateOnPath(){
    for(int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            int index = getArcIndex(a, d);
            if (assignmentMatrix_d[d][index] == true){
                (*vecOnPath[d])[a] = getToBeRouted_k(d).getId();
            }
            else{
                (*vecOnPath[d])[a] = -1;
            }
        }
    }
}

/* ***********************************************************************************************************************
*                                                         DISPLAYS
*********************************************************************************************************************** */

std::string lagFlow::getPathString(int d, DijkstraCost &path, const ListDigraph::Node &s, const ListDigraph::Node &t){
    std::string pathString = "";
    ListDigraph::Node n = t;
    while (n != s){
        pathString += "(" + std::to_string(getNodeLabel(n, d) + 1) + "," + std::to_string(getNodeSlice(n, d) + 1) + ")";
        pathString += "-";
        n = path.predNode(n);
    }
    pathString += "(" + std::to_string(getNodeLabel(s, d) + 1) + "," + std::to_string(getNodeSlice(s, d) + 1) + ")";
    return pathString;
}

void lagFlow::displayMultiplier(std::ostream & saida){
    std::string display = "Length Multiplier = [ \n";
    for (unsigned int i = 0; i < lagrangianMultiplierLength.size(); i++){
        display += "Demand " + std::to_string(i+1) + ": " + std::to_string(getLengthMultiplier_k(i)) + "\n"; 
    }
    display += "]";
    //std::cout << display << std::endl;
    saida << display << std::endl;

    display = "Overlap Multiplier = [ \n";
    for (int e = 0; e < instance.getNbEdges(); e++){
        display += "Edge " + std::to_string(e+1) + ":\n ";
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            display += "\t Slice " + std::to_string(s+1) + ": "+ std::to_string(getOverlapMultiplier_k( e, s)) + "\n"; 
        }
    }
    display += "]";
    //std::cout << display << std::endl;
    saida << display << std::endl;

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);

    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        display = "Max Used Slice Overall Multiplier = [ \n";
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            display += "Demand " + std::to_string(d+1) + ": " + std::to_string(getMaxUsedSliceOverallMultiplier_k(d)) + "\n"; 
        }
        display += "]";
        saida << display << std::endl;

        /*display = "Max Used Slice Overall 2 Multiplier = [ \n";
        for (int e = 0; e < instance.getNbEdges(); e++){
            maxUsedSliceOverallSlack2[e].resize(getNbSlicesLimitFromEdge(e));
            display += "Edge " + std::to_string(e+1) + ":\n ";
            for (int s = 0; s < getNbSlicesLimitFromEdge(e); s++){
                display += "\t Slice " + std::to_string(s+1) + ": "+ std::to_string(getMaxUsedSliceOverall2Multiplier_k(e, s)) + "\n"; 
            }
        }
        display += "]";
        saida << display << std::endl;

        display = "Max Used Slice Overall 3 Multiplier = [ \n";
        for (int v = 0; v < instance.getNbNodes(); v++){
            maxUsedSliceOverallSlack3[v].resize(getNbSlicesGlobalLimit());
            display += "Node " + std::to_string(v+1) + ":\n ";
            for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
                display += "\t Slice " + std::to_string(s+1) + ": "+ std::to_string(getMaxUsedSliceOverall3Multiplier_k(v, s)) + "\n"; 
            }
        }
        display += "]";
        saida << display << std::endl;
        */
    }
}

void lagFlow::displaySlack(std::ostream & saida){
    std::string display = "Length Slack = [ \n";
    //std::cout << "Length:" << std::endl;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        display += "Demand " + std::to_string(d+1) + ": " +std::to_string(-getLengthSlack_k(d)) + "\n"; 
        //std::cout << "Demand:" << d << " "<< getLengthSlack_k(d) << std::endl;
    }
    display += "]";
    //std::cout << display << std::endl;
    saida << display << std::endl;

    display = "Overlap Slack = [ \n";
    for (int e = 0; e < instance.getNbEdges(); e++){
        display += "Edge " + std::to_string(e+1) + "\n ";
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            display += "\t Slice " + std::to_string(s+1) + ": "+std::to_string(-getOverlapSlack_k( e, s)) + "\n"; 
        }
    }
    display += "]";
    //std::cout << display << std::endl;
    saida << display << std::endl;

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);

    if(chosenMetric == Input::OBJECTIVE_METRIC_8){

        display = "Max Used Slice Overall Slack = [ \n";
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            display += "Demand " + std::to_string(d+1) + ": " + std::to_string(getMaxUsedSliceOverallSlack_k(d)) + "\n"; 
        }
        display += "]";
        saida << display << std::endl;

        /*display = "Max Used Slice Overall 2 Slack = [ \n";
        for (int e = 0; e < instance.getNbEdges(); e++){
            maxUsedSliceOverallSlack2[e].resize(getNbSlicesLimitFromEdge(e));
            display += "Edge " + std::to_string(e+1) + ":\n ";
            for (int s = 0; s < getNbSlicesLimitFromEdge(e); s++){
                display += "\t Slice " + std::to_string(s+1) + ": "+ std::to_string(getMaxUsedSliceOverall2Slack_k(e, s)) + "\n"; 
            }
        }
        display += "]";
        saida << display << std::endl;

        display = "Max Used Slice Overall 3 Slack = [ \n";
        for (int v = 0; v < instance.getNbNodes(); v++){
            maxUsedSliceOverallSlack3[v].resize(getNbSlicesGlobalLimit());
            display += "Node " + std::to_string(v+1) + ":\n ";
            for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
                display += "\t Slice " + std::to_string(s+1) + ": "+ std::to_string(getMaxUsedSliceOverall3Slack_k(v, s)) + "\n"; 
            }
        }
        display += "]";
        saida << display << std::endl;
        */
    }
}