#include "lagFlow.h"

/* *******************************************************************************
*                             INITIALIZATION METHODS
******************************************************************************* */

void lagFlow::init(){
    initMultipliers();
    initSlacks();
    initSlacks_v2();
    initDirection();
    initCosts();
    initAssignmentMatrix();
    initStabilityCenter();
    initPrimalSolution();
}

/********************************* MULTIPLIERS ***********************************/

/* Sets the initial lagrangian multipliers for the subgradient to run. */
void lagFlow::initMultipliers(){
    double initialMultiplier = instance.getInput().getInitialLagrangianMultiplier();
    initializeLengthMultipliers(initialMultiplier);
    initializeOverlapMultipliers(initialMultiplier);

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        initializeMaxUsedSliceOverallMultipliers(initialMultiplier);
        initializeMaxUsedSliceOverall2Multipliers(initialMultiplier);
        initializeMaxUsedSliceOverall3Multipliers(initialMultiplier);
    }
    
    std::cout << "> Initial Lagrangian multipliers were defined. " << std::endl;
}

/**************************** Option 2: Warmstart ****************************/

/********************************** STABILITY CENTER ***************************************/

void lagFlow::initStabilityCenter(){
    initializeLengthSC();
    initializeOverlapSC();

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        initializeMaxUsedSliceOverallSC();
        initializeMaxUsedSliceOverall2SC();
        initializeMaxUsedSliceOverall3SC();
    }

    std::cout << "> Initial Stability Centers were defined. " << std::endl;
}

/********************************** SLACK ***************************************/

/* Initializes the slack of relaxed constraints. */
void lagFlow::initSlacks(){
    initializeLengthSlacks();
    initializeOverlapSlacks();

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        initializeMaxUsedSliceOverallSlacks();
        initializeMaxUsedSliceOverall2Slacks();
        initializeMaxUsedSliceOverall3Slacks();
    }
    
    std::cout << "> Initial Slacks were defined. " << std::endl;
}

/********************************** SLACK CONSIDERING PRIMAL VARIABLES ***************************************/

void lagFlow::initSlacks_v2(){
    initializeLengthSlacks_v2();
    initializeOverlapSlacks_v2();

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        initializeMaxUsedSliceOverallSlacks_v2();
        initializeMaxUsedSliceOverall2Slacks_v2();
        initializeMaxUsedSliceOverall3Slacks_v2();
    }

    std::cout << "> Initial Slacks 2 were defined. " << std::endl;
}

/******************************** DIRECTION *********************************/

void lagFlow::initDirection(){
    initializeLengthDirection();
    initializeOverlapDirection();

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        initializeMaxUsedSliceOverallDirection();
        initializeMaxUsedSliceOverall2Direction();
        initializeMaxUsedSliceOverall3Direction();
    }

    std::cout << "> Initial Direction were defined. " << std::endl;
}

/********************************** ASSIGNMENT MATRIX ***************************************/

/* Initializes the assignement matrix. */
void lagFlow::initAssignmentMatrix(){
    AbstractLagFormulation::initAssignmentMatrix(); 
    std::cout << "> Initial Assignment matrix was defined. " << std::endl;
}

/********************************** COST ***************************************/

/* Initializes the costs in the objective function. */
void lagFlow::initCosts(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        cost.emplace_back( std::make_shared<ArcCost>((*vecGraph[d]), 0.0)); 
        //Old version: If compile is ok, delete this line.
        //cost.emplace_back(new ArcCost((*vecGraph[d]), 0.0));
    }
    updateCosts();

    std::cout << "> Initial costs were defined. " << std::endl;
}

/* *******************************************************************************
*                             UPDATE METHODS
******************************************************************************* */

/********************************** MULTIPLIERS ***************************************/

/* Updates lagrangian multiplier with the rule: u[k+1] = u[k] + t[k]*violation */
void lagFlow::updateMultiplier(double step){
    updateLengthMultiplier(step);
    updateOverlapMultiplier(step);

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        updateMaxUsedSliceOverallMultiplier(step);
        updateMaxUsedSliceOverall2Multiplier(step);
        updateMaxUsedSliceOverall3Multiplier(step);
    }
  
    //displayMultiplier();
}

/************** MULTIPLIER CONSIDERING THE STABILITY CENTER ****************/

void lagFlow::updateMultiplier_v2(double step){
    updateLengthMultiplier_v2(step);
    updateOverlapMultiplier_v2(step);
    //displayMultiplier();

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        updateMaxUsedSliceOverallMultiplier_v2(step);
        updateMaxUsedSliceOverall2Multiplier_v2(step);
        updateMaxUsedSliceOverall3Multiplier_v2(step);
    }
}

/********************************** STABILITY CENTER ***************************************/

void lagFlow::updateStabilityCenter(){
    updateLengthSC();
    updateOverlapSC();

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        updateMaxUsedSliceOverallSC();
        updateMaxUsedSliceOverall2SC();
        updateMaxUsedSliceOverall3SC();
    }
}

/********************************** COSTS ***************************************/

/* Updates the arc costs according to the last lagrangian multiplier available. cost = c + u_k*length */
void lagFlow::updateCosts(){
    // Original costs
    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric != Input::OBJECTIVE_METRIC_8){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
                setArcCost(a, d, getCoeff(a, d));
            }
        }
    }else{
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
                setArcCost(a, d, 0.0);
            }
        }
    }

    // Length multipliers
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double multiplier =  getLengthMultiplier_k(d);
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            //double incrementValue = multiplier*getArcLength(a, d);
            double incrementValue = multiplier*getArcLength(a, d)/getToBeRouted_k(d).getMaxLength();
            incArcCost(a, d, incrementValue);
        }
    }

    // Overlap multipliers 
    for (int e = 0; e < instance.getNbEdges(); e++){
        int linkLabel = instance.getPhysicalLinkFromIndex(e).getId();
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            double multiplier = getOverlapMultiplier_k(e, s);
            for (int d = 0; d < getNbDemandsToBeRouted(); d++){
                int demandLoad = getToBeRouted_k(d).getLoad();
                for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
                    if( (getArcLabel(a, d) == linkLabel) && (getArcSlice(a, d) >= s)  && (getArcSlice(a, d) <= s + demandLoad - 1) ){
                        double incrementValue = multiplier;
                        incArcCost(a, d, incrementValue);
                    }
                }
            }
        }
    }

    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            double multiplier =  getMaxUsedSliceOverallMultiplier_k(d);
            for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
                if (getToBeRouted_k(d).getSource() == getNodeLabel((*vecGraph[d]).source(a), d)){
                    //double incrementValue = multiplier*getArcSlice(a, d)/getInstance().getMaxSlice();
                    double incrementValue = multiplier*getArcSlice(a, d);
                    incArcCost(a,d,incrementValue);
                }
            }
        }
        
        for (int linkLabel = 0; linkLabel < instance.getNbEdges(); linkLabel++){
            int sliceLimit = getNbSlicesLimitFromEdge(linkLabel);
            for (int s = 0; s < sliceLimit; s++){
                double multiplier =  getMaxUsedSliceOverall2Multiplier_k(linkLabel,s);
                for (int d = 0; d < getNbDemandsToBeRouted(); d++){
                    int demandLoad = getToBeRouted_k(d).getLoad();
                    for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
                        if ((getArcLabel(a, d) == linkLabel) && (getArcSlice(a, d) >= s) && (getArcSlice(a, d) <= s + demandLoad - 1)){
                            //double incrementValue = multiplier*getArcSlice(a, d)/getInstance().getMaxSlice();
                            double incrementValue = multiplier*getArcSlice(a, d);
                            incArcCost(a,d,incrementValue);
                        }
                    }
                }
            }
        }

        for (int nodeLabel = 0; nodeLabel < instance.getNbNodes(); nodeLabel++){
            for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
                double multiplier = getMaxUsedSliceOverall3Multiplier_k(nodeLabel,s);
                for (int d = 0; d < getNbDemandsToBeRouted(); d++){
                    int demandLoad = getToBeRouted_k(d).getLoad();
                    for (ListDigraph::NodeIt v(*vecGraph[d]); v != INVALID; ++v){
                        if ((getNodeLabel(v, d) == nodeLabel)){
                            for (ListDigraph::OutArcIt a(*vecGraph[d], v); a != INVALID; ++a){
                                if ( (getArcSlice(a, d) >= s) && (getArcSlice(a, d) <= s + demandLoad - 1) ){
                                    //double incrementValue = multiplier*getArcSlice(a, d)/getInstance().getMaxSlice();
                                    double incrementValue = multiplier*getArcSlice(a, d);
                                    incArcCost(a,d,incrementValue);
                                }
                            }
                        }
                    }
                }
            }
        }

    }
    
}

/********************************** SLACK ***************************************/

void lagFlow::updateSlack(){
    updateLengthSlack();
    updateOverlapSlack();

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        updateMaxUsedSliceOverallSlack();
        updateMaxUsedSliceOverall2Slack();
        updateMaxUsedSliceOverall3Slack();
    }
}

/************************* SLACK CONSIDERING THE PRIMAL VECTOR *******************/

void lagFlow::updateSlack_v2(){
    updateLengthSlack_v2();
    updateOverlapSlack_v2();

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        updateMaxUsedSliceOverallSlack_v2();
        updateMaxUsedSliceOverall2Slack_v2();
        updateMaxUsedSliceOverall3Slack_v2();
    }
}
        
/******************************** DIRECTION *********************************/

void lagFlow::updateDirection(){

    Input::DirectionMethod chosenDirectionMethod = getInstance().getInput().getChosenDirectionMethod();
    double theta= 0.0;
    if(chosenDirectionMethod == Input::NORMAL){
        theta = 0.0;
    }else if(chosenDirectionMethod == Input::CROWDER){
        theta = getInstance().getInput().getCrowderParameter();
    }else if(chosenDirectionMethod == Input::CARMERINI){
        if(getSlackDirectionProdNormal() < 0.0){
            theta = -(getInstance().getInput().getCarmeriniParameter()*getSlackDirectionProdNormal())/getDirectionModule();
        }else{
            theta = 0.0;
        }
    }else if(chosenDirectionMethod == Input::MODIFIED_CARMERINI){
        if(getSlackDirectionProdNormal() < 0.0){
            theta = std::sqrt(getSlackModule())/std::sqrt(getDirectionModule());
        }else{
            theta = 0.0;
        }
    }
    updateLengthDirection(theta);
    updateOverlapDirection(theta);

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        updateMaxUsedSliceOverallDirection(theta);
        updateMaxUsedSliceOverall2Direction(theta);
        updateMaxUsedSliceOverall3Direction(theta);
    }
}

/********************************** ASSIGNMENT MATRIX ***************************************/

/* Updates the assignment of a demand based on the a given path. */
void lagFlow::updateAssignment_k(int d, Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &SOURCE, const ListDigraph::Node &TARGET){
    std::fill(assignmentMatrix_d[d].begin(), assignmentMatrix_d[d].end(), false);
    ListDigraph::Node currentNode = TARGET;
    while (currentNode != SOURCE){
        ListDigraph::Arc arc = path.predArc(currentNode);
        int index = getArcIndex(arc, d);
        assignmentMatrix_d[d][index] = true;
        currentNode = path.predNode(currentNode);
    }
}

/********************************** CHECK FEASIBILITY ***************************************/

/* Checks if all slacks are non-negative. */
bool lagFlow::checkFeasibility(){
    if (checkLengthFeasibility() == false){
        return false;
    }

    if (checkOverlapFeasibility() == false){
        return false;
    }

    return true;
    
}

/* *******************************************************************************
*                                GET METHODS
******************************************************************************* */

double lagFlow::getRealCostFromPath(int d, Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &SOURCE, const ListDigraph::Node &TARGET){
    double total = 0.0;
    ListDigraph::Node currentNode = TARGET;
    while (currentNode != SOURCE){
        ListDigraph::Arc arc = path.predArc(currentNode);
        total += getCoeff(arc, d);
        currentNode = path.predNode(currentNode);
    }
    return total;
}

/* Returns |slack|^2 */
double lagFlow::getSlackModule() {

    double denominator = 0.0;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        denominator += std::pow(getLengthSlack_k(d), 2);
    }
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            denominator += std::pow(getOverlapSlack_k( e, s), 2);
        }
    }

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            denominator += std::pow(getMaxUsedSliceOverallSlack_k(d), 2);
        }

        for (int e = 0; e < instance.getNbEdges(); e++){
            int sliceLimit = getNbSlicesLimitFromEdge(e);
            for (int s = 0; s < sliceLimit; s++){
                denominator += std::pow(getMaxUsedSliceOverall2Slack_k(e,s), 2);
            }
        }

        for (int v = 0; v < instance.getNbNodes(); v++){
            for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
                denominator += std::pow(getMaxUsedSliceOverall3Slack_k(v,s), 2);
            }
        }
    }
    return denominator;
}

/* Returns |slack|^2, with the slack considering the primal variables*/
double lagFlow::getSlackModule_v2() {
    double denominator = 0.0;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        denominator += std::pow(getLengthSlack_v2_k(d), 2);
    }
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            denominator += std::pow(getOverlapSlack_v2_k( e, s), 2);
        }
    }

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            denominator += std::pow(getMaxUsedSliceOverallSlack_v2_k(d), 2);
        }

        for (int e = 0; e < instance.getNbEdges(); e++){
            int sliceLimit = getNbSlicesLimitFromEdge(e);
            for (int s = 0; s < sliceLimit; s++){
                denominator += std::pow(getMaxUsedSliceOverall2Slack_v2_k(e,s), 2);
            }
        }

        for (int v = 0; v < instance.getNbNodes(); v++){
            for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
                denominator += std::pow(getMaxUsedSliceOverall3Slack_v2_k(v,s), 2);
            }
        }
    }
    return denominator;
}

/* Returns |Direction|^2 */
double lagFlow::getDirectionModule(){

    double denominator = 0.0;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        denominator += std::pow(getLengthDirection_k(d), 2);
    }
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            denominator += std::pow(getOverlapDirection_k( e, s), 2);
        }
    }

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            denominator += std::pow(getMaxUsedSliceOverallDirection_k(d), 2);
        }

        for (int e = 0; e < instance.getNbEdges(); e++){
            int sliceLimit = getNbSlicesLimitFromEdge(e);
            for (int s = 0; s < sliceLimit; s++){
                denominator += std::pow(getMaxUsedSliceOverall2Direction_k(e,s), 2);
            }
        }

        for (int v = 0; v < instance.getNbNodes(); v++){
            for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
                denominator += std::pow(getMaxUsedSliceOverall3Direction_k(v,s), 2);
            }
        }
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
        }
        
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
            }

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
            }

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
        }

    }
    return denominator;
}

/* Returns mean of |slack|, slack considering the primal variables */
double lagFlow::getMeanSlackModule_v2(){
    double module = std::sqrt(getSlackModule_v2());
    double numRest = 0.0;
    for (int e = 0; e < instance.getNbEdges(); e++){
        numRest += instance.getPhysicalLinkFromIndex(e).getNbSlices();
    }
    numRest += getNbDemandsToBeRouted();

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        numRest = getNbDemandsToBeRouted();
        for (int e = 0; e < instance.getNbEdges(); e++){
            numRest += getNbSlicesLimitFromEdge(e);
        }
        for (int e = 0; e < instance.getNbNodes(); e++){
            numRest += getNbSlicesGlobalLimit();
        }
    }

    return module/numRest;
}

/* Returns (slack*slack_v2)*/ 
double lagFlow::get_prod_slack(){
    double prod = 0.0;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        prod += getLengthSlack_k(d)*getLengthSlack_v2_k(d);
    }
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            prod += getOverlapSlack_k( e, s)*getOverlapSlack_v2_k( e, s);
        }
    }

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            prod += getMaxUsedSliceOverallSlack_k(d)*getMaxUsedSliceOverallSlack_v2_k(d);
        }

        for (int e = 0; e < instance.getNbEdges(); e++){
            int sliceLimit = getNbSlicesLimitFromEdge(e);
            for (int s = 0; s < sliceLimit; s++){
                prod += getMaxUsedSliceOverall2Slack_k(e,s)*getMaxUsedSliceOverall2Slack_v2_k(e,s);
            }
        }

        for (int v = 0; v < instance.getNbNodes(); v++){
            for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
                prod += getMaxUsedSliceOverall3Slack_k(v,s)*getMaxUsedSliceOverall2Slack_v2_k(v,s);
            }
        }
    }


    return prod;
}

/* Returns the physical length of the path. */
double lagFlow::getPathLength(int d, Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &s, const ListDigraph::Node &t){
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
double lagFlow::getPathCost(int d, Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &s, const ListDigraph::Node &t){
    double pathCost = 0.0;
    ListDigraph::Node n = t;
    while (n != s){
        ListDigraph::Arc arc = path.predArc(n);
        n = path.predNode(n);
        pathCost += getCoeff(arc, d);
    }
    return pathCost;
}

/* *******************************************************************************
*                             RUNING METHODS
******************************************************************************* */

void lagFlow::run(){

    setCurrentLagrCost(0.0);
    setCurrentRealCost(0.0);
    if(getStatus() == STATUS_FEASIBLE){
        setStatus(STATUS_UNKNOWN);
    }

    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        const ListDigraph::Node SOURCE = getFirstNodeFromLabel(d, getToBeRouted_k(d).getSource());
        const ListDigraph::Node TARGET = getFirstNodeFromLabel(d, getToBeRouted_k(d).getTarget());
        // run a first shortest path not taking into account length (i.e., u=0)
        Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > shortestPath((*vecGraph[d]), (*cost[d]));
        shortestPath.run(SOURCE, TARGET);
        //displayPath(shortestPath, s, t);

        if (shortestPath.reached(TARGET) == false){
            //setIsFeasible(false);
            //setIsUnfeasible(true);
            setStatus(STATUS_INFEASIBLE);
            std::cout << "> RSA is infeasible because there is no path from " << getToBeRouted_k(d).getSource()+1 << " to " << getToBeRouted_k(d).getTarget()+1 << " required for routing demand " << getToBeRouted_k(d).getId()+1 << "." << std::endl;
            return;
        }
        updateAssignment_k(d, shortestPath, SOURCE, TARGET);
        if(chosenMetric != Input::OBJECTIVE_METRIC_8){
            incCurrentRealCost(getRealCostFromPath(d, shortestPath, SOURCE, TARGET));
        }
        //std::cout << "DIST " << shortestPath.dist(TARGET) << std::endl;
        incCurrentLagrCost(shortestPath.dist(TARGET));
    }
    subtractConstantValuesFromLagrCost();
    //std::cout << "DIST SUB " << getLagrCurrentCost() << std::endl;
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        solveProblemMaxUsedSliceOverall();
    }

    /**
    updateSlack();
    updateSlack_v2();
    updateDirection();
 
    if (checkFeasibility() == true){
        setStatus(STATUS_FEASIBLE);
    }
    **/

    int soma = 0;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            int index = getArcIndex(a,d);
            soma += assignmentMatrix_d[d][index];
        }
    }
    
    //std::cout << "SOMA " << soma << std::endl;
}

void lagFlow::subtractConstantValuesFromLagrCost(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        //double val = - getLengthMultiplier_k(d)*getToBeRouted_k(d).getMaxLength();
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
    //std::cout << "EXP1 " << exp <<std::endl;

    for (int e = 0; e < instance.getNbEdges(); e++){
        int sliceLimit = getNbSlicesLimitFromEdge(e);
        for (int s = 0; s < sliceLimit; s++){
            exp += getMaxUsedSliceOverall2Multiplier_k(e,s);
        }
    }

    for (int v = 0; v < instance.getNbNodes(); v++){
        for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
            exp += getMaxUsedSliceOverall3Multiplier_k(v,s);
        }
    }
    

    double aux = 0.0;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            if (getToBeRouted_k(d).getSource() == getNodeLabel((*vecGraph[d]).source(a), d)){
                int indexArc = getArcIndex(a,d);
                if(getArcSlice(a, d)*assignmentMatrix_d[d][indexArc] >aux){
                    aux = getArcSlice(a, d)*assignmentMatrix_d[d][indexArc];
                }
            }
        }
    }

    //if(exp/getInstance().getMaxSlice() > 1.000){
    if(exp > 1.000){
        maxUsedSliceOverall = getInstance().getMaxSlice();
        //maxUsedSliceOverall = aux;
        //maxUsedSliceOverall = maxUsedSliceOverall/getInstance().getMaxSlice();
    }else{
        maxUsedSliceOverall = 0.0;
    }
    incCurrentLagrCost(maxUsedSliceOverall*(1-exp));
    //incCurrentLagrCost(aux*(1-exp));
    //std::cout << "DIST3 " << getLagrCurrentCost() <<std::endl;
    //incCurrentRealCost(maxUsedSliceOverall);
    incCurrentRealCost(maxUsedSliceOverall);
}

/* Assigns length as the main cost of the arcs. */
void lagFlow::setLengthCost(){
    for(int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            (*cost[d])[a] = getArcLength(a, d) ;
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

/* *******************************************************************************
*                                 DISPLAYS
******************************************************************************* */

std::string lagFlow::getPathString(int d, Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &s, const ListDigraph::Node &t){
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

        display = "Max Used Slice Overall 2 Multiplier = [ \n";
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

        display = "Max Used Slice Overall 2 Slack = [ \n";
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
    }
}

/* *******************************************************************************
*                             DESTRUCTOR
******************************************************************************* */

lagFlow::~lagFlow(){
    
    cost.clear();
}
