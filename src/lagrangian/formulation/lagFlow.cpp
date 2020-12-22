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
}

/********************************* MULTIPLIERS ***********************************/

/* Sets the initial lagrangian multipliers for the subgradient to run. */
void lagFlow::initMultipliers(){
    initializeLengthMultipliers();
    initializeOverlapMultipliers();
    
    std::cout << "> Initial Lagrangian multipliers were defined. " << std::endl;
}

/* Sets the initial lagrangian multipliers associated with length constraints. */
void lagFlow::initializeLengthMultipliers(){
    double initialMultiplier = instance.getInput().getInitialLagrangianMultiplier();
    lagrangianMultiplierLength.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        lagrangianMultiplierLength[d] = initialMultiplier;
    }
}

/* Sets the initial lagrangian multipliers associated with non-overlapping constraints. */
void lagFlow::initializeOverlapMultipliers(){
    double initialMultiplier = instance.getInput().getInitialLagrangianMultiplier();
    lagrangianMultiplierOverlap.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        lagrangianMultiplierOverlap[e].resize(instance.getPhysicalLinkFromIndex(e).getNbSlices());
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            lagrangianMultiplierOverlap[e][s] = initialMultiplier;
        }
    }
}

/********************************** STABILITY CENTER ***************************************/

void lagFlow::initStabilityCenter(){
    initializeLengthSC();
    initializeOverlapSC();
}

/* Sets the initial lagrangian stability center associated with length constraints. */
void lagFlow::initializeLengthSC(){
    double initialSC = instance.getInput().getInitialLagrangianMultiplier();
    lagrangianSCLength.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        lagrangianSCLength[d] = initialSC;
    }
}

/* Sets the initial lagrangian stability center associated with non-overlapping constraints. */
void lagFlow::initializeOverlapSC(){
    double initialSC = instance.getInput().getInitialLagrangianMultiplier();
    lagrangianSCOverlap.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        lagrangianSCOverlap[e].resize(instance.getPhysicalLinkFromIndex(e).getNbSlices());
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            lagrangianSCOverlap[e][s] = initialSC;
        }
    }
}

/**************************** Option 2: Warmstart ****************************/



/********************************** ASSIGNMENT MATRIX ***************************************/

/* Initializes the assignement matrix. */
void lagFlow::initAssignmentMatrix(){
    assignmentMatrix_d.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        assignmentMatrix_d[d].resize(countArcs(*vecGraph[d]));
        std::fill(assignmentMatrix_d[d].begin(), assignmentMatrix_d[d].end(), false);
    }
    
    std::cout << "> Initial Assignment matrix was defined. " << std::endl;
}

/********************************** SLACK ***************************************/

/* Initializes the slack of relaxed constraints. */
void lagFlow::initSlacks(){
    initializeLengthSlacks();
    initializeOverlapSlacks();
    
    std::cout << "> Initial Slacks were defined. " << std::endl;
}

/* Initializes the slack of length constraints. */
void lagFlow::initializeLengthSlacks(){
    lengthSlack.resize(getNbDemandsToBeRouted(), 0.0);
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        //lengthSlack[d] = getToBeRouted_k(d).getMaxLength();
        lengthSlack[d] = 1;
    }
}

/* Initializes the slack of non-overlap constraints. */
void lagFlow::initializeOverlapSlacks(){
    overlapSlack.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        overlapSlack[e].resize(instance.getPhysicalLinkFromIndex(e).getNbSlices());
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            overlapSlack[e][s] = 1;
        }
    }
}


void lagFlow::initSlacks_v2(){
    initializeLengthSlacks_v2();
    initializeOverlapSlacks_v2();
}

void lagFlow::initializeLengthSlacks_v2(){
    lengthSlack_v2.resize(getNbDemandsToBeRouted(), 0.0);
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        //lengthSlack_v2[d] = getToBeRouted_k(d).getMaxLength();
        lengthSlack_v2[d] = 1;
    }
}
        
void lagFlow::initializeOverlapSlacks_v2(){
    overlapSlack_v2.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        overlapSlack_v2[e].resize(instance.getPhysicalLinkFromIndex(e).getNbSlices());
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            overlapSlack_v2[e][s] = 1;
        }
    }
}

/******************************** DIRECTION *********************************/

void lagFlow::initDirection(){
    lengthDirection.resize(getNbDemandsToBeRouted(), 0.0);
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        lengthDirection[d] = lengthSlack[d];
    }

    overlapDirection.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        overlapDirection[e].resize(instance.getPhysicalLinkFromIndex(e).getNbSlices());
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            overlapDirection[e][s] = overlapSlack[e][s];
        }
    }


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

/********************************** COSTS ***************************************/

/* Updates the arc costs according to the last lagrangian multiplier available. cost = c + u_k*length */
void lagFlow::updateCosts(){
    // Original costs
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            setArcCost(a, d, getCoeff(a, d));
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

}

/********************************** SLACK ***************************************/

void lagFlow::updateSlack(){
    updateLengthSlack();
    updateOverlapSlack();
}

void lagFlow::updateLengthSlack(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double exp = 0.0;

        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            int index = getArcIndex(a,d);
            //exp += assignmentMatrix_d[d][index]*getArcLength(a,d);
            exp += (assignmentMatrix_d[d][index]*getArcLength(a,d))/(getToBeRouted_k(d).getMaxLength());
        }
        //const int MAX_LENGTH = getToBeRouted_k(d).getMaxLength(); 
        const int MAX_LENGTH = 1;
        lengthSlack[d] = (MAX_LENGTH - exp);
    }
}

/* Updates the slack of non-overlap constraints. */
void lagFlow::updateOverlapSlack(){

    for (int e = 0; e < instance.getNbEdges(); e++){
        int linkLabel = instance.getPhysicalLinkFromIndex(e).getId();
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            double expr = 0.0;
             for (int d = 0; d < getNbDemandsToBeRouted(); d++){
                 int demandLoad = getToBeRouted_k(d).getLoad();
                for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
                    if( (getArcLabel(a, d) == linkLabel) && (getArcSlice(a, d) >= s)  && (getArcSlice(a, d) <= s + demandLoad - 1) ){
                        int id = getArcIndex(a, d);
                        expr += assignmentMatrix_d[d][id];
                    }
                }
            }
            overlapSlack[e][s] = 1 - expr;
        }
    }
}

/******** SLACK CONSIDERING THE PRIMAL VECTOR ************/

void lagFlow::updateSlack_v2(){
    updateLengthSlack_v2();
    updateOverlapSlack_v2();
}

void lagFlow::updateLengthSlack_v2(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double exp = 0.0;

        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            int index = getArcIndex(a,d);
            //exp += primal_linear_solution[d][index]*getArcLength(a,d);
            exp += (primal_linear_solution[d][index]*getArcLength(a,d))/(getToBeRouted_k(d).getMaxLength());
        }
        //const int MAX_LENGTH = getToBeRouted_k(d).getMaxLength(); 
        const int MAX_LENGTH = 1; 
        lengthSlack_v2[d] = (MAX_LENGTH - exp);
    }

}
        
void lagFlow::updateOverlapSlack_v2(){
    for (int e = 0; e < instance.getNbEdges(); e++){
        int linkLabel = instance.getPhysicalLinkFromIndex(e).getId();
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            double expr = 0.0;
             for (int d = 0; d < getNbDemandsToBeRouted(); d++){
                 int demandLoad = getToBeRouted_k(d).getLoad();
                for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
                    if( (getArcLabel(a, d) == linkLabel) && (getArcSlice(a, d) >= s)  && (getArcSlice(a, d) <= s + demandLoad - 1) ){
                        int id = getArcIndex(a, d);
                        expr += primal_linear_solution[d][id];
                    }
                }
            }
            overlapSlack_v2[e][s] = 1 - expr;
        }
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
}

void lagFlow::updateLengthDirection(double theta){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double value = getLengthSlack_k(d) + theta*getLengthDirection_k(d);
        setLengthDirection_k(d,value);
    }
}

void lagFlow::updateOverlapDirection(double theta){

    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            double value = getOverlapSlack_k(e,s) + theta*getOverlapDirection_k(e,s);
            setOverlapDirection_k(e,s,value);
        }
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

/********************************** MULTIPLIERS ***************************************/

/* Updates lagrangian multiplier with the rule: u[k+1] = u[k] + t[k]*violation */
void lagFlow::updateMultiplier(double step){
    //length
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double violation = -getLengthSlack_k(d);
        double new_multipliplier = getLengthMultiplier_k(d) + (step*violation);
        setLengthMultiplier_k(d, std::max(new_multipliplier, 0.0));
    }
  
    //overlap 
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            double violation = -getOverlapSlack_k( e, s);
            double new_multipliplier = getOverlapMultiplier_k( e, s) + (step*violation);
            setOverlapMultiplier_k( e, s, std::max(new_multipliplier, 0.0));
        }
    }  
    //displayMultiplier();
}

/******** MULTIPLIER CONSIDERING THE STABILITY CENTER ************/

void lagFlow::updateMultiplier_v2(double step){
    //length
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double violation = -getLengthSlack_v2_k(d);
        double new_multipliplier = getLengthSC_k(d) + (step*violation);
        setLengthMultiplier_k(d, std::max(new_multipliplier, 0.0));
    }
  
    //overlap 
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            double violation = -getOverlapSlack_v2_k( e, s);
            double new_multipliplier = getOverlapSC_k( e, s) + (step*violation);
            setOverlapMultiplier_k( e, s, std::max(new_multipliplier, 0.0));
        }
    }  
    //displayMultiplier();

}

/********************************** STABILITY CENTER ***************************************/

void lagFlow::updateStabilityCenter(){
    //length
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        setLengthSC_k(d, getLengthMultiplier_k(d));
    }
  
    //overlap 
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            setOverlapSC_k( e, s, getOverlapMultiplier_k( e, s));
        }
    }  

}



/********************************** CHECK FEASIBILITY ***************************************/


/* Checks if all length slacks are non-negative. */
bool lagFlow::checkLengthFeasibility(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        if (lengthSlack[d] < -DBL_EPSILON){
            return false;
        }
    }
    return true;
}

/* Checks if all overlap1 slacks are non-negative. */
bool lagFlow::checkOverlapFeasibility(){   
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            if (overlapSlack[e][s] < -DBL_EPSILON){
                return false;
            }
        }
    }
    return true;
}

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

/* Updates the step size with the rule: lambda*(UB - Z[u])/|slack| */
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
    return denominator;
}

/* Updates the step size with the rule: lambda*(UB - Z[u])/|slack| */
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
    return denominator;
}

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
    return denominator;

}

double lagFlow::getSlackDirectionProd(){

    Input::ProjectionType projection = getInstance().getInput().getChosenProjection();
    if((projection == Input::IMPROVED) || (projection ==Input::PROJECTED)){
        return getSlackDirectionProdProjected(projection);
    }   
    return getSlackDirectionProdNormal(); 
}

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
    return denominator;
}

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
    }
    return denominator;
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

double lagFlow::getMeanSlackModule_v2(){
    double module = std::sqrt(getSlackModule_v2());
    double numRest = 0.0;
    for (int e = 0; e < instance.getNbEdges(); e++){
        numRest += instance.getPhysicalLinkFromIndex(e).getNbSlices();
    }
    
    numRest += getNbDemandsToBeRouted();
    return module/numRest;
}

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
    return prod;
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
        incCurrentLagrCost(shortestPath.dist(TARGET));
        incCurrentRealCost(getRealCostFromPath(d, shortestPath, SOURCE, TARGET));
    }
    updateSlack();
    updateSlack_v2();
    updateDirection();

    int soma = 0;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            int index = getArcIndex(a,d);
            soma += assignmentMatrix_d[d][index];
        }
    }
    //std::cout << soma << std::endl;
    subtractConstantValuesFromLagrCost();

    if (checkFeasibility() == true){
        setStatus(STATUS_FEASIBLE);
    }
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
}

/* *******************************************************************************
*                             DESTRUCTOR
******************************************************************************* */

lagFlow::~lagFlow(){
    lagrangianMultiplierLength.clear();
    lagrangianMultiplierOverlap.clear();

    lengthSlack.clear();
    overlapSlack.clear();

    lengthSlack_v2.clear();
    overlapSlack_v2.clear();

    lagrangianSCLength.clear();
    lagrangianSCOverlap.clear();

    lengthDirection.clear();
    overlapDirection.clear();

    assignmentMatrix_d.clear();
    cost.clear();
}
