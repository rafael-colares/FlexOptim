#include "subgradient.h"


Subgradient::Subgradient(const Instance &inst) : RSA(inst), 
        MAX_NB_IT_WITHOUT_IMPROVEMENT(instance.getInput().getNbIterationsWithoutImprovement()), 
        MAX_NB_IT(instance.getInput().getMaxNbIterations()),
        MIN_STEPSIZE(0.00001){
    
    std::cout << "--- Subgradient was invoked ---" << std::endl;
    displayToBeRouted();

    initialization();
    std::cout << "> Subgradient was initialized. " << std::endl;

    bool STOP = false;
    while (!STOP){
        runIteration();
        if (getIsUnfeasible() == false){
            displayMainParameters();

            updateLambda();
            updateStepSize();
            updateMultiplier();
            updateCosts();
            
            if (getLB() >= getUB() - DBL_EPSILON){
                setIsOptimal(true);
                STOP = true;
            }
            if (getIteration() >= MAX_NB_IT){
                STOP = true;
            }
            if (getItWithoutImprovement() >= MAX_NB_IT_WITHOUT_IMPROVEMENT){
                STOP = true;
            }
            if (isGradientMoving() == false){
                STOP = true;
            }
        }
        else{
            STOP = true;
        }
    }
}

bool Subgradient::isGradientMoving(){
    if (getStepSize() >= MIN_STEPSIZE){
        return true;
    }
    return false;
}

/* Sets the initial lagrangian multipliers for the subgradient to run. */
void Subgradient::initMultipliers(){
    double initialMultiplier = instance.getInput().getInitialLagrangianMultiplier();
    
    lagrangianMultiplierLength.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        lagrangianMultiplierLength[d] = initialMultiplier;
    }

    lagrangianMultiplierOverlap.resize(getNbDemandsToBeRouted());
    for (int d1 = 0; d1 < getNbDemandsToBeRouted(); d1++){
        lagrangianMultiplierOverlap[d1].resize(getNbDemandsToBeRouted());
        for (int d2 = 0; d2 < getNbDemandsToBeRouted(); d2++){
            lagrangianMultiplierOverlap[d1][d2].resize(instance.getNbEdges());
            for (int e = 0; e < instance.getNbEdges(); e++){
                lagrangianMultiplierOverlap[d1][d2][e].resize(instance.getPhysicalLinkFromIndex(e).getNbSlices());
                for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
                    lagrangianMultiplierOverlap[d1][d2][e][s] = initialMultiplier;
                }
            }
        }
    }
}

/* Initializes the assignement matrix. */
void Subgradient::initAssignmentMatrix(){
    assignmentMatrix.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        assignmentMatrix[d].resize(countArcs(*vecGraph[d]));
        std::fill(assignmentMatrix[d].begin(), assignmentMatrix[d].end(), false);
    }
}

/* Initializes the slack of relaxed constraints. */
void Subgradient::initSlacks(){
    lengthSlack.resize(getNbDemandsToBeRouted(), 0.0);
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        lengthSlack[d] = getToBeRouted_k(d).getMaxLength();
    }

    overlapSlack.resize(getNbDemandsToBeRouted());
    for (int d1 = 0; d1 < getNbDemandsToBeRouted(); d1++){
        overlapSlack[d1].resize(getNbDemandsToBeRouted());
        for (int d2 = 0; d2 < getNbDemandsToBeRouted(); d2++){
            overlapSlack[d1][d2].resize(instance.getNbEdges());
            for (int e = 0; e < instance.getNbEdges(); e++){
                overlapSlack[d1][d2][e].resize(instance.getPhysicalLinkFromIndex(e).getNbSlices());
                for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
                    overlapSlack[d1][d2][e][s] = 1;
                }
            }
        }
    }
}

/* Initializes the costs in the objective function. */
void Subgradient::initCosts(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){ 
        cost.emplace_back(new ArcCost((*vecGraph[d]), 0.0));
    }
    updateCosts();
}

/* Sets the initial parameters for the subgradient to run. */
void Subgradient::initialization(){
    setIteration(0);
    setItWithoutImprovement(0);

    setLB(-__DBL_MAX__);
    setUB(__DBL_MAX__);

    setIsOptimal(false);
    setIsFeasible(false);
    setIsUnfeasible(false);

    initMultipliers();
    std::cout << "> Initial Lagrangian multipliers were defined. " << std::endl;

    setLambda(instance.getInput().getInitialLagrangianLambda());
    std::cout << "> Initial lambda was defined. " << std::endl;
    
    initAssignmentMatrix();
    std::cout << "> Initial Assignment matrix was defined. " << std::endl;

    initSlacks();
    std::cout << "> Initial Slacks were defined. " << std::endl;
    
    initCosts();
    std::cout << "> Initial costs were defined. " << std::endl;

}


/* Updates the arc costs according to the last lagrangian multiplier available. cost = c + u_k*length */
void Subgradient::updateCosts(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            setArcCost(a, d, getCoeff(a, d));
            double incrementValue = getLengthMultiplier_k(d)*getArcLength(a, d);
            incArcCost(a, d, incrementValue);
            for (int d2 = 0; d2 < getNbDemandsToBeRouted(); d2++){
                if (d != d2){
                    incrementValue = getOverlapMultiplier_k(d, d2, getArcLabel(a, d), getArcSlice(a, d));
                    incArcCost(a, d, incrementValue);
                }
            }
            /** @todo Add the second part of cost associated with non overlapping constraints. **/
        }
    }
}

/* Updates the assignment of a demand based on the a given path. */
void Subgradient::updateAssignment_k(int d, Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &SOURCE, const ListDigraph::Node &TARGET){
    std::fill(assignmentMatrix[d].begin(), assignmentMatrix[d].end(), false);
    ListDigraph::Node currentNode = TARGET;
    while (currentNode != SOURCE){
        ListDigraph::Arc arc = path.predArc(currentNode);
        int index = getArcIndex(arc, d);
        assignmentMatrix[d][index] = true;
        currentNode = path.predNode(currentNode);
    }
}

/* Checks if all slacks are non-negative. */
bool Subgradient::checkFeasibility(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        if (lengthSlack[d] < -DBL_EPSILON){
            return false;
        }
    }
    /*
    for (int d1 = 0; d1 < getNbDemandsToBeRouted(); d1++){
        for (int d2 = 0; d2 < getNbDemandsToBeRouted(); d2++){
            if (d1 != d2){
                for (int e = 0; e < instance.getNbEdges(); e++){
                    for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
                        if (overlapSlack[d1][d2][e][s] < -DBL_EPSILON){
                            return false;
                        }
                    }
                }
            }
        }
    }
    */
    return true;
}

double Subgradient::getRealCostFromPath(int d, Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &SOURCE, const ListDigraph::Node &TARGET){
    double total = 0.0;
    ListDigraph::Node currentNode = TARGET;
    while (currentNode != SOURCE){
        ListDigraph::Arc arc = path.predArc(currentNode);
        total += getCoeff(arc, d);
        currentNode = path.predNode(currentNode);
    }
    return total;
}

void Subgradient::runIteration(){
    incIteration();
    setCurrentLagrCost(0.0);
    setCurrentRealCost(0.0);
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        const ListDigraph::Node SOURCE = getFirstNodeFromLabel(d, getToBeRouted_k(d).getSource());
        const ListDigraph::Node TARGET = getFirstNodeFromLabel(d, getToBeRouted_k(d).getTarget());
        const int MAX_LENGTH = getToBeRouted_k(d).getMaxLength();
        // run a first shortest path not taking into account length (i.e., u=0)
        Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > shortestPath((*vecGraph[d]), (*cost[d]));
        shortestPath.run(SOURCE, TARGET);
        //displayPath(shortestPath, s, t);

        if (shortestPath.reached(TARGET) == false){
            setIsFeasible(false);
            setIsUnfeasible(true);
            std::cout << "> RSA is unfeasible because there is no path from " << getToBeRouted_k(d).getSource()+1 << " to " << getToBeRouted_k(d).getTarget()+1 << " required for routing demand " << getToBeRouted_k(d).getId()+1 << "." << std::endl;
            return;
        }

        updateAssignment_k(d, shortestPath, SOURCE, TARGET);
        incCurrentLagrCost(shortestPath.dist(TARGET) - (getLengthMultiplier_k(d)*MAX_LENGTH) ); /** @todo Subtract constant values. **/
        incCurrentRealCost(getRealCostFromPath(d, shortestPath, SOURCE, TARGET));
        const int PATH_LENGTH = getPathLength(d, shortestPath, SOURCE, TARGET);
        updateLengthSlack(d, PATH_LENGTH);
    }
    
    updateOverlapSlack();
    updateLB(getLagrCurrentCost());
    std::cout << "> Set LB: " << getLB() << std:: endl;
    if (checkFeasibility() == true){
        setIsFeasible(true);
        setIsUnfeasible(false);
        double feasibleSolutionCost = getRealCurrentCost();
        if (feasibleSolutionCost < getUB()){
            updateOnPath();
            updateUB(feasibleSolutionCost);
            std::cout << "> Set UB: " << getUB() << std:: endl;
        }
    }
}


/* Updates the known lower bound. */
void Subgradient::updateLB(double bound){
    if (bound > getLB()){
        setLB(bound);
        setItWithoutImprovement(0);
    }
    else{
        incItWithoutImprovement();
    }
}

/* Updates the known upper bound. */
void Subgradient::updateUB(double bound){
    if (bound < getUB()){
        setUB(bound);
    }
}

/* Updates lagrangian multiplier with the rule: u[k+1] = u[k] + t[k]*violation */
void Subgradient::updateMultiplier(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double violation = -getLengthSlack_k(d);
        double new_multipliplier = getLengthMultiplier_k(d) + (getStepSize()*violation);
        setLengthMultiplier_k(d, std::max(new_multipliplier, 0.0));
    }
    for (int d1 = 0; d1 < getNbDemandsToBeRouted(); d1++){
        for (int d2 = 0; d2 < getNbDemandsToBeRouted(); d2++){
            for (int e = 0; e < instance.getNbEdges(); e++){
                for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
                    double violation = -getOverlapSlack_k(d1, d2, e, s);
                    double new_multipliplier = getOverlapMultiplier_k(d1, d2, e, s) + (getStepSize()*violation);
                    setOverlapMultiplier_k(d1, d2, e, s, std::max(new_multipliplier, 0.0));
                }
            }
        }
    }
    //displayMultiplier();
}

/* Updates the step size with the rule: lambda*(UB - Z[u])/|slack| */
void Subgradient::updateStepSize(){
    double numerator = getLambda()*(getUB() - getLagrCurrentCost());
    double denominator = 0.0;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        denominator += std::pow(getLengthSlack_k(d), 2);
    }
    for (int d1 = 0; d1 < getNbDemandsToBeRouted(); d1++){
        for (int d2 = 0; d2 < getNbDemandsToBeRouted(); d2++){
            for (int e = 0; e < instance.getNbEdges(); e++){
                for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
                    denominator += std::pow(getOverlapSlack_k(d1, d2, e, s), 2);
                }
            }
        }
    }
    double new_stepSize = (numerator/denominator); 
    setStepSize(new_stepSize);
    //displayStepSize();
}

/* Updates the lambda used in the update of step size. Lambda is halved if LB has failed to increade in some fixed number of iterations */
void Subgradient::updateLambda(){
    double last_lambda = getLambda();
    if (getItWithoutImprovement() >= MAX_NB_IT_WITHOUT_IMPROVEMENT){
        setItWithoutImprovement(0);
        setLambda(last_lambda / 2.0);
    }
    //displayLambda();
}


/* Updates the slack of length constraint for a given path length */
void Subgradient::updateLengthSlack(int d, double pathLength){
    const int MAX_LENGTH = getToBeRouted_k(d).getMaxLength();
    lengthSlack[d] = (MAX_LENGTH - pathLength);
    //displaySlack();
}


/* Updates the slack of length constraint for a given path length */
void Subgradient::updateOverlapSlack(){
    for (int d1 = 0; d1 < getNbDemandsToBeRouted(); d1++){
        for (int d2 = 0; d2 < getNbDemandsToBeRouted(); d2++){
            for (int e = 0; e < instance.getNbEdges(); e++){
                int linkLabel = instance.getPhysicalLinkFromIndex(e).getId();
                for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
                    double expr = 0.0;
                    for (ListDigraph::ArcIt a(*vecGraph[d1]); a != INVALID; ++a){
                        if( (getArcLabel(a, d1) == linkLabel) && (getArcSlice(a, d1) == s) ){
                            int id = getArcIndex(a, d1);
                            expr += assignmentMatrix[d1][id];
                        }
                    }
                    for (ListDigraph::ArcIt a(*vecGraph[d2]); a != INVALID; ++a){
                        if( (getArcLabel(a, d2) == linkLabel) && (getArcSlice(a, d2) >= s - getToBeRouted_k(d1).getLoad() + 1) && (getArcSlice(a, d2) <= s + getToBeRouted_k(d2).getLoad() - 1) ){
                            int id = getArcIndex(a, d2);
                            expr += assignmentMatrix[d2][id];
                        }
                    }
                    overlapSlack[d1][d2][e][s] = 1 - expr;
                }
            }
        }
    }
}

/* Assigns length as the main cost of the arcs. */
void Subgradient::setLengthCost(){
    for(int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            (*cost[d])[a] = getArcLength(a, d) ;
        }
    }
}

/* Verifies if optimality condition has been achieved and update STOP flag. */
void Subgradient::updateStop(bool &STOP){
    if(getLB() >= getUB() - DBL_EPSILON){
        setIsOptimal(true);
        STOP = true;
    }
}

/* Tests if CSP is feasible by searching for a shortest path with arc costs based on their physical length. */
/*
bool Subgradient::testFeasibility(const ListDigraph::Node &s, const ListDigraph::Node &t){
    const int MAX_LENGTH = getToBeRouted()[0].getMaxLength();
    setLengthCost();
    Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > shortestLengthPath((*vecGraph[0]), cost);
    shortestLengthPath.run(s,t);
    //displayPath(shortestLengthPath, s, t);
    double smallestLength = getPathLength(shortestLengthPath, s, t);
    if( smallestLength >= MAX_LENGTH + DBL_EPSILON ){
        std::cout << "> CSP is unfeasiable." << std:: endl;
        return false;
    }
    else{
        updateUB(getPathCost(shortestLengthPath, s, t));
        std::cout << "> Set UB: " << getUB() << std:: endl;
    }
    return true;
}
*/
/* Stores the path found in the arcMap onPath. */
void Subgradient::updateOnPath(){
    for(int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            int index = getArcIndex(a, d);
            if (assignmentMatrix[d][index] == true){
                (*vecOnPath[d])[a] = getToBeRouted_k(d).getId();
            }
            else{
                (*vecOnPath[d])[a] = -1;
            }
        }
    }
}


/* Returns the physical length of the path. */
double Subgradient::getPathLength(int d, Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &s, const ListDigraph::Node &t){
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
double Subgradient::getPathCost(int d, Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &s, const ListDigraph::Node &t){
    double pathCost = 0.0;
    ListDigraph::Node n = t;
    while (n != s){
        ListDigraph::Arc arc = path.predArc(n);
        n = path.predNode(n);
        pathCost += getCoeff(arc, d);
    }
    return pathCost;
}
/*
void Subgradient::displayPath(Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &s, const ListDigraph::Node &t){
    std::cout << "The path found is:" << std::endl;
    ListDigraph::Node n = t;
    while (n != s){
        ListDigraph::Arc arc = path.predArc(n);
        n = path.predNode(n);
        displayArc(0, arc);
    }
    double pathLength = getPathLength(path, s, t);
    std::cout << "Path length: " << pathLength << std::endl;
    double pathCost = getPathCost(path, s, t);
    std::cout << "Path cost: " << pathCost << std::endl;
    double lagrangianObjFunc = path.dist(t);
    std::cout << "Lagrangian Objective Function: " << lagrangianObjFunc << std::endl;
}
*/
std::string Subgradient::getPathString(int d, Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &s, const ListDigraph::Node &t){
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
 

void Subgradient::displayMainParameters(){
    int k = getIteration();
    int sizeOfField[7];
    sizeOfField[0] = 6;
    sizeOfField[1] = 6;
    sizeOfField[2] = 6;
    sizeOfField[3] = 6;
    sizeOfField[4] = 6;
    sizeOfField[5] = 8;
    sizeOfField[6] = 10;
    char space = ' ';
    std::string field[8];
    if (k == 1){
        field[0] = "Iter";
        field[1] = "LB";
        field[2] = "f(x)";
        field[3] = "UB";
        field[4] = "Step";
        field[5] = "Lambda";
        field[6] = "Feasible";
        
        for (int i = 0; i < 7; i++){
            field[i].resize(sizeOfField[i], space);
            std::cout << field[i] << " | ";
        }
        std::cout << std::endl;
    }
    field[0] = std::to_string(k);
    field[1] = std::to_string(getLB());
    field[2] = std::to_string(getLagrCurrentCost());
    field[3] = std::to_string(getUB());
    field[4] = std::to_string(getStepSize());
    field[5] = std::to_string(getLambda());
    if (checkFeasibility()){
        field[6] = "YES";
    }
    else{
        field[6] = "NO";
    }
    for (int i = 0; i < 7; i++){
        field[i].resize(sizeOfField[i], space);
        std::cout << field[i] << " | ";
    }
    std::cout << std::endl;
}
/*
void Subgradient::displayMultiplier(){
    std::string display = "Multiplier = [ ";
    for (unsigned int i = 0; i < lagrangianMultiplier.size(); i++){
        display += std::to_string(getMultiplier_k(i)) + " "; 
    }
    display += "]";
    std::cout << display << std::endl;
}
void Subgradient::displayLambda(){
    std::string display = "Lambda = [ ";
    for (unsigned int i = 0; i < lambda.size(); i++){
        display += std::to_string(getLambda_k(i)) + " "; 
    }
    display += "]";
    std::cout << display << std::endl;
}
void Subgradient::displaySlack(){
    std::string display = "Slack = [ ";
    for (unsigned int i = 0; i < slack.size(); i++){
        display += std::to_string(getSlack_k(i)) + " "; 
    }
    display += "]";
    std::cout << display << std::endl;
}
void Subgradient::displayStepSize(){
    std::string display = "StepSize = [ ";
    for (unsigned int i = 0; i < stepSize.size(); i++){
        display += std::to_string(getStepSize_k(i)) + " "; 
    }
    display += "]";
    std::cout << display << std::endl;
}*/