#include "subgradient.h"


Subgradient::Subgradient(const Instance &inst) : RSA(inst), 
        SOURCE(getToBeRouted()[0].getSource()), TARGET(getToBeRouted()[0].getTarget()), 
        MAX_LENGTH(getToBeRouted()[0].getMaxLength()), 
        MAX_NB_IT_WITHOUT_IMPROVEMENT(inst.getInput().getNbIterationsWithoutImprovement()), 
        MAX_NB_IT(inst.getInput().getMaxNbIterations()), cost(*vecGraph[0]) {
    
    
    std::cout << "--- Subgradient was invoked ---" << std::endl;
    std::cout << "> Route demand " << SOURCE+1 << "->" << TARGET+1 << std::endl;
    std::cout << "> Max Length: " << MAX_LENGTH << std::endl << std::endl;

    initialization();
    std::cout << "> Subgradient was initialized. " << std::endl;

    run(getFirstNodeFromLabel(0, SOURCE), getFirstNodeFromLabel(0, TARGET));

}

/* Sets the initial parameters for the subgradient to run. */
void Subgradient::initialization(){
    setIteration(0);
    setItWithoutImprovement(0);
    setLB(-__DBL_MAX__);
    setUB(__DBL_MAX__);
    setIsOptimal(false);

    double initialMultiplier = instance.getInput().getInitialLagrangianMultiplier();
    lagrangianMultiplier.push_back(initialMultiplier);
    std::cout << "> Initial lagrangian multiplier was defined. " << std::endl;

    subgradientPreprocessing();
    updateCosts();

    double initialLambda = instance.getInput().getInitialLagrangianLambda();
    lambda.push_back(initialLambda);
    std::cout << "> Initial lambda was defined. " << std::endl;

}

/* Call preprocessing methods. */
void Subgradient::subgradientPreprocessing(){
    contractNodesFromLabel(0, SOURCE);
    contractNodesFromLabel(0, TARGET);
}

/* Updates the arc costs according to the last lagrangian multiplier available. cost = c + u_k*length */
void Subgradient::updateCosts(){
    for (ListDigraph::ArcIt a(*vecGraph[0]); a != INVALID; ++a){
        double arcCost = getCoeff(a, 0);
        cost[a] = arcCost + ( getLastMultiplier()*getArcLength(a, 0) );
    }
}

/* Solves the Constrained Shortest Path from node s to node t using the Subgradient Method. */
void Subgradient::run(const ListDigraph::Node &s, const ListDigraph::Node &t){
    bool STOP = false;
    
    // run a first shortest path not taking into account length (i.e., u=0)
    Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > shortestPath((*vecGraph[0]), cost);
    shortestPath.run(s,t);
    //displayPath(shortestPath, s, t);
    
    if (shortestPath.reached(t) == false){
        std::cout << "> CSP is not feasible. " << std::endl;
        exit(0);
    }
    
    setCurrentCost(shortestPath.dist(t) - (MAX_LENGTH*getMultiplier_k(getIteration())));
    updateLB(getCurrentCost());
    std::cout << "> Set LB: " << getLB() << std:: endl;
    updateSlack(getPathLength(shortestPath, s, t));

    // if the path found is unfeasible, check for problem's feasibility
    if( getPathLength(shortestPath, s, t) >= MAX_LENGTH + EPSILON ){
        if (testFeasibility(s,t) == false){
            STOP = true;
            setIsFeasible(false);
            std::cout << "> CSP is not feasible. " << std::endl;
            exit(0);
        }
    }
    else{
        double ub = getPathCost(shortestPath, s, t);
        updateUB(ub);
        std::cout << "> Set UB: " << getUB() << std:: endl;
        setIsFeasible(true);
        updateOnPath(shortestPath, s, t);
    }

    if (getLB() >= getUB() - EPSILON){
        setIsOptimal(true);
        STOP = true;
    }
    
    
    updateStepSize();
    displayMainParameters(shortestPath, s, t);
    while (!STOP){
        updateMultiplier();
        incIteration();

        updateCosts();
        shortestPath.lengthMap(cost);
        shortestPath.run(s,t);

        setCurrentCost(shortestPath.dist(t) - (MAX_LENGTH*getMultiplier_k(getIteration())));
        updateSlack(getPathLength(shortestPath, s, t));
        updateLB(getCurrentCost());
        double new_path_cost = getPathCost(shortestPath, s, t);
        if((getSlack_k(getIteration()) >= 0.0 - EPSILON) && (new_path_cost < getUB())){
            updateUB(new_path_cost);
            updateOnPath(shortestPath, s, t);
        }
        updateLambda();
        updateStepSize();

        displayMainParameters(shortestPath, s, t);
        updateStop(STOP);

        if (getIteration() >= MAX_NB_IT){
            STOP = true;
        }
        if(getLB() >= getUB() - EPSILON){
            setIsOptimal(true);
            STOP = true;
        }
        //displayPath(shortestPath, s, t);
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
    double violation = -getLastSlack();
    double new_multipliplier = getLastMultiplier() + (getLastStepSize()*violation);
    lagrangianMultiplier.push_back(std::max(new_multipliplier,0.0));

    //displayMultiplier();
}

/* Updates the step size with the rule: lambda*(UB - Z[u])/|slack| */
void Subgradient::updateStepSize(){
    int k = getIteration();
    double numerator = getLambda_k(k)*(getUB() - getCurrentCost());
    double denominator = std::pow(std::sqrt(getSlack_k(k)*getSlack_k(k)), 2);
    double new_stepSize = (numerator/denominator); 
    stepSize.push_back(new_stepSize);

    //displayStepSize();
}

/* Updates the lambda used in the update of step size. Lambda is halved if LB has failed to increade in some fixed number of iterations */
void Subgradient::updateLambda(){
    double new_lambda = getLastLambda();
    if (getItWithoutImprovement() >= MAX_NB_IT_WITHOUT_IMPROVEMENT){
        setItWithoutImprovement(0);
        new_lambda = new_lambda / 2;
    }
    lambda.push_back(new_lambda);

    //displayLambda();
}


/* Updates the slack of length constraint for a given path length */
void Subgradient::updateSlack(double pathLength){
    slack.push_back(MAX_LENGTH - pathLength);

    //displaySlack();
}

/* Assigns length as the main cost of the arcs. */
void Subgradient::setLengthCost(){
    for (ListDigraph::ArcIt a(*vecGraph[0]); a != INVALID; ++a){
        cost[a] = getArcLength(a, 0) ;
    }
}

/* Verifies if optimality condition has been achieved and update STOP flag. */
void Subgradient::updateStop(bool &STOP){
    if(getLB() >= getUB() - EPSILON){
        setIsOptimal(true);
        STOP = true;
    }
}

/* Tests if CSP is feasible by searching for a shortest path with arc costs based on their physical length. */
bool Subgradient::testFeasibility(const ListDigraph::Node &s, const ListDigraph::Node &t){
    setLengthCost();
    Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > shortestLengthPath((*vecGraph[0]), cost);
    shortestLengthPath.run(s,t);
    //displayPath(shortestLengthPath, s, t);
    double smallestLength = getPathLength(shortestLengthPath, s, t);
    if( smallestLength >= MAX_LENGTH + EPSILON ){
        std::cout << "> CSP is unfeasiable." << std:: endl;
        return false;
    }
    else{
        updateUB(getPathCost(shortestLengthPath, s, t));
        std::cout << "> Set UB: " << getUB() << std:: endl;
    }
    return true;
}

/* Stores the path found in the arcMap onPath. */
void Subgradient::updateOnPath(Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &s, const ListDigraph::Node &t){
    ListDigraph::Node n = t;
    for (ListDigraph::ArcIt a(*vecGraph[0]); a != INVALID; ++a){
        (*vecOnPath[0])[a] = -1;
    }
    while (n != s){
        ListDigraph::Arc arc = path.predArc(n);
        n = path.predNode(n);
        (*vecOnPath[0])[arc] = getToBeRouted()[0].getId();
    }
}


/* Returns the physical length of the path. */
double Subgradient::getPathLength(Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &s, const ListDigraph::Node &t){
    double pathLength = 0.0;
    ListDigraph::Node n = t;
    while (n != s){
        ListDigraph::Arc arc = path.predArc(n);
        n = path.predNode(n);
        pathLength += getArcLength(arc, 0);
    }
    return pathLength;
}

/* Returns the actual cost of the path according to the metric used. */
double Subgradient::getPathCost(Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &s, const ListDigraph::Node &t){
    double pathCost = 0.0;
    ListDigraph::Node n = t;
    while (n != s){
        ListDigraph::Arc arc = path.predArc(n);
        n = path.predNode(n);
        pathCost += getCoeff(arc, 0);
    }
    return pathCost;
}

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

std::string Subgradient::getPathString(Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &s, const ListDigraph::Node &t){
    std::string pathString = "";
    ListDigraph::Node n = t;
    while (n != s){
        pathString += "(" + std::to_string(getNodeLabel(n, 0) + 1) + "," + std::to_string(getNodeSlice(n, 0) + 1) + ")";
        pathString += "-";
        n = path.predNode(n);
    }
    pathString += "(" + std::to_string(getNodeLabel(s, 0) + 1) + "," + std::to_string(getNodeSlice(s, 0) + 1) + ")";
    return pathString;
}
 

void Subgradient::displayMainParameters(Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &s, const ListDigraph::Node &t){
    int k = getIteration();
    int sizeOfField[7];
    sizeOfField[0] = 4;
    sizeOfField[1] = 5;
    sizeOfField[2] = 5;
    sizeOfField[3] = 11;
    sizeOfField[4] = 7;
    sizeOfField[5] = 12;
    sizeOfField[6] = 13;
    char space = ' ';
    std::string field[8];
    if (k == 0){
        field[0] = "It k";
        field[1] = "LB[k]";
        field[2] = "UB[k]";
        field[3] = "Lagr(u[k])";
        field[4] = "u[k]";
        field[5] = "Slack s[k]";
        field[6] = "StepSize t[k]";
        field[7] = "Path P[k]";
        for (int i = 0; i < 7; i++){
            field[i].resize(sizeOfField[i], space);
            std::cout << field[i] << " | ";
        }
        std::cout << field[7] << std::endl;
    }
    field[0] = std::to_string(k);
    field[1] = std::to_string(getLB());
    field[2] = std::to_string(getUB());
    field[3] = std::to_string(getCurrentCost());
    field[4] = std::to_string(getMultiplier_k(k));
    field[5] = std::to_string(getSlack_k(k));
    field[6] = std::to_string(getStepSize_k(k));
    field[7] = getPathString(path, s, t);
    for (int i = 0; i < 7; i++){
        field[i].resize(sizeOfField[i], space);
        std::cout << field[i] << " | ";
    }
    std::cout << field[7] << " | " << std::endl;
}

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
}