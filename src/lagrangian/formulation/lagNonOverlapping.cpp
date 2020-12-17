#include "lagNonOverlapping.h"
#include <set>

/* *******************************************************************************
*                             INITIALIZATION METHODS
******************************************************************************* */
void lagNonOverlapping::init(){
    initMultipliers();
    initSlacks();
    build_Graphs_e();
    initCosts();
    initAssignmentMatrix();
}

/********************************* MULTIPLIERS ***********************************/

/* Sets the initial lagrangian multipliers for the subgradient to run. */
void lagNonOverlapping::initMultipliers(){
    initializeLengthMultipliers();
    initializeSourceTargetMultipliers();
    initializeFlowMultipliers();

    std::cout << "> Initial Lagrangian multipliers were defined. " << std::endl;
}

/* Sets the initial lagrangian multipliers associated with length constraints. */
void lagNonOverlapping::initializeLengthMultipliers(){
    double initialMultiplier = instance.getInput().getInitialLagrangianMultiplier();
    lagrangianMultiplierLength.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        lagrangianMultiplierLength[d] = initialMultiplier;
    }
}

/* Sets the initial lagrangian multipliers associated with source/target constraints*/
void lagNonOverlapping::initializeSourceTargetMultipliers(){
    double initialMultiplier = instance.getInput().getInitialLagrangianMultiplier();
    lagrangianMultiplierSourceTarget.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        lagrangianMultiplierSourceTarget[d].resize(instance.getNbNodes());
        for (int v = 0; v < instance.getNbNodes(); v++){
            if((v == getToBeRouted_k(d).getSource()) || v == getToBeRouted_k(d).getTarget()) {
                //int load = getToBeRouted_k(d).getLoad();
                //lagrangianMultiplierSourceTarget[d][v]= -load/2;
                lagrangianMultiplierSourceTarget[d][v]= initialMultiplier;
            }
            else{
                lagrangianMultiplierSourceTarget[d][v]= initialMultiplier;
            }
        }
    }
}

/* Sets the initial lagrangian multipliers associated with flow constraints */
void lagNonOverlapping::initializeFlowMultipliers(){
    double initialMultiplier = instance.getInput().getInitialLagrangianMultiplier();
    lagrangianMultiplierFlow.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        lagrangianMultiplierFlow[d].resize(instance.getNbNodes());
        for (int v = 0; v < instance.getNbNodes(); v++){
            /* The multiplier is not defined for the source and the target */
            if((v == getToBeRouted_k(d).getSource()) || v == getToBeRouted_k(d).getTarget()) {
                lagrangianMultiplierFlow[d][v]= 0.0;
            }else{
                lagrangianMultiplierFlow[d][v]= initialMultiplier;
            }
        }
    }
}

/********************************** SLACK ***************************************/

/* Analyse the signal!  -> feasibility
    -> Lenght slack has to be  positive to be feasible
    -> Flow slack has to be zero (equality)                
    -> Target and (if) of source slack has to be zero
    -> (otherwise) source slack has to be psitive
*/

/* Initializes the slack of relaxed constraints. */
void lagNonOverlapping::initSlacks(){
    initializeLengthSlacks();
    initializeSourceTargetSlacks();
    initializeFlowSlacks();
    
    std::cout << "> Initial Slacks were defined. " << std::endl;
}

/* Initializes the slack of Length constraints. */ 
void lagNonOverlapping::initializeLengthSlacks(){
    lengthSlack.resize(getNbDemandsToBeRouted(), 0.0);
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        //lengthSlack[d] = getToBeRouted_k(d).getMaxLength();
        //lengthSlack[d] = getToBeRouted_k(d).getMaxLength()/100;
        lengthSlack[d] = 1;
    }
}

/* Initializes the slack of Source/Target constraints. */
void lagNonOverlapping::initializeSourceTargetSlacks(){
    sourceTargetSlack.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        sourceTargetSlack[d].resize(instance.getNbNodes());
        for(int v = 0; v < instance.getNbNodes();v++){
            sourceTargetSlack[d][v] = 1;
        }
    }
}

/* Initializes the slack of Flow constraints. */
void lagNonOverlapping::initializeFlowSlacks(){
    flowSlack.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        flowSlack[d].resize(instance.getNbNodes());
        for(int v = 0; v < instance.getNbNodes();v++){
            /* The slack is not defined for the source and the target */
            if((v == getToBeRouted_k(d).getSource()) || v == getToBeRouted_k(d).getTarget()) {
                flowSlack[d][v] = 0;
            }else{
                flowSlack[d][v] = 0;
            } 
        }
    }
}

/***************************** BUILDING AUXILIARY GRAPH *******************************/

/* Builds an auxiliary graph for each edge */
void lagNonOverlapping::build_Graphs_e(){  
    vecESourceIndex.resize(instance.getNbEdges());
    vecEDestinationIndex.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        /* Initialization - memory allocation with std::make_shared - it allocates the memory and returs the pointer */
        /* Its a vector of graphs, for each edge, it includes its graph in the list (as push_back but do not need the constructor) */   
        vecEGraph.emplace_back(std::make_shared<ListDigraph>());
        vecENodeID.emplace_back(std::make_shared<NodeMap>((*vecEGraph[e])));
        vecENodeDemand.emplace_back(std::make_shared<NodeMap>((*vecEGraph[e])));
        vecENodeSlice.emplace_back(std::make_shared<NodeMap>((*vecEGraph[e])));
        vecENodeFirstConst.emplace_back(std::make_shared<NodeMap>((*vecEGraph[e])));
        vecENodeArc.emplace_back(std::make_shared<ListDigraph::NodeMap<std::shared_ptr<ListDigraph::Arc>>>((*vecEGraph[e])));
        vecEArcId.emplace_back(std::make_shared<ArcMap>((*vecEGraph[e])));
        vecECost.emplace_back(std::make_shared<ArcCost>((*vecEGraph[e])));
        
        int linklabel = instance.getPhysicalLinkFromIndex(e).getId(); /* label of link e */
        /** For all arcs in G(D) -> create a corresponding node **/
        for(int d = 0; d < getNbDemandsToBeRouted(); d++){
            int load = getToBeRouted_k(d).getLoad();
            for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){ /* vecGraph is protected in RSA **/
                /** if the arc corresponds to edge e **/  
                if((getArcLabel(a, d)  == linklabel)){
                    int slice = getArcSlice(a,d);
                    /** add a node - label(same as the arc), demand(analysed), slice(same as the arc), firstconstraint(slice-load[k]+1)**/
                    addENode(linklabel,d,slice,(slice-load+1),a); 
                }
            }
        }

        /** Adding two auxiliary nodes to do the routing: source and destination  **/
        addENode(e,-1,-1,-1,INVALID); //source
        addENode(e,-2,instance.getPhysicalLinkFromIndex(e).getNbSlices(),instance.getPhysicalLinkFromIndex(e).getNbSlices(),INVALID); //destination
        
        for(ListDigraph::NodeIt v(*vecEGraph[e]); v != INVALID; ++v){
            int lastConst = getNodeESlice(v,e);
            for(ListDigraph::NodeIt v2(*vecEGraph[e]); v2 != INVALID; ++v2){
                int firstConst = getNodeEFirstConst(v2,e);
                /** If they are different nodes, and the firstconst the v2 appears is greater than the last const v appears **/
                if((v2 != v) && (firstConst > lastConst)){
                    /** Adding an arc between v and v2 **/
                    ListDigraph::Arc a = vecEGraph[e]->addArc(v,v2);
                    int id = vecEGraph[e]->id(a);
                    setArcEId(a, e, id);
                    setArcECost(a, e, 0.0);
                }
            }
        }

        /* Sets node index. */
        vecENodeIndex.emplace_back(new NodeMap((*vecEGraph[e]), -1));
        int index=0;
        for (ListDigraph::NodeIt v(*vecEGraph[e]); v != INVALID; ++v){
            setNodeEIndex(v, e, index);
            index++;
        }
    }
    std::cout << "> Graphs per edge were defined. " << std::endl;
    //createGraphFile(label);
    //displayEGraph(label);
}

/* Function to add a new node to the auxiliary graph */
void lagNonOverlapping::addENode(int e, int demand, int slice, int firstConst,const ListDigraph::Arc & a){
    ListDigraph::Node node = vecEGraph[e]->addNode();
    int id = vecEGraph[e]->id(node);
    setNodeEId(node,e,id);
    setNodeEDemand(node,e,demand);
    setNodeESlice(node,e,slice);
    setNodeEFirstConst(node,e,firstConst);
    setNodeEArc(node,e,a);

    if(demand == -1){
        setIndexSource(id,e);
    }
    else if(demand == -2){
        setIndexDestination(id,e);
    }
}

/********************************** COST ***************************************/

/* Initializes the costs in the objective function - considering, firstly, the original variables. */
void lagNonOverlapping::initCosts(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        cost.emplace_back(std::make_shared<ArcCost>((*vecGraph[d]), 0.0)); 
    }
    updateCosts();
    std::cout << "> Initial costs were defined. " << std::endl;
}

/***************************** ASSIGMENT MATRIX *********************************/

/* Initializes the assignement matrix. */
void lagNonOverlapping::initAssignmentMatrix(){
    assignmentMatrix.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        assignmentMatrix[e].resize(countNodes(*vecEGraph[e]));
        std::fill(assignmentMatrix[e].begin(), assignmentMatrix[e].end(), false);
    }
    std::cout << "> Initial Assignment matrix was defined. " << std::endl;
}

/* *******************************************************************************
*                             RUNING METHODS
******************************************************************************* */

void lagNonOverlapping::run(){
    setCurrentLagrCost(0.0);
    setCurrentRealCost(0.0);
    //if(getStatus() == STATUS_FEASIBLE){
    //    setStatus(STATUS_UNKNOWN);
    //}
    for (int e = 0; e < instance.getNbEdges(); e++){
        const ListDigraph::Node SOURCE = getNodeFromIndex(e, getIndexSource(e));
        //std::cout << "SOurce " << getNodeEDemand(SOURCE,e); //ok
        const ListDigraph::Node TARGET = getNodeFromIndex(e, getIndexDestination(e));
        //std::cout << "DEst " << getNodeEDemand(TARGET,e);
        /* Solving a shortest path for each edge considering the auxiliary graph */
        /* From the artificial source to the artificial target*/
        //Dijkstra<ListDigraph,ListDigraph::ArcMap<double>> shortestPath((*vecEGraph[e]), (*vecECost[e]));
        //shortestPath.run(SOURCE, TARGET);
        BellmanFord<ListDigraph,ListDigraph::ArcMap<double>> shortestPath((*vecEGraph[e]), (*vecECost[e]));
        shortestPath.run(SOURCE);

        /* I think there is always a path, analysing the auxiliary graph */
        if(shortestPath.reached(TARGET) == false){
            setStatus(STATUS_INFEASIBLE);
            std::cout << "> RSA is infeasible because there is no path from the artificial source to the artificial destination. Edge " << e << "." << std::endl;
            return;
        }
        updateAssignment_e(e, shortestPath, SOURCE, TARGET);
        incCurrentLagrCost(shortestPath.dist(TARGET));    
        incCurrentRealCost(getRealCostFromPath(e, shortestPath, SOURCE, TARGET));
    }
    //displayEGraph(0);
    //createGraphFile(0);

    int soma = 0;
    for (int ee = 0; ee < instance.getNbEdges(); ee++){
        for(int vv = 0;vv < countNodes(*vecEGraph[ee]); vv++){
            soma += assignmentMatrix[ee][vv];
        }
    }
    //std::cout << "Sum: " << soma << std::endl;
    /* Update Slacks */ 
    updateLengthSlack();
    updateSourceTargetSlack();
    updateFlowSlack();

    subtractConstantValuesFromLagrCost();

    if (checkFeasibility() == true){
        setStatus(STATUS_FEASIBLE);
    }
}

/* Checks with the slacks if the solution is feasible. */
bool lagNonOverlapping::checkFeasibility(){
    if (checkLengthFeasibility() == false){
        return false;
    }
    if (checkSourceTargetFeasibility() == false){
        return false;
    }
    if (checkFlowFeasibility() == false){
        return false;
    }
    return true;
}

/** Checks if all slacks are non-negative. **/
bool lagNonOverlapping::checkLengthFeasibility(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        if (lengthSlack[d] < -DBL_EPSILON){
            return false;
        }
    }
    return true;
}

/*  For the source and destination : its an equality, it must be 0, for the rest it must be non negative*/
bool lagNonOverlapping::checkSourceTargetFeasibility(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){ 
            if((v == getToBeRouted_k(d).getSource()) ||  (v == getToBeRouted_k(d).getTarget())){
                if((sourceTargetSlack[d][v] < -DBL_EPSILON) || (sourceTargetSlack[d][v] > DBL_EPSILON)){
                    return false;
                }
            }else{
                if(sourceTargetSlack[d][v] < -DBL_EPSILON){
                    return false;
                }
            }
        }
    }
    return true;
}

/* Equalities: it most be 0.*/
bool lagNonOverlapping::checkFlowFeasibility(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){ 
            if((v != getToBeRouted_k(d).getSource()) &&  (v != getToBeRouted_k(d).getTarget())){
                if((flowSlack[d][v] < -DBL_EPSILON) || (flowSlack[d][v] > DBL_EPSILON)){
                    return false;
                }
            }
        }
    }
    return true;
}

void lagNonOverlapping::subtractConstantValuesFromLagrCost(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        //double val = - getLengthMultiplier_k(d)*getToBeRouted_k(d).getMaxLength();
        //double val = - getLengthMultiplier_k(d)*getToBeRouted_k(d).getMaxLength()/100;
        double val = -getLengthMultiplier_k(d);
        incCurrentLagrCost(val);
    }

    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double val = -getSourceTargetMultiplier_k(d,getToBeRouted_k(d).getSource()) - getSourceTargetMultiplier_k(d,getToBeRouted_k(d).getTarget());
        incCurrentLagrCost(val);
    }

    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){ 
            if((v != getToBeRouted_k(d).getSource()) &&  (v != getToBeRouted_k(d).getTarget())){
                double val = -getSourceTargetMultiplier_k(d,v);
                incCurrentLagrCost(val);
            }

        }
    }
}

/* *******************************************************************************
*                                GET METHODS
******************************************************************************* */

ListDigraph::Node lagNonOverlapping::getNodeFromIndex(int e, int id){
    for (ListDigraph::NodeIt n(*vecEGraph[e]); n != INVALID; ++n){
        if(getNodeEId(n,e) == id){
            return n;
        }
    }
    return INVALID;
}

double lagNonOverlapping::getRealCostFromPath(int e, BellmanFord< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &SOURCE, const ListDigraph::Node &TARGET){
    double total = 0.0;
    ListDigraph::Node currentNode = TARGET;
    while (currentNode != SOURCE){
        /* Find the correspondent arc (graph for each demand) of the currentNode (auxiliary graph) */
        const ListDigraph::Arc arc = getNodeEArc(currentNode,e); 
        if(arc != INVALID){ /* do not consider artificial source and destination */
            total += getCoeff(arc, getNodeEDemand(currentNode,e)); // olhar para pe e p (???)   
            //std::cout << getCoeff(arc, getNodeEDemand(currentNode,e)) << std::endl;
        }
        currentNode = path.predNode(currentNode);
    }
    return total;
}

/** Returns the constraints slack module **/
double lagNonOverlapping::getSlackModule(){
    double denominator = 0.0;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        //if(!((-getLengthSlack_k(d) < - DBL_EPSILON) && (getLengthMultiplier_k(d) > -DBL_EPSILON && getLengthMultiplier_k(d) < DBL_EPSILON))){
            denominator += std::pow(getLengthSlack_k(d), 2);
        //}
    }

    for(int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){
            if((v == getToBeRouted_k(d).getSource()) || v == getToBeRouted_k(d).getTarget()){   
                //if(!((getSourceTargetSlack_k(d,v) > - DBL_EPSILON && getSourceTargetSlack_k(d,v) < DBL_EPSILON ) && (getSourceTargetMultiplier_k(d,v) > -DBL_EPSILON  && getSourceTargetMultiplier_k(d,v) < DBL_EPSILON))){
                    denominator += std::pow(getSourceTargetSlack_k(d, v), 2);
                //}
            }else{
                //if(!((-getSourceTargetSlack_k(d,v) < - DBL_EPSILON ) &&  (getSourceTargetMultiplier_k(d,v) > -DBL_EPSILON && getSourceTargetMultiplier_k(d,v) < DBL_EPSILON))){
                    denominator += std::pow(getSourceTargetSlack_k(d, v), 2);
                //}
            } 
        }
    }

    for (int d= 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){
            if((v != getToBeRouted_k(d).getSource()) && v != getToBeRouted_k(d).getTarget()){   
                //if(!((getFlowSlack_k(d,v) > - DBL_EPSILON && getFlowSlack_k(d,v) < DBL_EPSILON ) &&  (getFlowMultiplier_k(d,v) > -DBL_EPSILON && getFlowMultiplier_k(d,v) < DBL_EPSILON))){   
                    denominator += std::pow(getFlowSlack_k(d, v), 2);
                //}
            }
        }
    }
    
    //std::cout << denominator << std::endl;
    return denominator;
}

/* *******************************************************************************
*                             UPDATE METHODS
******************************************************************************* */

/********************************** SLACK ***************************************/
/* b -Ax */
void lagNonOverlapping::updateLengthSlack(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double exp = 0.0;
        for (int e = 0; e < instance.getNbEdges(); e++){
            for(ListDigraph::NodeIt n((*vecEGraph[e])); n!= INVALID; ++n){
                /* for every variable in the auxiliary graph */
                int index = getNodeEIndex(n,e);
                if(getNodeEDemand(n,e) == d){
                    const ListDigraph::Arc arc = getNodeEArc(n,e);
                    /* We do not consider the artificial source and destination */
                    if(arc != INVALID){
                        //exp += getArcLength(arc, d)*assignmentMatrix[e][index];
                        //exp += getArcLength(arc, d)*assignmentMatrix[e][index]/100;
                        exp += getArcLength(arc, d)*assignmentMatrix[e][index]/getToBeRouted_k(d).getMaxLength();
                    }
                }
            }
        }
        //lengthSlack[d] = getToBeRouted_k(d).getMaxLength() - exp;
        //lengthSlack[d] = getToBeRouted_k(d).getMaxLength()/100 - exp;
        lengthSlack[d] = 1 - exp;
        //lengthSlack[d] = 0;
    }

    //std::cout << "\t> Length slack was updated. " << std::endl;
}

/* b -Ax */
void lagNonOverlapping::updateSourceTargetSlack(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){ 
            double exp = 0.0;
            for (int e = 0; e < instance.getNbEdges(); e++){
                for(ListDigraph::NodeIt n((*vecEGraph[e])); n!= INVALID; ++n){
                    /* for every variable in the auxiliary graph */
                    int index = getNodeEIndex(n,e);
                    if(getNodeEDemand(n,e) == d){
                        const ListDigraph::Arc arc = getNodeEArc(n,e);
                        /* We do not consider the artificial source and destination */
                        if(arc != INVALID){
                            if(v == getToBeRouted_k(d).getTarget()){
                                if(getNodeLabel((*vecGraph[d]).target(arc),d) == v ){
                                    exp += assignmentMatrix[e][index];
                                }
                            }
                            else{
                                if(getNodeLabel((*vecGraph[d]).source(arc),d) == v ){
                                    exp += assignmentMatrix[e][index];
                                }
                            }
                        }
                    }
                }
            }
            sourceTargetSlack[d][v] = 1 - exp;
        }
    }
    //std::cout << "\t> Source/Target slack was updated. " << std::endl;
}

/* b -Ax */
void lagNonOverlapping::updateFlowSlack(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){
            /* The multiplier is not defined for the source and the target */
            if((v != getToBeRouted_k(d).getSource()) && v != getToBeRouted_k(d).getTarget()) {
                double exp = 0.0;
                for (int e = 0; e < instance.getNbEdges(); e++){
                    for(ListDigraph::NodeIt n((*vecEGraph[e])); n!= INVALID; ++n){
                        /* for every variable in the auxiliary graph */
                        int index = getNodeEIndex(n,e);
                        if(getNodeEDemand(n,e) == d){
                            const ListDigraph::Arc arc = getNodeEArc(n,e);
                            /* We do not consider the artificial source and destination */
                            if(arc != INVALID){
                                if( getNodeLabel((*vecGraph[d]).source(arc),d) == v ){
                                    exp -= assignmentMatrix[e][index];
                                }
                                if( getNodeLabel((*vecGraph[d]).target(arc),d) == v ){
                                    exp += assignmentMatrix[e][index];
                                }
                            }
                        }
                    }
                }
                flowSlack[d][v] = 0 + exp;
            }
        }
    } 
    //std::cout << "\t> Flow slack was updated. " << std::endl; 
}

/********************************** COSTS ***************************************/

/* Updates the arc costs according to the last lagrangian multiplier available. cost = c + u_k*length */
/* It calcules the cost for the original variables, and then it passes the cost to the auxiliary graph */
void lagNonOverlapping::updateCosts(){
    /********* Original Variables *********/
    /* Original costs */
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            setArcCost(a, d, getCoeff(a, d));
        }
    }

    /* Length multipliers */
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double multiplier =  getLengthMultiplier_k(d);
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            //double incrementValue = multiplier*getArcLength(a, d);
            //double incrementValue = multiplier*getArcLength(a, d)/100;
            double incrementValue = multiplier*getArcLength(a, d)/getToBeRouted_k(d).getMaxLength();
            incArcCost(a, d, incrementValue);
        }
    }

    /* Source/Target multipliers */
    for(int d =0; d < getNbDemandsToBeRouted(); d++){
        for (int label = 0; label < instance.getNbNodes(); label++){
            double multiplier = getSourceTargetMultiplier_k(d,label);
            for (ListDigraph::NodeIt v(*vecGraph[d]); v != INVALID; ++v){
                if((getNodeLabel(v, d) == label) && (label != getToBeRouted_k(d).getTarget()) && (label != getToBeRouted_k(d).getSource())){
                    for (ListDigraph::OutArcIt a((*vecGraph[d]), v); a != INVALID; ++a){
                       double incrementValue = multiplier;
                       incArcCost(a, d, incrementValue);
                    }
                }
                else if((getNodeLabel(v, d) == label) && (label == getToBeRouted_k(d).getSource())){
                    for (ListDigraph::OutArcIt a((*vecGraph[d]), v); a != INVALID; ++a){
                       double incrementValue =  multiplier;
                       incArcCost(a, d, incrementValue);
                    }
                }
                else if((getNodeLabel(v, d) == label) && (label == getToBeRouted_k(d).getTarget())){
                    for (ListDigraph::InArcIt a((*vecGraph[d]), v); a != INVALID; ++a){
                       double incrementValue =  multiplier;
                       incArcCost(a, d, incrementValue);
                    }
                }
            }
        }
    }

    /* Flow constraints */
    for(int d =0; d < getNbDemandsToBeRouted(); d++){
        for (int label = 0; label < instance.getNbNodes(); label++){
            /* It is not defined for the source and destination */
            if((label != getToBeRouted_k(d).getTarget()) && (label != getToBeRouted_k(d).getSource())){
                double multiplier = getFlowMultiplier_k(d,label);
                for (ListDigraph::NodeIt v(*vecGraph[d]); v != INVALID; ++v){
                    if(getNodeLabel(v, d) == label){
                        for (ListDigraph::OutArcIt a((*vecGraph[d]), v); a != INVALID; ++a){
                            double incrementValue =  multiplier;
                            incArcCost(a, d, incrementValue);
                        }
                        for (ListDigraph::InArcIt a((*vecGraph[d]), v); a != INVALID; ++a){
                            double incrementValue = - multiplier;
                            incArcCost(a, d, incrementValue);
                        }
                    }
                }
            }
        }
    }

    /* Passing the costs to the auxiliary graph */
    for (int e = 0; e < instance.getNbEdges(); e++){
        /* for all arcs in the auxiliary graph */
        for(ListDigraph::ArcIt a(*vecEGraph[e]); a != INVALID; ++a){
            const ListDigraph::Node DEST = (*vecEGraph[e]).target(a);
            const ListDigraph::Arc arc = getNodeEArc(DEST,e);
            /* If it is not the final destination */
            if(arc != INVALID){
                setArcECost(a,e,getArcCost(arc, getNodeEDemand(DEST,e)));
            }
        }
    }
}

/********************************** ASSIGNMENT MATRIX ***************************************/

/* Updates the assignment of a edge based on the a given path. */
void lagNonOverlapping::updateAssignment_e(int e, BellmanFord< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &SOURCE, const ListDigraph::Node &TARGET){
    std::fill(assignmentMatrix[e].begin(), assignmentMatrix[e].end(), false); 
    ListDigraph::Node currentNode = path.predNode(TARGET);

    //std::cout << "> Assignment Matrix"<< std::endl;

    while (currentNode != SOURCE){
        /* Find the correspondent arc (graph for each demand) of the currentNode (auxiliary graph) */
        int index = getNodeEIndex(currentNode, e);
        assignmentMatrix[e][index] = true;
        currentNode = path.predNode(currentNode);
        //std::cout << getNodeEDemand(currentNode,e) << std::endl;
    }
    //displayAssignmentMatrix(e);
    //std::cout << "\t> Assigment Matrix was updated. " << std::endl;
}

/********************************** MULTIPLIERS ***************************************/

/* Updates lagrangian multiplier with the rule: u[k+1] = u[k] + t[k]*violation */
void lagNonOverlapping::updateMultiplier(double step){
    updateLengthMultiplier(step);
    updateSourceTargetMultiplier(step);
    updateFlowMultiplier(step);
}

/** Update length multipliers **/
void lagNonOverlapping::updateLengthMultiplier(double step){
    /* length */
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double violation = -getLengthSlack_k(d); /* the constraint is <=, we have to pass to >= as it is a minimization problem*/
        double new_multipliplier = getLengthMultiplier_k(d) + (step*violation);
        setLengthMultiplier_k(d, std::max(new_multipliplier, 0.0)); /* the multiplier is >=0 */
    }
}

/** update source target multipliers **/
void lagNonOverlapping::updateSourceTargetMultiplier(double step){
    /* source target */
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){ 
            if((v == getToBeRouted_k(d).getSource()) || v == getToBeRouted_k(d).getTarget()){
                double violation = - getSourceTargetSlack_k(d,v);
                //std::cout << "violation " << violation << std::endl;
                double new_multipliplier = getSourceTargetMultiplier_k(d,v) + (step*violation); /* equality */
                setSourceTargetMultiplier_k(d,v,new_multipliplier); /* multiplier is a real number*/

            }else{
                double violation = - getSourceTargetSlack_k(d,v); /* the original constraints are <=, we have to change to >= because it is a min problem*/
                double new_multipliplier = getSourceTargetMultiplier_k(d,v) + (step*violation);
                setSourceTargetMultiplier_k(d,v,std::max(new_multipliplier, 0.0)); /* multiplier >= 0*/

            }
        }
    }
}

void lagNonOverlapping::updateFlowMultiplier(double step){
    /* flow */
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){
            if((v != getToBeRouted_k(d).getSource()) && v != getToBeRouted_k(d).getTarget()){
                double violation = - getFlowSlack_k(d,v);
                //std::cout << "violation " << violation << std::endl;
                double new_multipliplier = getFlowMultiplier_k(d,v) + (step*violation); /* equality */
                setFlowMultiplier_k(d,v,new_multipliplier); /* multiplier is a real number*/
            }
        }
    }
}

/* *******************************************************************************
*                                 DISPLAYS
******************************************************************************* */

/* Displays an arc from the graph #e. */
void lagNonOverlapping::displayEArc(int e, const ListDigraph::Arc &a){
    /* nodes from the auxiliary graph */
    const ListDigraph::Node SOURCE = (*vecEGraph[e]).source(a);
    const ListDigraph::Node DEST = (*vecEGraph[e]).target(a);

    /* arcs it represents in the original graph */
    const ListDigraph::Arc ARC_SOURCE = getNodeEArc(SOURCE,e);
    const ListDigraph::Arc ARC_DEST = getNodeEArc(DEST,e);

    /* if its artificial nodes source/ destination in the auxiliary graph */
    //if(getArcECost(a,e) < 0){
    //    std::cout << "edge:" << e;
    if((ARC_SOURCE == INVALID) || (ARC_DEST == INVALID)){
        if((ARC_SOURCE == INVALID) && (ARC_DEST == INVALID)){
            if((getNodeEDemand(SOURCE,e) == -1) && (getNodeEDemand(DEST,e) ==-2)){
                std::cout << "(Node source: SOURCE) \t (Node destination: DEST)";
                std::cout << "(Cost: " << getArcECost(a,e) << ")"<< std::endl;
            }else if((getNodeEDemand(SOURCE,e) == -2) && (getNodeEDemand(DEST,e) ==-1)){
                std::cout << "(Node source: DEST) \t (Node destination: SOURCE)";
                std::cout << "(Cost: " << getArcECost(a,e) << ")"<< std::endl;
            }
        }else if(ARC_SOURCE == INVALID){
            if(getNodeEDemand(SOURCE,e) == -1){
                std::cout << "(Node source: SOURCE) \t ";
                int d = getNodeEDemand(DEST,e);
                std::cout << "(Node destination : " << getNodeLabel((*vecGraph[d]).source(ARC_DEST), d) + 1;
                std::cout << "--";
                std::cout <<  getNodeLabel((*vecGraph[d]).target(ARC_DEST), d) + 1 << ", " << getArcSlice(ARC_DEST, d) + 1 <<", dem:"<< d +1 <<  ")";
                std::cout << "(Cost: " << getArcECost(a,e) << ")"<< std::endl;
            }else if(getNodeEDemand(SOURCE,e) == -2){
                std::cout << "(Node source: DEST) \t ";
                int d = getNodeEDemand(DEST,e);
                std::cout << "(Node destination : " << getNodeLabel((*vecGraph[d]).source(ARC_DEST), d) + 1;
                std::cout << "--";
                std::cout <<  getNodeLabel((*vecGraph[d]).target(ARC_DEST), d) + 1 << ", " << getArcSlice(ARC_DEST, d) + 1 <<", dem:"<< d+1 << ")";
                std::cout << "(Cost: " << getArcECost(a,e) << ")"<< std::endl;
            }
        }else if(ARC_DEST == INVALID){
            if(getNodeEDemand(DEST,e) == -1){
                int d = getNodeEDemand(SOURCE,e);
                std::cout << "(Node source : " << getNodeLabel((*vecGraph[d]).source(ARC_SOURCE) , d) + 1;
                std::cout << "--";
                std::cout <<  getNodeLabel((*vecGraph[d]).target(ARC_SOURCE), d) + 1 << ", " << getArcSlice(ARC_SOURCE, d) + 1 << ", dem:"<< d+1 << ")\t";
                std::cout << "(Node destination : SOURCE)";
                std::cout << "(Cost: " << getArcECost(a,e) << ")"<< std::endl;

            }else if(getNodeEDemand(DEST,e) == -2){
                int d = getNodeEDemand(SOURCE,e);
                std::cout << "(Node source : " << getNodeLabel((*vecGraph[d]).source(ARC_SOURCE) , d) + 1;
                std::cout << "--";
                std::cout <<  getNodeLabel((*vecGraph[d]).target(ARC_SOURCE), d) + 1 << ", " << getArcSlice(ARC_SOURCE, d) + 1 << ", dem:"<< d+1 << ")\t";
                std::cout << "(Node destination : DEST)";
                std::cout << "(Cost: " << getArcECost(a,e) << ")"<< std::endl;
            }
        }
    }else{
        /* Printing the conection */
        int d = getNodeEDemand(SOURCE,e);
        std::cout << "(Node source : " << getNodeLabel((*vecGraph[d]).source(ARC_SOURCE) , d) + 1;
        std::cout << "--";
        std::cout <<  getNodeLabel((*vecGraph[d]).target(ARC_SOURCE), d) + 1 << ", " << getArcSlice(ARC_SOURCE, d) + 1 << ", dem:"<< d+1 <<")\t";

        d = getNodeEDemand(DEST,e);
        std::cout << "(Node destination : " << getNodeLabel((*vecGraph[d]).source(ARC_DEST) , d) + 1;
        std::cout << "--";
        std::cout <<  getNodeLabel((*vecGraph[d]).target(ARC_DEST), d) + 1 << ", " << getArcSlice(ARC_DEST, d) + 1 << ", dem:"<< d+1 << ")";
        std::cout << "(Cost: " << getArcECost(a,e) << ")"<< std::endl;
    }
    //}

}

/* Display all arcs from the graph #e. */
void lagNonOverlapping::displayEGraph(int e){
    for (ListDigraph::ArcIt a(*vecEGraph[e]); a != INVALID; ++a){
        displayEArc(e, a);
    }
}

void lagNonOverlapping::displayAssignmentMatrix(int e){
    for (ListDigraph::NodeIt v(*vecEGraph[e]); v != INVALID; ++v){
        displayENode(v,e);
    }
}
/* Displays a node from the graph #d. */
void lagNonOverlapping::displayENode(const ListDigraph::Node &n, int e){
    const ListDigraph::Arc arc = getNodeEArc(n,e);
    if(arc != INVALID){
        int d = getNodeEDemand(n,e);
        int index = getNodeEIndex(n,e);
        std::cout << "(Node : " << getNodeLabel((*vecGraph[d]).source(arc) , d) + 1;
        std::cout << "--";
        std::cout  << getNodeLabel((*vecGraph[d]).target(arc), d)+1 << "," << getArcSlice(arc, d)+1 << ")";
        std::cout << "(Assigment: " << assignmentMatrix[e][index] << " )"<<std::endl;
    }
}

void lagNonOverlapping::displaySlack(std::ostream & saida){
    std::string display = "Length Slack = [ \n";
    //std::cout << "Length:" << std::endl;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        display += "Demand " + std::to_string(d+1) + ": " +std::to_string(-getLengthSlack_k(d)) + "\n"; 
        //std::cout << "Demand:" << d << " "<< getLengthSlack_k(d) << std::endl;
    }
    display += "]";
    //std::cout << display << std::endl;
    saida << display << std::endl;

    display = "Source Target Slack = [ \n";
    //std::cout << std::endl <<"Source/Target:" << std::endl;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        display += "Demand " + std::to_string(d+1) + "\n ";
        for (int v = 0; v < instance.getNbNodes(); v++){ 
            display += "\t Node " + std::to_string(v+1) + ": "+std::to_string(-getSourceTargetSlack_k(d,v)) + "\n"; 
            //std::cout << "Demand:" << d << " Node:" << v << " "<< getSourceTargetSlack_k(d,v) << std::endl;
        }
    }
    display += "]";
    //std::cout << display << std::endl;
    saida << display << std::endl;

    display = "Flow Slack = [ \n";
    //std::cout << std::endl <<"Source/Target:" << std::endl;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        display += "Demand " + std::to_string(d+1) + "\n ";
        for (int v = 0; v < instance.getNbNodes(); v++){ 
            //if((v != getToBeRouted_k(d).getSource()) && v != getToBeRouted_k(d).getTarget()){
                //std::cout << "Demand:" << d << " Node:" << v << " "<< getFlowSlack_k(d,v) << std::endl;
            //}
            display += "\t Node " + std::to_string(v+1) + ": "+ std::to_string(-getFlowSlack_k(d,v)) + "\n"; 
        }
    }

    display += "]";
    saida << display << std::endl;
    //std::cout << display << std::endl;

}

void lagNonOverlapping::displayMultiplier(std::ostream & saida){
    std::string display = "Length Multiplier = [ \n";
    for (unsigned int i = 0; i < lagrangianMultiplierLength.size(); i++){
        display += "Demand " + std::to_string(i+1) + ": " + std::to_string(getLengthMultiplier_k(i)) + "\n"; 
    }
    display += "]";
    //std::cout << display << std::endl;
    saida << display << std::endl;

    display = "Source Target Multiplier = [ \n";
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        display += "Demand " + std::to_string(d+1) + ":\n ";
        for (int v = 0; v < instance.getNbNodes(); v++){ 
            display += "\t Node " + std::to_string(v+1) + ": "+ std::to_string(getSourceTargetMultiplier_k(d,v)) + "\n"; 

        }
    }
    display += "]";
    //std::cout << display << std::endl;
    saida << display << std::endl;
    //overlap 2
    
    display = "Flow Multiplier = [ \n";
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        display += "Demand " + std::to_string(d+1) + "\n ";
        for (int v = 0; v < instance.getNbNodes(); v++){ 
            display += "\t Node "  + std::to_string(v+1) +": " + std::to_string(getFlowMultiplier_k(d,v)) + "\n"; 

        }
    }
    display += "]";
    //std::cout << display << std::endl;
    saida << display << std::endl;
}

void lagNonOverlapping::createGraphFile(int it){
    //std::cout<<"oi" << std::endl;
    for (int e = 0; e < instance.getNbEdges(); e++){
        std::string nom = "outputs/graph_edge_";
        nom.append(std::to_string(e+1));
        nom.append("_it_");
        nom.append(std::to_string(it));
        nom.append(".dot");
        std::ofstream fichier(nom);
        if (!fichier.fail()) {
            fichier << "digraph G {" << std::endl;
            for (ListDigraph::ArcIt a(*vecEGraph[e]); a != INVALID; ++a){

                /* nodes from the auxiliary graph */
                const ListDigraph::Node SOURCE = (*vecEGraph[e]).source(a);
                const ListDigraph::Node DEST = (*vecEGraph[e]).target(a);

                /* arcs it represents in the original graph */
                const ListDigraph::Arc ARC_SOURCE = getNodeEArc(SOURCE,e);
                const ListDigraph::Arc ARC_DEST = getNodeEArc(DEST,e);

                /* if its artificial nodes source/ destination in the auxiliary graph */
                //if(getArcECost(a,e) < 0){
                //    std::cout << "edge:" << e;
                if((ARC_SOURCE == INVALID) || (ARC_DEST == INVALID)){
                    if((ARC_SOURCE == INVALID) && (ARC_DEST == INVALID)){
                        if((getNodeEDemand(SOURCE,e) == -1) && (getNodeEDemand(DEST,e) ==-2)){
                            fichier << " \" SOURCE \" -> \" DEST \" [label=\"" << getArcECost(a,e) << "\"];" << std::endl;
                        }else if((getNodeEDemand(SOURCE,e) == -2) && (getNodeEDemand(DEST,e) ==-1)){
                            fichier << " \" DEST \" -> \" SOURCE \" [label=\"" << getArcECost(a,e) << "\"];" << std::endl;
                        }
                    }else if(ARC_SOURCE == INVALID){
                        if(getNodeEDemand(SOURCE,e) == -1){
                            int d = getNodeEDemand(DEST,e);
                            fichier << " \" SOURCE \" -> \" ";
                            fichier << (getNodeLabel((*vecGraph[d]).source(ARC_DEST), d) + 1) ;
                            fichier << "--";
                            fichier << (getNodeLabel((*vecGraph[d]).target(ARC_DEST), d) + 1) << ", ";
                            fichier << "s:" << (getArcSlice(ARC_DEST, d) + 1) << " d:" << (d+1) << " \" " ;
                            fichier << "[label=\"" << getArcECost(a,e) << "\"];" << std::endl;
                        }else if(getNodeEDemand(SOURCE,e) == -2){
                            int d = getNodeEDemand(DEST,e);
                            fichier << " \" DEST \" -> \" ";
                            fichier << getNodeLabel((*vecGraph[d]).source(ARC_DEST), d) + 1;
                            fichier << "--";
                            fichier << getNodeLabel((*vecGraph[d]).target(ARC_DEST), d) + 1 << ", ";
                            fichier << "s:" << getArcSlice(ARC_DEST, d) + 1 << " d:" << d+1 << " \" ";
                            fichier << "[label=\"" << getArcECost(a,e) << "\"];" << std::endl;
                        }
                    }else if(ARC_DEST == INVALID){
                        if(getNodeEDemand(DEST,e) == -1){
                            int d = getNodeEDemand(SOURCE,e);
                            fichier << " \" ";
                            fichier << getNodeLabel((*vecGraph[d]).source(ARC_SOURCE) , d) + 1;
                            fichier << "--";
                            fichier << getNodeLabel((*vecGraph[d]).target(ARC_SOURCE), d) + 1;
                            fichier << ", ";
                            fichier << "s:" << getArcSlice(ARC_SOURCE, d) + 1 << " d:" << d+1;
                            fichier << " \" ";
                            fichier << " -> ";
                            fichier << " \" SOURCE \" ";
                            fichier << "[label=\"" << getArcECost(a,e) << "\"];" << std::endl;
                        }else if(getNodeEDemand(DEST,e) == -2){
                            int d = getNodeEDemand(SOURCE,e);
                            fichier << " \" ";
                            fichier << getNodeLabel((*vecGraph[d]).source(ARC_SOURCE) , d) + 1;
                            fichier << "--";
                            fichier << getNodeLabel((*vecGraph[d]).target(ARC_SOURCE), d) + 1;
                            fichier << ", ";
                            fichier << "s:" << getArcSlice(ARC_SOURCE, d) + 1 << " d:" << d+1;
                            fichier << " \" ";
                            fichier << " -> ";
                            fichier << " \" DEST \" ";
                            fichier << "[label=\"" << getArcECost(a,e) << "\"];" << std::endl;
                        }
                    }
                }else{
                    /* Printing the conection */
                    int d = getNodeEDemand(SOURCE,e);
                    fichier << " \" ";
                    fichier << getNodeLabel((*vecGraph[d]).source(ARC_SOURCE) , d) + 1;
                    fichier << "--";
                    fichier << getNodeLabel((*vecGraph[d]).target(ARC_SOURCE), d) + 1;
                    fichier << ", ";
                    fichier << "s:" << getArcSlice(ARC_SOURCE, d) + 1 << " d:"<< d+1;
                    fichier << " \" ";
                    fichier << " -> ";
                    fichier << " \" ";
                    d = getNodeEDemand(DEST,e);
                    fichier << getNodeLabel((*vecGraph[d]).source(ARC_DEST) , d) + 1;
                    fichier << "--";
                    fichier << getNodeLabel((*vecGraph[d]).target(ARC_DEST), d) + 1;
                    fichier << ", ";
                    fichier << "s:" << getArcSlice(ARC_DEST, d) + 1 << " d:"<< d+1;
                    fichier << " \" ";
                    fichier << "[label=\"" << getArcECost(a,e) << "\"];" << std::endl;
                }
            }
            fichier << "}" << std::endl;
        }
        std::string command = "dot -Tjpg -o ";
        std::string nom2 = "outputs/graph_edge_";
        nom2.append(std::to_string(e+1));
        nom2.append("_it_");
        nom2.append(std::to_string(it));
        nom2.append(".jpg");
        command.append(nom2);
        command.append(" ");
        command.append(nom);
        system(command.c_str());
    }
}

/* *******************************************************************************
*                             DESTRUCTOR
******************************************************************************* */

lagNonOverlapping::~lagNonOverlapping(){
    lagrangianMultiplierLength.clear();
    lagrangianMultiplierSourceTarget.clear();
    lagrangianMultiplierFlow.clear();

    lengthSlack.clear();
    sourceTargetSlack.clear();
    flowSlack.clear();

    assignmentMatrix.clear();

    vecEGraph.clear();
    vecENodeID.clear();
    vecENodeDemand.clear();
    vecENodeSlice.clear();
    vecENodeFirstConst.clear();
    vecEArcId.clear();
    vecECost.clear();
    cost.clear();
    vecENodeArc.clear();
    vecESourceIndex.clear();
    vecEDestinationIndex.clear();

}

