#include "lagNonOverlap.h"
#include <set>
#include <lemon/core.h>
#include <lemon/capacity_scaling.h>

/* *******************************************************************************
*                             INITIALIZATION METHODS
******************************************************************************* */
void lagNonOverlap::init(){
    initMultipliers();
    initStabilityCenter();
    initSlacks();
    initSlacks_v2();
    initDirection();
    build_Graph_E();
    initCosts();
    initAssignmentMatrix();
}

/********************************* MULTIPLIERS ***********************************/

/*********************** Option 1: A fixed value, e.g.: 0 ************************/

/* Sets the initial lagrangian multipliers for the subgradient to run. */
void lagNonOverlap::initMultipliers(){

    bool warmstart = getInstance().getInput().getWarmstart();

    if(warmstart){
        /**** Initializing all with 0.0 *****/
        initializeLengthMultipliers(0.0);
        initializeSourceTargetMultipliers(0.0);
        initializeFlowMultipliers(0.0);
        initMultipliersWarmstart();
    }else{
        double initialMultiplier = instance.getInput().getInitialLagrangianMultiplier();
        initializeLengthMultipliers(initialMultiplier);
        initializeSourceTargetMultipliers(initialMultiplier);
        initializeFlowMultipliers(initialMultiplier);
    }

    std::cout << "> Initial Lagrangian multipliers were defined. " << std::endl;
}

/* Sets the initial lagrangian multipliers associated with length constraints. */
void lagNonOverlap::initializeLengthMultipliers(double initialMultiplier){  
    lagrangianMultiplierLength.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        lagrangianMultiplierLength[d] = initialMultiplier;
    }
}

/* Sets the initial lagrangian multipliers associated with source/target constraints*/
void lagNonOverlap::initializeSourceTargetMultipliers(double initialMultiplier){
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
void lagNonOverlap::initializeFlowMultipliers(double initialMultiplier){
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

/********************************* STABILITY CENTER ***********************************/

void lagNonOverlap::initStabilityCenter(){
    initializeLengthSC();
    initializeSourceTargetSC();
    initializeFlowSC();
}

/** Sets the initial lagrangian stability center associated with length constraints. **/
void lagNonOverlap::initializeLengthSC(){
    double initialSC = instance.getInput().getInitialLagrangianMultiplier();
    lagrangianSCLength.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        lagrangianSCLength[d] = initialSC;
    }

}

/** Sets the initial lagrangian stability center associated with Source/Target constraints **/
void lagNonOverlap::initializeSourceTargetSC(){
    double initialSC = instance.getInput().getInitialLagrangianMultiplier();
    lagrangianSCSourceTarget.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        lagrangianSCSourceTarget[d].resize(instance.getNbNodes());
        for (int v = 0; v < instance.getNbNodes(); v++){
            if((v == getToBeRouted_k(d).getSource()) || v == getToBeRouted_k(d).getTarget()) {
                //int load = getToBeRouted_k(d).getLoad();
                //lagrangianSCSourceTarget[d][v]= -load/2;
                lagrangianSCSourceTarget[d][v]= initialSC;
            }
            else{
                lagrangianMultiplierSourceTarget[d][v]= initialSC;
            }
        }
    }
}
        
/** Sets the initial lagrangian stability center associated with flow constraints **/
void lagNonOverlap::initializeFlowSC(){
    double initialSC = instance.getInput().getInitialLagrangianMultiplier();
    lagrangianSCFlow.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        lagrangianSCFlow[d].resize(instance.getNbNodes());
        for (int v = 0; v < instance.getNbNodes(); v++){
            /* The multiplier is not defined for the source and the target */
            if((v == getToBeRouted_k(d).getSource()) || v == getToBeRouted_k(d).getTarget()) {
                lagrangianSCFlow[d][v]= 0.0;
            }else{
                lagrangianSCFlow[d][v]= initialSC;
            }
        }
    }

}

/**************************** Option 2: Warmstart ****************************/

void lagNonOverlap::initMultipliersWarmstart(){
    
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        /** Copying the graph **/
        std::shared_ptr<ListDigraph> vecGraphAux = std::make_shared<ListDigraph>();
        DigraphCopy<ListDigraph, ListDigraph> cg(*vecGraph[d], *vecGraphAux);
        ListDigraph::NodeMap<ListDigraph::Node> nodes(*vecGraph[d]);
        cg.nodeRef(nodes);
        ListDigraph::ArcMap<ListDigraph::Arc> arcs(*vecGraphAux);
        cg.arcCrossRef(arcs);
        //ListDigraph::ArcMap<ListDigraph::Arc> arcs(*vecGraph[d]);
        //cg.arcRef(arcs);
        NodeMap new_map(*vecGraphAux);
        cg.nodeMap(*vecNodeLabel[d], new_map);
        cg.run();
        

        /** Modyfing the graph : contracting all the nodes with the same label **/
        for (int label = 0; label < instance.getNbNodes(); label++){
            if(label != getToBeRouted_k(d).getSource() && label!=  getToBeRouted_k(d).getTarget()){
                ListDigraph::NodeIt previousNode(*vecGraph[d]);
                ListDigraph::Node n = getFirstNodeFromLabel(d, label);
                ListDigraph::NodeIt v(*vecGraph[d]);
                ListDigraph::NodeIt currentNode(*vecGraph[d], v);
                if(n != INVALID){
                    //(*vecNodeSlice[d])[n] = -1;
                    while (v != INVALID){
                        currentNode = v;
                        ListDigraph::NodeIt nextNode(*vecGraph[d], ++currentNode);
                        currentNode = v;
                        if ( (getNodeLabel(v, d) == label) && ((*vecGraph[d]).id(n) != (*vecGraph[d]).id(v)) ){
                            (*vecGraphAux).contract(nodes[n],nodes[v]);
                            //(*vecGraph[d]).contract(n, v);
                        }
                        v = nextNode;
                    }
                }
            }
        }
        /** Creating maps for the COST and CAPACITY**/
        ListDigraph::ArcMap<double> COST(*vecGraphAux);
        ListDigraph::ArcMap<int> CAPACITY(*vecGraphAux);
        ListDigraph::NodeMap<int> SUPPLY(*vecGraphAux);
        for(ListDigraph::ArcIt a(*vecGraphAux); a != INVALID; ++a){
            COST[a] = getCoeff(arcs[a],d);
        }
        for(ListDigraph::ArcIt a(*vecGraphAux); a != INVALID; ++a){
            CAPACITY[a] = 1;
        }
        for(ListDigraph::NodeIt v(*vecGraphAux); v != INVALID; ++v){
            SUPPLY[v] = 0;
        }

        /** Finding the source and destination */
        const ListDigraph::Node SOURCE_ORIGINAL = getFirstNodeFromLabel(d, getToBeRouted_k(d).getSource());
        const ListDigraph::Node TARGET_ORIGINAL = getFirstNodeFromLabel(d, getToBeRouted_k(d).getTarget());
        const ListDigraph::Node SOURCE = nodes[SOURCE_ORIGINAL];
        const ListDigraph::Node TARGET = nodes[TARGET_ORIGINAL];

        SUPPLY[SOURCE] = 1;
        SUPPLY[TARGET] = -1;
        
        /** Running the CapacityScaling problem **/

        //CapacityScaling<ListDigraph,ListDigraph::ArcMap<int>,ListDigraph::ArcMap<double>,ListDigraph::NodeMap<int>> capScale(*vecGraphAux,CAPACITY,COST,SUPPLY);	
        CapacityScaling<ListDigraph,int,double> capScale(*vecGraphAux);//,CAPACITY,COST,SUPPLY);	
        capScale.upperMap(CAPACITY);
        capScale.costMap(COST);
        capScale.supplyMap(SUPPLY);
        capScale.run();

        /** Finding the node potentials **/
        for(ListDigraph::NodeIt v(*vecGraphAux); v != INVALID; ++v){
            int label = new_map[v];
            if(label == getToBeRouted_k(d).getSource() || label == getToBeRouted_k(d).getTarget()){
                lagrangianMultiplierSourceTarget[d][label] = capScale.potential(v);
            }else{
                lagrangianMultiplierFlow[d][label] = capScale.potential(v);
            }

        }

    }    
}

/* Contract nodes with the same given label from the given graph, that is connected with the d-graph. */
void lagNonOverlap::contractNodesFromLabel_v2(std::shared_ptr<ListDigraph> vecGraphAux, ListDigraph::NodeMap<ListDigraph::Node> nodes, int d, int label){
    int nb = 0;
    ListDigraph::NodeIt previousNode(*vecGraph[d]);
    ListDigraph::Node n = getFirstNodeFromLabel(d, label);
    ListDigraph::NodeIt v(*vecGraph[d]);
    ListDigraph::NodeIt currentNode(*vecGraph[d], v);
    if(n != INVALID){
        //(*vecNodeSlice[d])[n] = -1;
        while (v != INVALID){
            currentNode = v;
            ListDigraph::NodeIt nextNode(*vecGraph[d], ++currentNode);
            currentNode = v;
            if ( (getNodeLabel(v, d) == label) && ((*vecGraph[d]).id(n) != (*vecGraph[d]).id(v)) ){
                (*vecGraphAux).contract(nodes[n],nodes[v]);
                //(*vecGraph[d]).contract(n, v);
                nb++;
            }
            v = nextNode;
        }
    }
    //std::cout << "> Number of nodes with label " << label << " contracted: " << nb << std::endl; 
}

/********************************** SLACK ***************************************/

/* Analyse the signal!  -> feasibility
    -> Lenght slack has to be  positive to be feasible
    -> Flow slack has to be zero (equality)                
    -> Target and (if) of source slack has to be zero
    -> (otherwise) source slack has to be psitive
*/

/* Initializes the slack of relaxed constraints. */
void lagNonOverlap::initSlacks(){
    initializeLengthSlacks();
    initializeSourceTargetSlacks();
    initializeFlowSlacks();
    
    std::cout << "> Initial Slacks were defined. " << std::endl;
}

/* Initializes the slack of Length constraints. */ 
void lagNonOverlap::initializeLengthSlacks(){
    lengthSlack.resize(getNbDemandsToBeRouted(), 0.0);
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        //lengthSlack[d] = getToBeRouted_k(d).getMaxLength();
        //lengthSlack[d] = getToBeRouted_k(d).getMaxLength()/100;
        lengthSlack[d] = 1;
    }
}

/* Initializes the slack of Source/Target constraints. */
void lagNonOverlap::initializeSourceTargetSlacks(){
    sourceTargetSlack.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        sourceTargetSlack[d].resize(instance.getNbNodes());
        for(int v = 0; v < instance.getNbNodes();v++){
            sourceTargetSlack[d][v] = 1;
        }
    }
}

/* Initializes the slack of Flow constraints. */
void lagNonOverlap::initializeFlowSlacks(){
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

/********************************** SLACK CONSIDERING PRIMAL VARIABLES ***************************************/

void lagNonOverlap::initSlacks_v2(){
    initializeLengthSlacks_v2();
    initializeSourceTargetSlacks_v2();
    initializeFlowSlacks_v2();
}

void lagNonOverlap::initializeLengthSlacks_v2(){
    lengthSlack_v2.resize(getNbDemandsToBeRouted(), 0.0);
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        //lengthSlack[d] = getToBeRouted_k(d).getMaxLength();
        //lengthSlack[d] = getToBeRouted_k(d).getMaxLength()/100;
        lengthSlack_v2[d] = 1;
    }
}
        
void lagNonOverlap::initializeSourceTargetSlacks_v2(){
    sourceTargetSlack_v2.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        sourceTargetSlack_v2[d].resize(instance.getNbNodes());
        for(int v = 0; v < instance.getNbNodes();v++){
            sourceTargetSlack_v2[d][v] = 1;
        }
    }
}

void lagNonOverlap::initializeFlowSlacks_v2(){
    flowSlack_v2.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        flowSlack_v2[d].resize(instance.getNbNodes());
        for(int v = 0; v < instance.getNbNodes();v++){
            /* The slack is not defined for the source and the target */
            if((v == getToBeRouted_k(d).getSource()) || v == getToBeRouted_k(d).getTarget()) {
                flowSlack_v2[d][v] = 0;
            }else{
                flowSlack_v2[d][v] = 0;
            } 
        }
    }
}

/************************************ DIRECTION **************************************/

void lagNonOverlap::initDirection(){
    lengthDirection.resize(getNbDemandsToBeRouted(), 0.0);
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        lengthDirection[d] = lengthSlack[d];
    }

    sourceTargetDirection.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        sourceTargetDirection[d].resize(instance.getNbNodes());
        for(int v = 0; v < instance.getNbNodes();v++){
            sourceTargetDirection[d][v] = sourceTargetSlack[d][v];
        }
    }

    flowDirection.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        flowDirection[d].resize(instance.getNbNodes());
        for(int v = 0; v < instance.getNbNodes();v++){
            flowDirection[d][v] = flowSlack[d][v];
        }
    }
}

/***************************** BUILDING AUXILIARY GRAPH *******************************/

/* Builds an auxiliary graph the EDGES problems */
void lagNonOverlap::build_Graph_E(){   

    /* Initialization - memory allocation with std::make_shared - it allocates the memory and returs the pointer */
    /* Its a vector of graphs, for each edge, it includes its graph in the list (as push_back but do not need the constructor) */          
    EGraph          = std::make_shared<ListDigraph>();
    ENodeID         = std::make_shared<NodeMap>((*EGraph));
    ENodeDemand     = std::make_shared<NodeMap>((*EGraph));
    ENodeDirection  = std::make_shared<NodeMap>((*EGraph));
    ENodeSlice      = std::make_shared<NodeMap>((*EGraph));
    ENodeFirstConst = std::make_shared<NodeMap>((*EGraph));
    EArcId          = std::make_shared<ArcMap>((*EGraph));
    ECost           = std::make_shared<ArcCost>((*EGraph));
        
    /* The higthes value of slices of an edge */
    int highestSlice = 0;
    for (int e = 0; e < instance.getNbEdges(); e++){
        int nbSlices = instance.getPhysicalLinkFromIndex(e).getNbSlices();
        if(nbSlices > highestSlice){
            highestSlice = nbSlices;
        }
    }

    for(int d = 0; d < getNbDemandsToBeRouted(); d++){
        int load = getToBeRouted_k(d).getLoad();
        for(int slice = load; slice < highestSlice; slice++){
            addENode(d,slice-1,(slice-load+1),-1); 
            addENode(d,slice-1,(slice-load+1),1); 
        }
    }

    /** Adding two auxiliary nodes to do the routing: source and destination  **/
    addENode(-1,-1,-1,0); //source
    addENode(-2,highestSlice,highestSlice,0); //destination

    for(ListDigraph::NodeIt v(*EGraph); v != INVALID; ++v){
        int lastConst = getNodeESlice(v);
        for(ListDigraph::NodeIt v2(*EGraph); v2 != INVALID; ++v2){
            int firstConst = getNodeEFirstConst(v2);
            /** If they are different nodes, and the firstconst the v2 appears is greater than the last const v appears **/
            if((v2 != v) && (firstConst > lastConst)){
                /** Adding an arc between v and v2 **/
                ListDigraph::Arc a = EGraph->addArc(v,v2);
                int id = EGraph->id(a);
                setArcEId(a, id);
                setArcECost(a, __DBL_MAX__);
            }
        }
    }
    /* Sets node index. */
    ENodeIndex      = std::make_shared<NodeMap>((*EGraph)); 
    int index=0;
    for (ListDigraph::NodeIt v(*EGraph); v != INVALID; ++v){
        setNodeEIndex(v, index);
        index++;
    }

    std::cout << "> Graph edge was defined. " << std::endl;
}

/* Function to add a new node to the auxiliary graph */
void lagNonOverlap::addENode(int demand, int slice, int firstConst, int direction){

    ListDigraph::Node node = EGraph->addNode();
    int id = EGraph->id(node);

    setNodeEId(node,id);
    setNodeEDemand(node,demand);
    setNodeEDirection(node,direction);
    setNodeESlice(node,slice);
    setNodeEFirstConst(node,firstConst);

    if(demand == -1){
        setIdSource(id);
    }
    else if(demand == -2){
        setIdDestination(id);
    }
}

/********************************** COST ***************************************/

/* Initializes the costs in the objective function - considering, firstly, the original variables. */
void lagNonOverlap::initCosts(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        cost.emplace_back(std::make_shared<ArcCost>((*vecGraph[d]), 0.0)); 
    }
    updateCosts();
    std::cout << "> Initial costs were defined. " << std::endl;
}

/***************************** ASSIGMENT MATRIX *********************************/

/* Initializes the assignement matrix. */
void lagNonOverlap::initAssignmentMatrix(){
    assignmentMatrix.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        assignmentMatrix[e].resize(countNodes(*EGraph));
        std::fill(assignmentMatrix[e].begin(), assignmentMatrix[e].end(), false);
    }

    assignmentMatrix_d.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        assignmentMatrix_d[d].resize(countArcs(*vecGraph[d]));
        std::fill(assignmentMatrix_d[d].begin(), assignmentMatrix_d[d].end(), false);
    }

    std::cout << "> Initial Assignment matrix was defined. " << std::endl;
}

/* *******************************************************************************
*                             RUNING METHODS
******************************************************************************* */

void lagNonOverlap::run(){
    setCurrentLagrCost(0.0);
    setCurrentRealCost(0.0);
    if(getStatus() == STATUS_FEASIBLE){
        setStatus(STATUS_UNKNOWN);
    }

    const ListDigraph::Node SOURCE = getNodeFromIndex(getIdSource());
    const ListDigraph::Node TARGET = getNodeFromIndex(getIdDestination());

    for (int e = 0; e < instance.getNbEdges(); e++){
        updateCost_E(e);

        /* Solving a shortest path for each edge considering the auxiliary graph */
        /* From the artificial source to the artificial target*/

        //Dijkstra<ListDigraph,ListDigraph::ArcMap<double>> shortestPath((*vecEGraph[e]), (*vecECost[e]));
        //shortestPath.run(SOURCE, TARGET);

        BellmanFord<ListDigraph,ListDigraph::ArcMap<double>> shortestPath((*EGraph), (*ECost));
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
        //createGraphFile(e,0);
    }
    updateAssignment_d();

    updateSlack();
    updateSlack_v2();
    updateDirection();

    int soma = 0;
    for (int ee = 0; ee < instance.getNbEdges(); ee++){
        for(int vv = 0;vv < countNodes(*EGraph); vv++){
            soma += assignmentMatrix[ee][vv];
        }
    }
    //std::cout << soma << std::endl;

    subtractConstantValuesFromLagrCost();

    if (checkFeasibility() == true){
        setStatus(STATUS_FEASIBLE);
    }
}

void lagNonOverlap::subtractConstantValuesFromLagrCost(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        //double val = -getLengthMultiplier_k(d)*getToBeRouted_k(d).getMaxLength();
        //double val = -getLengthMultiplier_k(d)*getToBeRouted_k(d).getMaxLength()/100;
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

/* Checks with the slacks if the solution is feasible. */
bool lagNonOverlap::checkFeasibility(){
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
bool lagNonOverlap::checkLengthFeasibility(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        if (lengthSlack[d] < -DBL_EPSILON){
            return false;
        }
    }
    return true;
}

/*  For the source and destination : its an equality, it must be 0, for the rest it must be non negative*/
bool lagNonOverlap::checkSourceTargetFeasibility(){
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
bool lagNonOverlap::checkFlowFeasibility(){
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

/* *******************************************************************************
*                                GET METHODS
******************************************************************************* */

ListDigraph::Node lagNonOverlap::getNodeFromIndex(int id){
    for (ListDigraph::NodeIt n(*EGraph); n != INVALID; ++n){
        if(getNodeEId(n) == id){
            return n;
        }
    }
    return INVALID;
}

double lagNonOverlap::getRealCostFromPath(int e, BellmanFord< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &SOURCE, const ListDigraph::Node &TARGET){
    double total = 0.0;
    ListDigraph::Node currentNode = path.predNode(TARGET);
    while (currentNode != SOURCE){
        /* Find the correspondent arc (graph for each demand) of the currentNode (auxiliary graph) */
        int d = getNodeEDemand(currentNode);
        int slice = getNodeESlice(currentNode);
        int direction = getNodeEDirection(currentNode);
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            if(getArcSlice(a,d)==slice){
                if((direction==1) && (getNodeLabel((*vecGraph[d]).target(a),d) > getNodeLabel((*vecGraph[d]).source(a),d))){
                    total += getCoeff(a, d); // olhar para pe e p (???)   
                }else if((direction==-1) && (getNodeLabel((*vecGraph[d]).target(a),d) < getNodeLabel((*vecGraph[d]).source(a),d))){
                    total += getCoeff(a, d); // olhar para pe e p (???)   
                }
            }
        }
        currentNode = path.predNode(currentNode);
    }
    return total;
}

/** Returns the constraints slack module **/
double lagNonOverlap::getSlackModule() {
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
    return denominator;
}

/** Returns the constraints direction module **/
double lagNonOverlap::getDirectionModule(){
    double denominator = 0.0;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        //if(!((-getLengthSlack_k(d) < - DBL_EPSILON) && (getLengthMultiplier_k(d) > -DBL_EPSILON && getLengthMultiplier_k(d) < DBL_EPSILON))){
            denominator += std::pow(getLengthDirection_k(d), 2);
        //}
    }
    for(int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){
            if((v == getToBeRouted_k(d).getSource()) || v == getToBeRouted_k(d).getTarget()){   
                //if(!((getSourceTargetSlack_k(d,v) > - DBL_EPSILON && getSourceTargetSlack_k(d,v) < DBL_EPSILON ) && (getSourceTargetMultiplier_k(d,v) > -DBL_EPSILON  && getSourceTargetMultiplier_k(d,v) < DBL_EPSILON))){
                    denominator += std::pow(getSourceTargetDirection_k(d, v), 2);
                //}
            }else{
                //if(!((-getSourceTargetSlack_k(d,v) < - DBL_EPSILON ) &&  (getSourceTargetMultiplier_k(d,v) > -DBL_EPSILON && getSourceTargetMultiplier_k(d,v) < DBL_EPSILON))){
                    denominator += std::pow(getSourceTargetDirection_k(d, v), 2);
                //}
            } 
        }
    }

    for (int d= 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){
            if((v != getToBeRouted_k(d).getSource()) && v != getToBeRouted_k(d).getTarget()){   
                //if(!((getFlowSlack_k(d,v) > - DBL_EPSILON && getFlowSlack_k(d,v) < DBL_EPSILON ) &&  (getFlowMultiplier_k(d,v) > -DBL_EPSILON && getFlowMultiplier_k(d,v) < DBL_EPSILON))){   
                    denominator += std::pow(getFlowDirection_k(d, v), 2);
                //}
            }
        }
    } 
    return denominator;
}

/** Returns the scalar product between the slack (gradient) and the direction **/
double lagNonOverlap::getSlackDirectionProd(){

    Input::ProjectionType projection = getInstance().getInput().getChosenProjection();
    if((projection == Input::IMPROVED) ||(projection == Input::PROJECTED)){
        return getSlackDirectionProdProjected(projection);
    }   
    return getSlackDirectionProdNormal();  
}

double lagNonOverlap::getSlackDirectionProdNormal(){
    double denominator = 0.0;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        denominator += getLengthSlack_k(d)*getLengthDirection_k(d);
    }
    for(int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){
            if((v == getToBeRouted_k(d).getSource()) || v == getToBeRouted_k(d).getTarget()){   
                denominator += getSourceTargetSlack_k(d, v)*getSourceTargetDirection_k(d,v);
            }else{
                denominator += getSourceTargetSlack_k(d, v)*getSourceTargetDirection_k(d,v);
            } 
        }
    }
    for (int d= 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){
            if((v != getToBeRouted_k(d).getSource()) && v != getToBeRouted_k(d).getTarget()){   
                denominator += getFlowSlack_k(d, v)*getFlowDirection_k(d,v);
            }
        }
    }
    return denominator;
}

double lagNonOverlap::getSlackDirectionProdProjected(Input::ProjectionType projection){
    double denominator = 0.0;
    if(projection == Input::IMPROVED){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            if(!((-getLengthDirection_k(d) < - DBL_EPSILON) && (getLengthMultiplier_k(d) > -DBL_EPSILON && getLengthMultiplier_k(d) < DBL_EPSILON))){ // if non negative direction or multiplier different from zero
                denominator += getLengthSlack_k(d)*getLengthDirection_k(d);
            }
        }
        for(int d = 0; d < getNbDemandsToBeRouted(); d++){
            for (int v = 0; v < instance.getNbNodes(); v++){
                if((v == getToBeRouted_k(d).getSource()) || v == getToBeRouted_k(d).getTarget()){   
                    //if(!((getSourceTargetDirection_k(d,v) > - DBL_EPSILON && getSourceTargetDirection_k(d,v) < DBL_EPSILON ) && (getSourceTargetMultiplier_k(d,v) > -DBL_EPSILON  && getSourceTargetMultiplier_k(d,v) < DBL_EPSILON))){
                    denominator += getSourceTargetSlack_k(d, v)*getSourceTargetDirection_k(d,v);
                    //}
                }else{
                    if(!((-getSourceTargetDirection_k(d,v) < - DBL_EPSILON ) &&  (getSourceTargetMultiplier_k(d,v) > -DBL_EPSILON && getSourceTargetMultiplier_k(d,v) < DBL_EPSILON))){ // if non negative direction or multiplier different from zero
                        denominator += getSourceTargetSlack_k(d, v)*getSourceTargetDirection_k(d,v);
                    }
                } 
            }
        }
        for (int d= 0; d < getNbDemandsToBeRouted(); d++){
            for (int v = 0; v < instance.getNbNodes(); v++){
                if((v != getToBeRouted_k(d).getSource()) && v != getToBeRouted_k(d).getTarget()){   
                    //if(!((getFlowDirection_k(d,v) > - DBL_EPSILON && getFlowDirection_k(d,v) < DBL_EPSILON ) &&  (getFlowMultiplier_k(d,v) > -DBL_EPSILON && getFlowMultiplier_k(d,v) < DBL_EPSILON))){   
                    denominator += getFlowSlack_k(d, v)*getFlowDirection_k(d,v);
                    //}
                }
            }
        }
    }else if(projection == Input::PROJECTED){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            if(!(-getLengthDirection_k(d) < - DBL_EPSILON)){ // if non negative direction 
                denominator += getLengthSlack_k(d)*getLengthDirection_k(d);
            }
        }
        for(int d = 0; d < getNbDemandsToBeRouted(); d++){
            for (int v = 0; v < instance.getNbNodes(); v++){
                if((v == getToBeRouted_k(d).getSource()) || v == getToBeRouted_k(d).getTarget()){   
                    //if(!((getSourceTargetDirection_k(d,v) > - DBL_EPSILON && getSourceTargetDirection_k(d,v) < DBL_EPSILON ) && (getSourceTargetMultiplier_k(d,v) > -DBL_EPSILON  && getSourceTargetMultiplier_k(d,v) < DBL_EPSILON))){
                    denominator += getSourceTargetSlack_k(d, v)*getSourceTargetDirection_k(d,v);
                    //}
                }else{
                    if(!(-getSourceTargetDirection_k(d,v) < - DBL_EPSILON )){ // if non negative direction
                        denominator += getSourceTargetSlack_k(d, v)*getSourceTargetDirection_k(d,v);
                    }
                } 
            }
        }
        for (int d= 0; d < getNbDemandsToBeRouted(); d++){
            for (int v = 0; v < instance.getNbNodes(); v++){
                if((v != getToBeRouted_k(d).getSource()) && v != getToBeRouted_k(d).getTarget()){   
                    //if(!((getFlowDirection_k(d,v) > - DBL_EPSILON && getFlowDirection_k(d,v) < DBL_EPSILON ) &&  (getFlowMultiplier_k(d,v) > -DBL_EPSILON && getFlowMultiplier_k(d,v) < DBL_EPSILON))){   
                    denominator += getFlowSlack_k(d, v)*getFlowDirection_k(d,v);
                    //}
                }
            }
        }
    }
    return denominator;
}

/** Returns the constraints slack module (slack considering the primal variables) **/
double lagNonOverlap::getSlackModule_v2(){
    double denominator = 0.0;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        denominator += std::pow(getLengthSlack_v2_k(d), 2);
    }

    for(int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){
            if((v == getToBeRouted_k(d).getSource()) || v == getToBeRouted_k(d).getTarget()){   
                denominator += std::pow(getSourceTargetSlack_v2_k(d, v), 2);
            }else{
                denominator += std::pow(getSourceTargetSlack_v2_k(d, v), 2);
            } 
        }
    }

    for (int d= 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){
            if((v != getToBeRouted_k(d).getSource()) && v != getToBeRouted_k(d).getTarget()){   
                denominator += std::pow(getFlowSlack_v2_k(d, v), 2);
            }
        }
    }

    return denominator;
}

/** Returns the module **/
double lagNonOverlap::getMeanSlackModule_v2(){
    double module = std::sqrt(getSlackModule_v2());
    double numRest = getNbDemandsToBeRouted() + getNbDemandsToBeRouted()*instance.getNbNodes() + (getNbDemandsToBeRouted()*(instance.getNbNodes()-2));
    return module/numRest;
}

/** Returns the scalar product of the normal slack with the slack considering the primal variables **/
double lagNonOverlap::get_prod_slack(){

    double prod = 0.0;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        prod += getLengthSlack_v2_k(d)*getLengthSlack_k(d);
    }

    for(int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){
            if((v == getToBeRouted_k(d).getSource()) || v == getToBeRouted_k(d).getTarget()){   
                prod += getSourceTargetSlack_v2_k(d, v)*getSourceTargetSlack_k(d, v);
            }else{
                prod += getSourceTargetSlack_v2_k(d, v)*getSourceTargetSlack_k(d, v);
            } 
        }
    }

    for (int d= 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){
            if((v != getToBeRouted_k(d).getSource()) && v != getToBeRouted_k(d).getTarget()){   
                prod += getFlowSlack_v2_k(d, v)*getFlowSlack_k(d, v);
            }
        }
    }
    return prod;
}

/* *******************************************************************************
*                             UPDATE METHODS
******************************************************************************* */

/********************************** SLACK ***************************************/

void lagNonOverlap::updateSlack(){
    /* Update Slacks */ 
    updateLengthSlack();
    updateSourceTargetSlack();
    updateFlowSlack();

}
/* b -Ax */
void lagNonOverlap::updateLengthSlack(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double exp = 0.0;
        for(ListDigraph::ArcIt a((*vecGraph[d])); a!= INVALID; ++a){
            int label = getArcLabel(a,d);
            for(ListDigraph::NodeIt n((*EGraph)); n!= INVALID; ++n){
                int index = getNodeEIndex(n);
                int demand = getNodeEDemand(n);
                int slice = getNodeESlice(n);
                int direction = getNodeEDirection(n);
                if((getArcSlice(a,d) == slice) && (d == demand)){
                    if(((direction==1) && (getNodeLabel((*vecGraph[d]).target(a),d) > getNodeLabel((*vecGraph[d]).source(a),d))) || ((direction==-1) && (getNodeLabel((*vecGraph[d]).target(a),d) < getNodeLabel((*vecGraph[d]).source(a),d)))){
                        exp += (getArcLength(a, d)*assignmentMatrix[label][index])/(getToBeRouted_k(d).getMaxLength());
                    }
                }
            }
        }
        //lengthSlack[d] = getToBeRouted_k(d).getMaxLength() - exp;
        //lengthSlack[d] = getToBeRouted_k(d).getMaxLength()/100 - exp;
        //lengthSlack[d] = 0;
        lengthSlack[d] = 1 - exp;
    }
    //std::cout << "\t> Length slack was updated. " << std::endl;   
}

/* b -Ax */
void lagNonOverlap::updateSourceTargetSlack(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){ 
            double exp = 0.0;
            for(ListDigraph::ArcIt a((*vecGraph[d])); a!= INVALID; ++a){
                int label = getArcLabel(a,d);
                for(ListDigraph::NodeIt n((*EGraph)); n!= INVALID; ++n){
                    int index = getNodeEIndex(n);
                    int demand = getNodeEDemand(n);
                    int slice = getNodeESlice(n);
                    int direction = getNodeEDirection(n);
                    if((getArcSlice(a,d) == slice) && (d == demand)){
                        if(((direction==1) && (getNodeLabel((*vecGraph[d]).target(a),d) > getNodeLabel((*vecGraph[d]).source(a),d))) || ((direction==-1) && (getNodeLabel((*vecGraph[d]).target(a),d) < getNodeLabel((*vecGraph[d]).source(a),d)))){
                            if(v == getToBeRouted_k(d).getTarget()){
                                if(getNodeLabel((*vecGraph[d]).target(a),d) == v ){
                                    exp += assignmentMatrix[label][index];
                                }
                            }
                            else{
                                if(getNodeLabel((*vecGraph[d]).source(a),d) == v ){
                                    exp += assignmentMatrix[label][index];
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
void lagNonOverlap::updateFlowSlack(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){
            /* The multiplier is not defined for the source and the target */
            if((v != getToBeRouted_k(d).getSource()) && v != getToBeRouted_k(d).getTarget()) {
                double exp = 0.0;
                for(ListDigraph::ArcIt a((*vecGraph[d])); a!= INVALID; ++a){
                    int label = getArcLabel(a,d);
                    for(ListDigraph::NodeIt n((*EGraph)); n!= INVALID; ++n){
                        int index = getNodeEIndex(n);
                        int demand = getNodeEDemand(n);
                        int slice = getNodeESlice(n);
                        int direction = getNodeEDirection(n);
                        if((getArcSlice(a,d) == slice) && (d == demand)){
                            if(((direction==1) && (getNodeLabel((*vecGraph[d]).target(a),d) > getNodeLabel((*vecGraph[d]).source(a),d))) || ((direction==-1) && (getNodeLabel((*vecGraph[d]).target(a),d) < getNodeLabel((*vecGraph[d]).source(a),d)))){
                                if( getNodeLabel((*vecGraph[d]).source(a),d) == v ){
                                    exp -= assignmentMatrix[label][index];
                                }
                                if( getNodeLabel((*vecGraph[d]).target(a),d) == v ){
                                    exp += assignmentMatrix[label][index];
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


/************************ SLACK CONSIDERING THE PRIMAL VECTOR ********************/

void lagNonOverlap::updateSlack_v2(){
    updateLengthSlack_v2();
    updateSourceTargetSlack_v2();
    updateFlowSlack_v2();
}

void lagNonOverlap::updateLengthSlack_v2(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double exp = 0.0;
        for(ListDigraph::ArcIt a((*vecGraph[d])); a!= INVALID; ++a){
            int index = getArcIndex(a,d);
            exp += (getArcLength(a, d)*primal_linear_solution[d][index])/(getToBeRouted_k(d).getMaxLength());
        }
        lengthSlack_v2[d] = 1 - exp;
    }
}

void lagNonOverlap::updateSourceTargetSlack_v2(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){ 
            double exp = 0.0;
            for(ListDigraph::ArcIt a((*vecGraph[d])); a!= INVALID; ++a){
                int index = getArcIndex(a,d);
                if(v == getToBeRouted_k(d).getTarget()){
                    if(getNodeLabel((*vecGraph[d]).target(a),d) == v ){
                        exp += primal_linear_solution[d][index];
                    }
                }
                else{
                    if(getNodeLabel((*vecGraph[d]).source(a),d) == v ){
                        exp += primal_linear_solution[d][index];
                    }
                }
            }
            sourceTargetSlack_v2[d][v] = 1 - exp;
        }
    }
}

void lagNonOverlap::updateFlowSlack_v2(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){
            /* The multiplier is not defined for the source and the target */
            if((v != getToBeRouted_k(d).getSource()) && v != getToBeRouted_k(d).getTarget()) {
                double exp = 0.0;
                for(ListDigraph::ArcIt a((*vecGraph[d])); a!= INVALID; ++a){
                    int index = getArcIndex(a,d);
                    if( getNodeLabel((*vecGraph[d]).source(a),d) == v ){
                        exp -= primal_linear_solution[d][index];
                    }
                    if( getNodeLabel((*vecGraph[d]).target(a),d) == v ){
                        exp += primal_linear_solution[d][index];
                    }
                }
                flowSlack_v2[d][v] = 0 + exp;
            }
        }
    }

}

/********************************* DIRECTION ***********************************/

void lagNonOverlap::updateDirection(){
    Input::DirectionMethod chosenDirectionMethod = getInstance().getInput().getChosenDirectionMethod();
    double theta= 0.0;;
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
    updateSourceTargetDirection(theta);
    updateFlowDirection(theta);
}

void lagNonOverlap::updateLengthDirection(double theta){

    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double value = getLengthSlack_k(d) + theta*getLengthDirection_k(d);
        setLengthDirection_k(d,value);
    }

}

void lagNonOverlap::updateSourceTargetDirection(double theta){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){ 
            double value = getSourceTargetSlack_k(d,v) + theta*getSourceTargetDirection_k(d,v);
            setSourceTargetDirection_k(d,v,value);
        }
    }
}

void lagNonOverlap::updateFlowDirection(double theta){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){
            /* The multiplier is not defined for the source and the target */
            if((v != getToBeRouted_k(d).getSource()) && v != getToBeRouted_k(d).getTarget()) {
                double value = getFlowSlack_k(d,v) + theta*getFlowDirection_k(d,v);
                setFlowDirection_k(d,v,value);
            }
        }
    }
    
}

/********************************** COSTS ***************************************/

/* Updates the arc costs according to the last lagrangian multiplier available. cost = c + u_k*length */
/* It calcules the cost for the original variables */
void lagNonOverlap::updateCosts(){
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

    
}

/*  Update the cost of the auxiliary graph using the cost of the original variables */
void lagNonOverlap::updateCost_E(int label){

    for(ListDigraph::ArcIt a(*EGraph); a != INVALID; ++a){
        setArcECost(a,__DBL_MAX__);
    }

    /* Passing the costs to the auxiliary graph */
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for(ListDigraph::ArcIt a((*vecGraph[d])); a!= INVALID; ++a){
            if(label == getArcLabel(a,d)){
                for(ListDigraph::NodeIt n((*EGraph)); n!= INVALID; ++n){
                    //int index = getNodeEIndex(n);
                    int demand = getNodeEDemand(n);
                    int slice = getNodeESlice(n);
                    int direction = getNodeEDirection(n);
                    if((getArcSlice(a,d) == slice) && (d == demand)){
                        if(((direction==1) && (getNodeLabel((*vecGraph[d]).target(a),d) > getNodeLabel((*vecGraph[d]).source(a),d))) || ((direction==-1) && (getNodeLabel((*vecGraph[d]).target(a),d) < getNodeLabel((*vecGraph[d]).source(a),d)))){
                            for(ListDigraph::InArcIt arc((*EGraph), n); arc != INVALID; ++arc){
                                setArcECost(arc,getArcCost(a, d));
                            }     
                        }
                    }
                }
            }
        }
    }
    // Include objective function pe
    const ListDigraph::Node TARGET = getNodeFromIndex(getIdDestination());
    for(ListDigraph::InArcIt arc((*EGraph), TARGET); arc != INVALID; ++arc){
        setArcECost(arc,0.0);
    }  
    
}

/********************************** ASSIGNMENT MATRIX ***************************************/

/* Updates the assignment of a edge based on the a given path. */
void lagNonOverlap::updateAssignment_e(int e, BellmanFord< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &SOURCE, const ListDigraph::Node &TARGET){
    std::fill(assignmentMatrix[e].begin(), assignmentMatrix[e].end(), false); 
    ListDigraph::Node currentNode = path.predNode(TARGET);

    //std::cout << "> Assignment Matrix"<< std::endl;

    while (currentNode != SOURCE){
        /* Find the correspondent arc (graph for each demand) of the currentNode (auxiliary graph) */
        int index = getNodeEIndex(currentNode);
        assignmentMatrix[e][index] = true;
        currentNode = path.predNode(currentNode);
        //std::cout << getNodeEDemand(currentNode,e) << std::endl;
    }
    //displayAssignmentMatrix(e);
    //std::cout << "\t> Assigment Matrix was updated. " << std::endl;
}

void lagNonOverlap::updateAssignment_d(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        std::fill(assignmentMatrix_d[d].begin(), assignmentMatrix_d[d].end(), false);
    }

    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for(ListDigraph::ArcIt a((*vecGraph[d])); a!= INVALID; ++a){
            int indexArc = getArcIndex(a,d);
            int label = getArcLabel(a,d);
            for(ListDigraph::NodeIt n((*EGraph)); n!= INVALID; ++n){
                int indexNode = getNodeEIndex(n);
                int demand = getNodeEDemand(n);
                int slice = getNodeESlice(n);
                int direction = getNodeEDirection(n);
                if((getArcSlice(a,d) == slice) && (d == demand)){
                    if(((direction==1) && (getNodeLabel((*vecGraph[d]).target(a),d) > getNodeLabel((*vecGraph[d]).source(a),d))) || ((direction==-1) && (getNodeLabel((*vecGraph[d]).target(a),d) < getNodeLabel((*vecGraph[d]).source(a),d)))){
                        assignmentMatrix_d[d][indexArc] = assignmentMatrix[label][indexNode];
                    }
                }
            }
        }
    }

}

/********************************** MULTIPLIERS ***************************************/

/* Updates lagrangian multiplier with the rule: u[k+1] = u[k] + t[k]*violation */
void lagNonOverlap::updateMultiplier(double step){
    updateLengthMultiplier(step);
    updateSourceTargetMultiplier(step);
    updateFlowMultiplier(step);
}

/** Update length multipliers **/
void lagNonOverlap::updateLengthMultiplier(double step){
    /* length */
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double violation = -getLengthSlack_k(d); /* the constraint is <=, we have to pass to >= as it is a minimization problem*/
        double new_multipliplier = getLengthMultiplier_k(d) + (step*violation);
        setLengthMultiplier_k(d, std::max(new_multipliplier, 0.0)); /* the multiplier is >=0 */
    }
}

/** update source target multipliers **/
void lagNonOverlap::updateSourceTargetMultiplier(double step){
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

void lagNonOverlap::updateFlowMultiplier(double step){
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

/******** MULTIPLIER CONSIDERING THE STABILITY CENTER ************/

void lagNonOverlap::updateMultiplier_v2(double step){
    updateLengthMultiplier_v2(step);
    updateSourceTargetMultiplier_v2(step);
    updateFlowMultiplier_v2(step);
}

void lagNonOverlap::updateLengthMultiplier_v2(double step){
    /* length */
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double violation = -getLengthSlack_v2_k(d); /* the constraint is <=, we have to pass to >= as it is a minimization problem*/
        double new_multipliplier = getLengthSC_k(d) + (step*violation);
        setLengthMultiplier_k(d, std::max(new_multipliplier, 0.0)); /* the multiplier is >=0 */
    }
}

void lagNonOverlap::updateSourceTargetMultiplier_v2(double step){
    /* source target */
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){ 
            if((v == getToBeRouted_k(d).getSource()) || v == getToBeRouted_k(d).getTarget()){
                double violation = - getSourceTargetSlack_v2_k(d,v);
                //std::cout << "violation " << violation << std::endl;
                double new_multipliplier = getSourceTargetSC_k(d,v) + (step*violation); /* equality */
                setSourceTargetMultiplier_k(d,v,new_multipliplier); /* multiplier is a real number*/

            }else{
                double violation = - getSourceTargetSlack_v2_k(d,v); /* the original constraints are <=, we have to change to >= because it is a min problem*/
                double new_multipliplier = getSourceTargetSC_k(d,v) + (step*violation);
                setSourceTargetMultiplier_k(d,v,std::max(new_multipliplier, 0.0)); /* multiplier >= 0*/

            }
        }
    }
}

void lagNonOverlap::updateFlowMultiplier_v2(double step){
    /* flow */
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){
            if((v != getToBeRouted_k(d).getSource()) && v != getToBeRouted_k(d).getTarget()){
                double violation = - getFlowSlack_v2_k(d,v);
                //std::cout << "violation " << violation << std::endl;
                double new_multipliplier = getFlowSC_k(d,v) + (step*violation); /* equality */
                setFlowMultiplier_k(d,v,new_multipliplier); /* multiplier is a real number*/
            }
        }
    }
}

/********************************** STABILITY CENTER ***************************************/
void lagNonOverlap::updateStabilityCenter(){
    updateLengthSC();
    updateSourceTargetSC();
    updateFlowSC();

}

void lagNonOverlap::updateLengthSC(){
    /* length */
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        setLengthSC_k(d, getLengthMultiplier_k(d)); 
    }
}

void lagNonOverlap::updateSourceTargetSC(){
     /* source target */
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){ 
            setSourceTargetSC_k(d,v,getSourceTargetMultiplier_k(d,v)); 
        }
    }

}

void lagNonOverlap::updateFlowSC(){
    /* flow */
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){
            setFlowSC_k(d,v, getFlowMultiplier_k(d,v)); 
        }
    }
}

/* *******************************************************************************
*                                 DISPLAYS
******************************************************************************* */

/* Displays an arc from the graph #e. */
void lagNonOverlap::displayEArc(int e, const ListDigraph::Arc &a){

    /* nodes from the auxiliary graph */
    const ListDigraph::Node SOURCE = (*EGraph).source(a);
    const ListDigraph::Node DEST = (*EGraph).target(a);
        
    if((getNodeEDemand(SOURCE) == -1) && (getNodeEDemand(DEST) ==-2)){
        std::cout << "(Node source: SOURCE) \t (Node destination: DEST)";
        std::cout << "(Cost: " << getArcECost(a) << ")"<< std::endl;

    }else if((getNodeEDemand(SOURCE) == -2) && (getNodeEDemand(DEST) ==-1)){
        std::cout << "(Node source: DEST) \t (Node destination: SOURCE)";
        std::cout << "(Cost: " << getArcECost(a) << ")"<< std::endl;

    }else if(getNodeEDemand(SOURCE) == -1){
        std::cout << "(Node source: SOURCE) \t ";
        int d = getNodeEDemand(DEST);
        std::cout << "(Node destination : ";

        if(getNodeEDirection(DEST) == 1){
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getSource() + 1);
            std::cout  << "--";
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getTarget() + 1);
        }else{
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getTarget() + 1);
            std::cout  << "--";
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getSource() + 1);
        }

        std::cout << ", " << getNodeESlice(DEST) + 1 <<", dem:"<< d +1 <<  ")";
        std::cout << "(Cost: " << getArcECost(a) << ")"<< std::endl;

    }else if(getNodeEDemand(SOURCE) == -2){
        std::cout << "(Node source: DEST) \t ";
        int d = getNodeEDemand(DEST);
        std::cout << "(Node destination : ";

        if(getNodeEDirection(DEST) == 1){
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getSource() + 1);
            std::cout  << "--";
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getTarget() + 1);
        }else{
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getTarget() + 1);
            std::cout  << "--";
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getSource() + 1);
        }
         
        std::cout <<  ", " << getNodeESlice(DEST) + 1 <<", dem:"<< d+1 << ")";
        std::cout << "(Cost: " << getArcECost(a) << ")"<< std::endl;
        
    }else if(getNodeEDemand(DEST) == -1){
        int d = getNodeEDemand(SOURCE);
        std::cout << "(Node source : " ;

        if(getNodeEDirection(SOURCE) == 1){
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getSource() + 1);
            std::cout  << "--";
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getTarget() + 1);
        }else{
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getTarget() + 1);
            std::cout  << "--";
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getSource() + 1);
        }

        std::cout <<   ", " << getNodeESlice(SOURCE) + 1 << ", dem:"<< d+1 << ")\t";
        std::cout << "(Node destination : SOURCE)";
        std::cout << "(Cost: " << getArcECost(a) << ")"<< std::endl;

    }else if(getNodeEDemand(DEST) == -2){
        int d = getNodeEDemand(SOURCE);
        std::cout << "(Node source : ";

        if(getNodeEDirection(SOURCE) == 1){
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getSource() + 1);
            std::cout  << "--";
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getTarget() + 1);
        }else{
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getTarget() + 1);
            std::cout  << "--";
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getSource() + 1);
        }

        std::cout <<  ", " << getNodeESlice(SOURCE) + 1 << ", dem:"<< d+1 << ")\t";
        std::cout << "(Node destination : DEST)";
        std::cout << "(Cost: " << getArcECost(a) << ")"<< std::endl;
            
        
    }else{
        /* Printing the conection */
        int d = getNodeEDemand(SOURCE);
        std::cout << "(Node source : ";

        if(getNodeEDirection(SOURCE) == 1){
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getSource() + 1);
            std::cout  << "--";
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getTarget() + 1);
        }else{
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getTarget() + 1);
            std::cout  << "--";
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getSource() + 1);
        }
        
        std::cout <<  ", " << getNodeESlice(SOURCE) + 1 << ", dem:"<< d+1 <<")\t";

        d = getNodeEDemand(DEST);
        std::cout << "(Node destination : ";
        if(getNodeEDirection(DEST) == 1){
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getSource() + 1);
            std::cout  << "--";
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getTarget() + 1);
        }else{
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getTarget() + 1);
            std::cout  << "--";
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getSource() + 1);
        }
        
        std::cout <<  ", " << getNodeESlice(DEST) + 1 << ", dem:"<< d+1 << ")";
        std::cout << "(Cost: " << getArcECost(a) << ")"<< std::endl;
    }
}

/* Display all arcs from the graph #e. */
void lagNonOverlap::displayEGraph(int e){
    for (ListDigraph::ArcIt a(*EGraph); a != INVALID; ++a){
        displayEArc(e, a);
    }
}

void lagNonOverlap::displayAssignmentMatrix(int e){
    for (ListDigraph::NodeIt v(*EGraph); v != INVALID; ++v){
        displayENode(v,e);
    }
}
/* Displays a node from the graph #d. */
void lagNonOverlap::displayENode(const ListDigraph::Node &n, int e){
    
    if((getNodeEDemand(n) == -1) && (getNodeEDemand(n) ==-2)){
        int d = getNodeEDemand(n);
        int index = getNodeEIndex(n);
        std::cout << "(Node : ";

        if(getNodeEDirection(n) == 1){
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getSource() + 1);
            std::cout  << "--";
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getTarget() + 1);
        }else{
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getTarget() + 1);
            std::cout  << "--";
            std::cout  << (instance.getPhysicalLinkFromIndex(e).getSource() + 1);
        }
        
        std::cout << "," << getNodeESlice(n)+1 << "dem: "<< d<< ")";
        std::cout << "(Assigment: " << assignmentMatrix[e][index] << " )"<<std::endl;
    }
}

void lagNonOverlap::displaySlack(std::ostream & saida){
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

void lagNonOverlap::displayMultiplier(std::ostream & saida){
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

void lagNonOverlap::createGraphFile(int e, int it){
    
    std::string nom = "outputs/graph_edge_";
    nom.append(std::to_string(e+1));
    nom.append("_it_");
    nom.append(std::to_string(it));
    nom.append(".dot");

    std::ofstream fichier(nom);
    if (!fichier.fail()) {
        fichier << "digraph G {" << std::endl;

        for (ListDigraph::ArcIt a(*EGraph); a != INVALID; ++a){

            /* nodes from the auxiliary graph */
            const ListDigraph::Node SOURCE = (*EGraph).source(a);
            const ListDigraph::Node DEST = (*EGraph).target(a);
     
            if((getNodeEDemand(SOURCE) == -1) && (getNodeEDemand(DEST) == -2)){
                fichier << " \" SOURCE \" -> \" DEST \" [label=\"" << getArcECost(a) << "\"];" << std::endl;
            }
            else if((getNodeEDemand(SOURCE) == -2) && (getNodeEDemand(DEST) ==-1)){
                fichier << " \" DEST \" -> \" SOURCE \" [label=\"" << getArcECost(a) << "\"];" << std::endl;
            }
            else if(getNodeEDemand(SOURCE) == -1){
                fichier << " \" SOURCE \" -> \" ";
                int d = getNodeEDemand(DEST);
                if(getNodeEDirection(DEST) == 1){
                    fichier << (instance.getPhysicalLinkFromIndex(e).getSource() + 1);
                    fichier << "--";
                    fichier << (instance.getPhysicalLinkFromIndex(e).getTarget() + 1);
                }else{
                    fichier << (instance.getPhysicalLinkFromIndex(e).getTarget() + 1);
                    fichier << "--";
                    fichier << (instance.getPhysicalLinkFromIndex(e).getSource() + 1);
                }
                fichier << ", ";
                fichier << "s:" << (getNodeESlice(DEST) + 1) << " d:" << (d+1) << " \" " ;
                fichier << "[label=\"" << getArcECost(a) << "\"];" << std::endl;

            }else if(getNodeEDemand(SOURCE) == -2 ){
                fichier << " \" DEST \" -> \" ";
                int d = getNodeEDemand(DEST);
                if(getNodeEDirection(DEST) == 1){
                    fichier << (instance.getPhysicalLinkFromIndex(e).getSource() + 1);
                    fichier << "--";
                    fichier << (instance.getPhysicalLinkFromIndex(e).getTarget() + 1);
                }else{
                    fichier << (instance.getPhysicalLinkFromIndex(e).getTarget() + 1);
                    fichier << "--";
                    fichier << (instance.getPhysicalLinkFromIndex(e).getSource() + 1);
                }
                fichier << ", ";
                fichier << "s:" << (getNodeESlice(DEST) + 1) << " d:" << (d+1) << " \" " ;
                fichier << "[label=\"" << getArcECost(a) << "\"];" << std::endl;
                                          
            }else if(getNodeEDemand(DEST) == -1){
                int d = getNodeEDemand(SOURCE);
                fichier << " \" ";
                if(getNodeEDirection(SOURCE) == 1){
                    fichier << (instance.getPhysicalLinkFromIndex(e).getSource() + 1);
                    fichier << "--";
                    fichier << (instance.getPhysicalLinkFromIndex(e).getTarget() + 1);
                }else{
                    fichier << (instance.getPhysicalLinkFromIndex(e).getTarget() + 1);
                    fichier << "--";
                    fichier << (instance.getPhysicalLinkFromIndex(e).getSource() + 1);
                }
                fichier << ", ";
                fichier << "s:" << getNodeESlice(SOURCE) + 1 << " d:" << d+1;
                fichier << " \" ";
                fichier << " -> ";
                fichier << " \" SOURCE \" ";
                fichier << "[label=\"" << getArcECost(a) << "\"];" << std::endl;

            }else if(getNodeEDemand(DEST) == -2){
                int d = getNodeEDemand(SOURCE);
                fichier << " \" ";
                if(getNodeEDirection(SOURCE) == 1){
                    fichier << (instance.getPhysicalLinkFromIndex(e).getSource() + 1);
                    fichier << "--";
                    fichier << (instance.getPhysicalLinkFromIndex(e).getTarget() + 1);
                }else{
                    fichier << (instance.getPhysicalLinkFromIndex(e).getTarget() + 1);
                    fichier << "--";
                    fichier << (instance.getPhysicalLinkFromIndex(e).getSource() + 1);
                }
                fichier << ", ";
                fichier << "s:" << getNodeESlice(SOURCE) + 1 << " d:" << d+1;
                fichier << " \" ";
                fichier << " -> ";
                fichier << " \" DEST \" ";
                fichier << "[label=\"" << getArcECost(a) << "\"];" << std::endl;
            }else{
                /* Printing the conection */
                int d = getNodeEDemand(SOURCE);
                fichier << " \" ";
                if(getNodeEDirection(SOURCE) == 1){
                    fichier << (instance.getPhysicalLinkFromIndex(e).getSource() + 1);
                    fichier << "--";
                    fichier << (instance.getPhysicalLinkFromIndex(e).getTarget() + 1);
                }else{
                    fichier << (instance.getPhysicalLinkFromIndex(e).getTarget() + 1);
                    fichier << "--";
                    fichier << (instance.getPhysicalLinkFromIndex(e).getSource() + 1);
                }
                fichier << ", ";
                fichier << "s:" << getNodeESlice(SOURCE) + 1 << " d:"<< d+1;
                fichier << " \" ";
                fichier << " -> ";
                fichier << " \" ";

                d = getNodeEDemand(DEST);
                if(getNodeEDirection(DEST) == 1){
                    fichier << (instance.getPhysicalLinkFromIndex(e).getSource() + 1);
                    fichier << "--";
                    fichier << (instance.getPhysicalLinkFromIndex(e).getTarget() + 1);
                }else{
                    fichier << (instance.getPhysicalLinkFromIndex(e).getTarget() + 1);
                    fichier << "--";
                    fichier << (instance.getPhysicalLinkFromIndex(e).getSource() + 1);
                }
                fichier << ", ";
                fichier << "s:" << getNodeESlice(DEST) + 1 << " d:"<< d+1;
                fichier << " \" ";
                fichier << "[label=\"" << getArcECost(a) << "\"];" << std::endl;
            }   
        }

        fichier << "}" << std::endl;
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

lagNonOverlap::~lagNonOverlap(){
    lagrangianMultiplierLength.clear();
    lagrangianMultiplierSourceTarget.clear();
    lagrangianMultiplierFlow.clear();

    lengthSlack.clear();
    sourceTargetSlack.clear();
    flowSlack.clear();

    lengthSlack_v2.clear();
    sourceTargetSlack_v2.clear();
    flowSlack_v2.clear();

    lagrangianSCLength.clear();
    lagrangianSCSourceTarget.clear();
    lagrangianSCFlow.clear();

    lengthDirection.clear();
    sourceTargetDirection.clear();
    flowDirection.clear();

    assignmentMatrix.clear();
    cost.clear();
}
