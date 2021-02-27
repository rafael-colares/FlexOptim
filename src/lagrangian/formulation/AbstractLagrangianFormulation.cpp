#include "AbstractLagrangianFormulation.h"

/************************************************************************************************************/
/*			                                     CONSTRUCTORS	       	                                    */
/************************************************************************************************************/
AbstractLagFormulation::AbstractLagFormulation(const Instance &instance): FlowForm(instance), time(ClockTime::getTimeNow()){
    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            mapItLabel.emplace_back(std::make_shared<IterableIntMap<ListDigraph, ListDigraph::Node>>((*vecGraph[d])) );
            mapCopy<ListDigraph,NodeMap,IterableIntMap<ListDigraph, ListDigraph::Node>>((*vecGraph[d]),(*vecNodeLabel[d]),(*mapItLabel[d]));
        }
    }
}

/* **************************************************************************************************************
*                                                   GETTERS 
*************************************************************************************************************** */

/* Returns the scalar used to update the direction according to the chosen direction method */
double AbstractLagFormulation::getDirectionMult(){

    Input::DirectionMethod chosenDirectionMethod = getInstance().getInput().getChosenDirectionMethod();
    double theta = 0.0;

    if(chosenDirectionMethod == Input::NORMAL){
        theta = 0.0;
    }
    else if(chosenDirectionMethod == Input::CROWDER){
        theta = getInstance().getInput().getCrowderParameter();
    }
    else if(chosenDirectionMethod == Input::CARMERINI){
        
        if(getSlackDirectionProdNormal() < 0.0){
            double carmerini = getInstance().getInput().getCarmeriniParameter();
            theta = -(carmerini*getSlackDirectionProdNormal())/getDirectionModule();
        }
        else{
            theta = 0.0;
        }
    }
    else if(chosenDirectionMethod == Input::MODIFIED_CARMERINI){
        
        if(getSlackDirectionProdNormal() < 0.0){
            theta = std::sqrt(getSlackModule())/std::sqrt(getDirectionModule());
        }
        else{
            theta = 0.0;
        }
    }
    return theta;
}

double AbstractLagFormulation::initialUBValue(){
    double value = 0.0;
    switch (getInstance().getInput().getChosenObj_k(0)){
        case Input::OBJECTIVE_METRIC_1:
        {
            value = initialUBValueObj1();
            break;
        }
        case Input::OBJECTIVE_METRIC_2:
        {
            value = initialUBValueObj2();
            break;
        }
        case Input::OBJECTIVE_METRIC_4:
        {
            value = initialUBValueObj4();
            break;
        }
        case Input::OBJECTIVE_METRIC_8:
        {
            value = initialUBValueObj8();
            break;
        }
        default:
        {
            std::cerr << "Objective metric out of range.\n";
            exit(0);
            break;
        }
    }
    
    return value;
} 

double AbstractLagFormulation::initialUBValueObj1(){
    double value = getNbDemandsToBeRouted()*getNbSlicesGlobalLimit();
    return value;
}

double AbstractLagFormulation::initialUBValueObj2(){
    double value = 0.0;
    std::cout << "oi" << std::endl;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        const ListDigraph::Node SOURCE = getFirstNodeFromLabel(d, getToBeRouted_k(d).getSource());
        const ListDigraph::Node TARGET = getFirstNodeFromLabel(d, getToBeRouted_k(d).getTarget());

        ScaleMapCost scaleMap((*coeff[d]),-1.0);

        BellmanFord<ListDigraph,ScaleMapCost>  shortestPath((*vecGraph[d]), scaleMap);
        shortestPath.run(SOURCE);
        //Dijkstra<ListDigraph,ScaleMapCost> shortestPath((*vecGraph[d]), scaleMap);
        //shortestPath.run(SOURCE, TARGET);

        double auxiliary = -1.0*shortestPath.dist(TARGET);
        value += auxiliary;
    }
    return value;
}

double AbstractLagFormulation::initialUBValueObj4(){
    double value = 0.0;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        const ListDigraph::Node SOURCE = getFirstNodeFromLabel(d, getToBeRouted_k(d).getSource());
        const ListDigraph::Node TARGET = getFirstNodeFromLabel(d, getToBeRouted_k(d).getTarget());

        ScaleMapCost scaleMap((*coeff[d]),-1.0);

        //BellmanFord<ListDigraph,ScaleMapCost>  shortestPath((*vecGraph[d]), scaleMap);
        //shortestPath.run(SOURCE);
        Dijkstra<ListDigraph,ScaleMapCost> shortestPath((*vecGraph[d]), scaleMap);
        shortestPath.run(SOURCE, TARGET);

        double auxiliary = std::min((-1.0*shortestPath.dist(TARGET)),getToBeRouted_k(d).getMaxLength());
        value += auxiliary;
    }
    return value;
}

double AbstractLagFormulation::initialUBValueObj8(){
    return getNbSlicesGlobalLimit();
}

void AbstractLagFormulation::updateLowerUpperBound(double *lower, double *upper){
    /* Refers to the upper bound of the arcs depending on the chosen objective function. */
    VarMatrix x = getMatrixX();
    operatorLowerUpperBound operLB(x,lower);
    operatorLowerUpperBound operUB(x,upper);
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        operLB.setDemand(d);
        operUB.setDemand(d);
        lowerBound.emplace_back(std::make_shared<ArcMap>((*vecGraph[d])));
        upperBound.emplace_back(std::make_shared<ArcMap>((*vecGraph[d])));
        CombineArcMapArcMapLowerUpperBound combLB((*vecArcIndex[d]),(*vecArcIndex[d]),operLB);
        CombineArcMapArcMapLowerUpperBound combUB((*vecArcIndex[d]),(*vecArcIndex[d]),operUB);
        mapCopy<ListDigraph,CombineArcMapArcMapLowerUpperBound,ArcMap>((*vecGraph[d]),combLB,(*lowerBound[d]));
        mapCopy<ListDigraph,CombineArcMapArcMapLowerUpperBound,ArcMap>((*vecGraph[d]),combUB,(*upperBound[d]));
    }
    if(instance.getInput().isObj8(0)){
        int id = getMaxSliceOverall().getId();
        maxUsedSliceOverallUpperBound = upper[id];
        maxUsedSliceOverallLowerBound = upper[id];
    }

    /*for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            int index = getArcIndex(a);
            int id = x[d][index].getId();
            (*upperBound[d])[a] = upper[id];
            (*lowerBound[d])[a] = upper[id];
        }
    }*/
}

void AbstractLagFormulation::getPrimalSolution(double * colsol){
    std::cout<< "ola" << std::endl;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for(int index = 0; index < assignmentMatrix_d[d].size();index++){
            int id = getVariableX_d(d,index).getId();
            colsol[id] = assignmentMatrix_d[d][index];
        }
    }
    std::cout<< "ola2" << std::endl;
    if(instance.getInput().isObj8(0)){
        int id = getMaxSliceOverall().getId();
        colsol[id] = maxUsedSliceOverall;
    }

}

void AbstractLagFormulation::getPrimalAppSolution(double * colsol){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for(int index = 0; index < primal_linear_solution[d].size();index++){
            int id = getVariableX_d(d,index).getId();
            colsol[id] = primal_linear_solution[d][index];
        }
    }
     if(instance.getInput().isObj8(0)){
        int id = getMaxSliceOverall().getId();
        colsol[id] = maxUsedSliceOverall;
    }
}
/* **************************************************************************************************************
*                                            INITIALIZATION METHODS
*************************************************************************************************************** */

/************************************************ MULTIPLIERS ***************************************************/

/* Sets the initial lagrangian multipliers associated with length constraints. */
void AbstractLagFormulation::initializeLengthMultipliers(double initialMultiplier){  
    lagrangianMultiplierLength.resize(getNbDemandsToBeRouted(),initialMultiplier);
}

/* Sets the initial lagrangian multipliers associated with source/target constraints*/
void AbstractLagFormulation::initializeSourceTargetMultipliers(double initialMultiplier){
    lagrangianMultiplierSourceTarget.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        lagrangianMultiplierSourceTarget[d].resize(instance.getNbNodes(),initialMultiplier);
    }
}

/* Sets the initial lagrangian multipliers associated with flow constraints */
void AbstractLagFormulation::initializeFlowMultipliers(double initialMultiplier){
    lagrangianMultiplierFlow.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        lagrangianMultiplierFlow[d].resize(countNodes(*vecGraph[d]),initialMultiplier);
        int indexSource = getSourceNodeIndex(d);
        int indexTarget = getTargetNodeIndex(d);
        lagrangianMultiplierFlow[d][indexSource] = 0.0;
        lagrangianMultiplierFlow[d][indexTarget] = 0.0;
    }
}

/* Sets the initial lagrangian multipliers associated with non-overlapping constraints. */
void AbstractLagFormulation::initializeOverlapMultipliers(double initialMultiplier){
    lagrangianMultiplierOverlap.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        lagrangianMultiplierOverlap[e].resize(instance.getPhysicalLinkFromIndex(e).getNbSlices(),initialMultiplier);
    }
}

/* Sets the nitial lagrangian multipliers associated with one slice per demand constraints. */
void AbstractLagFormulation::initializeOneSlicePerDemandMultipliers(double initialMultiplier){
    lagrangianMultiplierOneSlicePerDemand.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        lagrangianMultiplierOneSlicePerDemand[e].resize(getNbDemandsToBeRouted(),initialMultiplier);
    }
}

/* Sets the initial lagrangian multipliers associated with maximum slice overall 1 constraints. */
void AbstractLagFormulation::initializeMaxUsedSliceOverallMultipliers(double initialMultiplier){
    lagrangianMultiplierMaxUsedSliceOverall.resize(getNbDemandsToBeRouted(),initialMultiplier);
}

/* Sets the initial lagrangian multipliers associated with maximum slice overall 1 (auxiliary) constraints. */
void AbstractLagFormulation::initializeMaxUsedSliceOverallAuxMultipliers(double initialMultiplier){
    lagrangianMultiplierMaxUsedSliceOverallAux.resize(getNbDemandsToBeRouted(),initialMultiplier);
}

/* Sets the initial lagrangian multipliers associated with maximum slice overall 2 constraints. */
void AbstractLagFormulation::initializeMaxUsedSliceOverall2Multipliers(double initialMultiplier){
    lagrangianMultiplierMaxUsedSliceOverall2.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        lagrangianMultiplierMaxUsedSliceOverall2[e].resize(getNbSlicesLimitFromEdge(e),initialMultiplier);
    }
}

/* Sets the initial lagrangian multipliers associated with maximum slice overall 3 constraints. */
void AbstractLagFormulation::initializeMaxUsedSliceOverall3Multipliers(double initialMultiplier){
    lagrangianMultiplierMaxUsedSliceOverall3.resize(instance.getNbNodes());
    for (int v = 0; v < instance.getNbNodes(); v++){
        lagrangianMultiplierMaxUsedSliceOverall3[v].resize(getNbSlicesGlobalLimit(),initialMultiplier);
    }
}

/************************************************ STABILITY CENTER ***************************************************/

/** Sets the initial lagrangian stability center associated with length constraints. **/
void AbstractLagFormulation::initializeLengthSC(){
    std::copy(lagrangianMultiplierLength.begin(), lagrangianMultiplierLength.end(),std::back_inserter(lagrangianSCLength));
}

/** Sets the initial lagrangian stability center associated with Source/Target constraints **/
void AbstractLagFormulation::initializeSourceTargetSC(){
    lagrangianSCSourceTarget.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        std::copy(lagrangianMultiplierSourceTarget[d].begin(), lagrangianMultiplierSourceTarget[d].end(),std::back_inserter(lagrangianSCSourceTarget[d]));
    }
}
        
/** Sets the initial lagrangian stability center associated with flow constraints **/
void AbstractLagFormulation::initializeFlowSC(){
    lagrangianSCFlow.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        std::copy(lagrangianMultiplierFlow[d].begin(), lagrangianMultiplierFlow[d].end(),std::back_inserter(lagrangianSCFlow[d]));
    }
}

/* Sets the initial lagrangian stability center associated with non-overlapping constraints. */
void AbstractLagFormulation::initializeOverlapSC(){
    lagrangianSCOverlap.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        std::copy(lagrangianMultiplierOverlap[e].begin(), lagrangianMultiplierOverlap[e].end(),std::back_inserter(lagrangianSCOverlap[e]));
    }
}

void AbstractLagFormulation::initializeOneSlicePerDemandSC(){
    lagrangianSCOneSlicePerDemand.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        std::copy(lagrangianMultiplierOneSlicePerDemand[e].begin(), lagrangianMultiplierOneSlicePerDemand[e].end(),std::back_inserter(lagrangianSCOneSlicePerDemand[e]));
    }
}

/* Sets the initial lagrangian stability center associated with maximum slice overall 1 constraints. */
void AbstractLagFormulation::initializeMaxUsedSliceOverallSC(){
    std::copy(lagrangianMultiplierMaxUsedSliceOverall.begin(), lagrangianMultiplierMaxUsedSliceOverall.end(),std::back_inserter(lagrangianSCMaxUsedSliceOverall));
}

/* Sets the initial lagrangian stability center associated with maximum slice overall 1 constraints. */
void AbstractLagFormulation::initializeMaxUsedSliceOverallAuxSC(){
    std::copy(lagrangianMultiplierMaxUsedSliceOverallAux.begin(), lagrangianMultiplierMaxUsedSliceOverallAux.end(),std::back_inserter(lagrangianSCMaxUsedSliceOverallAux));
}

/* Sets the initial lagrangian stability center associated with maximum slice overall 2 constraints. */
void AbstractLagFormulation::initializeMaxUsedSliceOverall2SC(){
    lagrangianSCMaxUsedSliceOverall2.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        std::copy(lagrangianSCMaxUsedSliceOverall2[e].begin(), lagrangianSCMaxUsedSliceOverall2[e].end(),std::back_inserter(lagrangianSCMaxUsedSliceOverall2[e]));
    }
}

/* Sets the initial lagrangian stability center associated with maximum slice overall 3 constraints. */
void AbstractLagFormulation::initializeMaxUsedSliceOverall3SC(){
    lagrangianSCMaxUsedSliceOverall3.resize(instance.getNbNodes());
    for (int v = 0; v < instance.getNbNodes(); v++){
        std::copy(lagrangianSCMaxUsedSliceOverall3[v].begin(), lagrangianSCMaxUsedSliceOverall3[v].end(),std::back_inserter(lagrangianSCMaxUsedSliceOverall3[v]));
    }
}

/**************************************************** SLACKS *******************************************************/

/* Initializes the slack of Length constraints. */ 
void AbstractLagFormulation::initializeLengthSlacks(){
    lengthSlack.resize(getNbDemandsToBeRouted(), 1.0);
}

/* Initializes the slack of Source/Target constraints. */
void AbstractLagFormulation::initializeSourceTargetSlacks(){
    sourceTargetSlack.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        sourceTargetSlack[d].resize(instance.getNbNodes(),1.0);
    }
}

/* Initializes the slack of Flow constraints. */
void AbstractLagFormulation::initializeFlowSlacks(){
    flowSlack.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        flowSlack[d].resize(countNodes(*vecGraph[d]),0.0);
    }
}

/* Initializes the slack of non-overlap constraints. */
void AbstractLagFormulation::initializeOverlapSlacks(){
    overlapSlack.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        overlapSlack[e].resize(instance.getPhysicalLinkFromIndex(e).getNbSlices(),1.0);
    }
}

/* Initializes the slack of one slice per demand constraints. */
void AbstractLagFormulation::initializeOneSlicePerDemandSlacks(){
    oneSlicePerDemandSlack.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        oneSlicePerDemandSlack[e].resize(getNbDemandsToBeRouted(),1.0);
    }
}

/* Initializes the slack of maximum used slice overall 1 constraints. */
void AbstractLagFormulation::initializeMaxUsedSliceOverallSlacks(){
    maxUsedSliceOverallSlack.resize(getNbDemandsToBeRouted(),0.0);
}

/* Initializes the slack of maximum used slice overall 1 (auxiliary) constraints. */
void AbstractLagFormulation::initializeMaxUsedSliceOverallAuxSlacks(){
    maxUsedSliceOverallAuxSlack.resize(getNbDemandsToBeRouted(),0.0);
}

/* Initializes the slack of maximum used slice overall 2 constraints. */
void AbstractLagFormulation::initializeMaxUsedSliceOverall2Slacks(){
    maxUsedSliceOverallSlack2.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        maxUsedSliceOverallSlack2[e].resize(getNbSlicesLimitFromEdge(e),0.0);
    }
}

/* Initializes the slack of maximum used slice overall 3 constraints. */
void AbstractLagFormulation::initializeMaxUsedSliceOverall3Slacks(){
    maxUsedSliceOverallSlack3.resize(instance.getNbNodes());
    for (int v = 0; v < instance.getNbNodes(); v++){
        maxUsedSliceOverallSlack3[v].resize(getNbSlicesGlobalLimit(),0.0);
    }
}

/* Reset the slack of length constraints. */
void AbstractLagFormulation::resetLengthSlacks(){
    std::fill(lengthSlack.begin(),lengthSlack.end(),1.0);
}

/* Reset the slack of source target constraints. */
void AbstractLagFormulation::resetSourceTargetSlacks(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        std::fill(sourceTargetSlack[d].begin(),sourceTargetSlack[d].end(),1.0);
    }
}

/* Reset the slack of source target constraints. */
void AbstractLagFormulation::resetFlowSlacks(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        std::fill(flowSlack[d].begin(),flowSlack[d].end(),0.0);
    } 
}

/* Reset the slack of non-overlap constraints. */
void AbstractLagFormulation::resetOverlapSlacks(){
    for (int e = 0; e < instance.getNbEdges(); e++){
        std::fill(overlapSlack[e].begin(),overlapSlack[e].end(),1.0);
    }
}

/* Resets the slack of one slice per demand constraints. */
void AbstractLagFormulation::resetOneSlicePerDemandSlacks(){
    for (int e = 0; e < instance.getNbEdges(); e++){
        std::fill(oneSlicePerDemandSlack[e].begin(),oneSlicePerDemandSlack[e].end(),1.0);
    }
}

/* Resets the slack of max used slice overall constraints. */
void AbstractLagFormulation::resetMaxUsedSliceOverallSlacks(){
    std::fill(maxUsedSliceOverallSlack.begin(),maxUsedSliceOverallSlack.end(),0.0);
}

/* Resets the slack of max used slice overall (auxiliary) constraints. */
void AbstractLagFormulation::resetMaxUsedSliceOverallAuxSlacks(){
    std::fill(maxUsedSliceOverallAuxSlack.begin(),maxUsedSliceOverallAuxSlack.end(),0.0);
}

/* Resets the slack of max used slice overall 2 constraints. */
void AbstractLagFormulation::resetMaxUsedSliceOverall2Slacks(){
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < getNbSlicesLimitFromEdge(e); s++){
            maxUsedSliceOverallSlack2[e][s]    = maxUsedSliceOverall;
            maxUsedSliceOverallSlack2_v2[e][s] = maxUsedSliceOverall;
        }
    }
}

/* Resets the slack of max used slice overall 3 constraints. */
void AbstractLagFormulation::resetMaxUsedSliceOverall3Slacks(){
    for (int node = 0; node < instance.getNbNodes(); node++){
        for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
            int degree = 0;
            for (ListGraph::NodeIt v(compactGraph); v != INVALID; ++v){
                if ((getCompactNodeLabel(v) == node)){
                    for (ListGraph::IncEdgeIt a(compactGraph, v); a != INVALID; ++a){
                        degree++;
                    }
                }
            }
            maxUsedSliceOverallSlack3[node][s] =     maxUsedSliceOverall*degree;
            maxUsedSliceOverallSlack3_v2[node][s] =  maxUsedSliceOverall*degree;
        }
    }
}


/********************************************* SLACKS CONSIDERING PRIMAL VARIABLES ********************************************/

/* Initializes the slack (primal approximation) of Length constraints. */ 
void AbstractLagFormulation::initializeLengthPrimalSlacks(){
    std::copy(lengthSlack.begin(),lengthSlack.end(),std::back_inserter(lengthSlack_v2));
}

/* Initializes the slack (primal approximation) of Source/Target constraints. */
void AbstractLagFormulation::initializeSourceTargetPrimalSlacks(){
    sourceTargetSlack_v2.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        std::copy(sourceTargetSlack[d].begin(),sourceTargetSlack[d].end(),std::back_inserter(sourceTargetSlack_v2[d]));
    }
}

/* Initializes the slack (primal approximation) of Flow constraints. */
void AbstractLagFormulation::initializeFlowPrimalSlacks(){
    flowSlack_v2.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        std::copy(flowSlack[d].begin(),flowSlack[d].end(),std::back_inserter(flowSlack_v2[d]));
    }
}

/* Initializes the slack (primal approximation) of non-overlap constraints. */
void AbstractLagFormulation::initializeOverlapPrimalSlacks(){
    overlapSlack_v2.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        std::copy(overlapSlack[e].begin(),overlapSlack[e].end(),std::back_inserter(overlapSlack_v2[e]));
    }
}

/* Initializes the slack (primal approximation) of one slice per demand constraints. */
void AbstractLagFormulation::initializeOneSlicePerDemandPrimalSlacks(){
    oneSlicePerDemandSlack_v2.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        std::copy(oneSlicePerDemandSlack[e].begin(),oneSlicePerDemandSlack[e].end(),std::back_inserter(oneSlicePerDemandSlack_v2[e]));
    }
}

/* Initializes the slack (primal approximation) of maximum used slice overall 1 constraints. */
void AbstractLagFormulation::initializeMaxUsedSliceOverallPrimalSlacks(){
    std::copy(maxUsedSliceOverallSlack.begin(),maxUsedSliceOverallSlack.end(),std::back_inserter(maxUsedSliceOverallSlack_v2));
}

/* Initializes the slack (primal approximation) of maximum used slice overall 1 (auxiliary) constraints. */
void AbstractLagFormulation::initializeMaxUsedSliceOverallAuxPrimalSlacks(){
    std::copy(maxUsedSliceOverallAuxSlack.begin(),maxUsedSliceOverallAuxSlack.end(),std::back_inserter(maxUsedSliceOverallAuxSlack_v2));
}

/* Initializes the slack (primal approximation) of maximum used slice overall 2 constraints. */
void AbstractLagFormulation::initializeMaxUsedSliceOverall2PrimalSlacks(){
    maxUsedSliceOverallSlack2_v2.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        std::copy(maxUsedSliceOverallSlack2[e].begin(),maxUsedSliceOverallSlack2[e].end(),std::back_inserter(maxUsedSliceOverallSlack2_v2[e]));
    }
}

/* Initializes the slack (primal approximation) of maximum used slice overall 3 constraints. */
void AbstractLagFormulation::initializeMaxUsedSliceOverall3PrimalSlacks(){
    maxUsedSliceOverallSlack3.resize(instance.getNbNodes());
    for (int v = 0; v < instance.getNbNodes(); v++){
        std::copy(maxUsedSliceOverallSlack3[v].begin(),maxUsedSliceOverallSlack3[v].end(),std::back_inserter(maxUsedSliceOverallSlack3_v2[v]));
    }
}

/************************************************ DIRECTION ***********************************************/

/* Initializes the direction of Length constraints. */ 
void AbstractLagFormulation::initializeLengthDirection(){
    std::copy(lengthSlack.begin(), lengthSlack.end(),std::back_inserter(lengthDirection));
}

/* Initializes the direction of Source/Target constraints. */
void AbstractLagFormulation::initializeSourceTargetDirection(){
    sourceTargetDirection.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        std::copy(sourceTargetSlack[d].begin(), sourceTargetSlack[d].end(),std::back_inserter(sourceTargetDirection[d]));
    }
}

/* Initializes the direction of flow constraints. */
void AbstractLagFormulation::initializeFlowDirection(){
    flowDirection.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        std::copy(flowSlack[d].begin(), flowSlack[d].end(),std::back_inserter(flowDirection[d]));
    }
}

/* Initializes the direction of overlap constraints. */
void AbstractLagFormulation::initializeOverlapDirection(){
    overlapDirection.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        std::copy(overlapSlack[e].begin(), overlapSlack[e].end(),std::back_inserter(overlapDirection[e]));
    }
}

/* Initializes the direction of one slice per demand constraints. */
void AbstractLagFormulation::initializeOneSlicePerDemandDirection(){
    oneSlicePerDemandDirection.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        std::copy(oneSlicePerDemandSlack[e].begin(), oneSlicePerDemandSlack[e].end(),std::back_inserter(oneSlicePerDemandDirection[e]));
    }
}

/* Initializes the direction of maximum used slice overall 1 constraints. */
void AbstractLagFormulation::initializeMaxUsedSliceOverallDirection(){
    std::copy(maxUsedSliceOverallSlack.begin(), maxUsedSliceOverallSlack.end(),std::back_inserter(maxUsedSliceOverallDirection));
}

/* Initializes the direction of maximum used slice overall 1 (auxiliary) constraints. */
void AbstractLagFormulation::initializeMaxUsedSliceOverallAuxDirection(){
    std::copy(maxUsedSliceOverallAuxSlack.begin(), maxUsedSliceOverallAuxSlack.end(),std::back_inserter(maxUsedSliceOverallAuxDirection));
}

/* Initializes the direction of maximum used slice overall 2 constraints. */
void AbstractLagFormulation::initializeMaxUsedSliceOverall2Direction(){
    maxUsedSliceOverall2Direction.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        std::copy(maxUsedSliceOverallSlack2[e].begin(), maxUsedSliceOverallSlack2[e].end(),std::back_inserter(maxUsedSliceOverall2Direction[e]));
    }
}

/* Initializes the direction of maximum used slice overall 3 constraints. */
void AbstractLagFormulation::initializeMaxUsedSliceOverall3Direction(){
    maxUsedSliceOverall3Direction.resize(instance.getNbNodes());
    for (int v = 0; v < instance.getNbNodes(); v++){
        std::copy(maxUsedSliceOverallSlack3[v].begin(), maxUsedSliceOverallSlack3[v].end(),std::back_inserter(maxUsedSliceOverall3Direction[v]));
    }
}

/************************************************** ASSIGNMENT MATRIX *************************************************/

/* Initializes the assignment matrix (the sub problem solution) and the maxUsedSlice Overall **/
/* The assignment matrix represents the x variables of the model */
void AbstractLagFormulation::initAssignmentMatrix(){
    assignmentMatrix_d.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        assignmentMatrix_d[d].resize(countArcs(*vecGraph[d]),false);
    }
    maxUsedSliceOverall = 0.0;
    varAuxZ.resize(getNbDemandsToBeRouted(),false);

    std::cout << "> Initial Assignment matrix was defined. " << std::endl;
}

/************************************************** PRIMAL SOLUTION **************************************************/

/* Initializes the primal solution - primal approximation given by the volume method */
void AbstractLagFormulation::initPrimalSolution(){
    primal_linear_solution.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        std::copy(assignmentMatrix_d[d].begin(),assignmentMatrix_d[d].end(),std::back_inserter(primal_linear_solution[d]));   
    } 

    primalVarAuxZ.resize(getNbDemandsToBeRouted());
}

/************************************************** PRIMAL APPROXIMATION **********************************************/

void AbstractLagFormulation::initPrimalApproximation(){
    initStabilityCenter();
    initPrimalSolution();
    initPrimalSlacks();
    setCurrentPrimalCost(getLagrCurrentCost());
}

/****************************************************** COEFF *********************************************************/

/* Initializes the coefficients in the objective function ( according to chosen objective function). */
void AbstractLagFormulation::initCoeff(){
    operatorSliceCoefficient operSliceCoeff; // operator to modify keep the slice value just when arc leaves the source
    operatorCoefficient operCoeff(getInstance().getInput().getChosenObj_k(0)); // operator to return the objective function coefficient
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        operSliceCoeff.setSource(getFirstNodeFromLabel(d, getToBeRouted_k(d).getSource()));
        operCoeff.setDemandLoad(getToBeRouted_k(d).getLoad());
        SourceMap<ListDigraph> sourceMap((*vecGraph[d])); // Map of the source nodes of the arcs
        CombineMapCoeffSlice combineMapCS((*vecArcSlice[d]),sourceMap,operSliceCoeff); 
        CombineMapCoeff combineMap((*vecArcLength[d]),combineMapCS,operCoeff);
        coeff.emplace_back(std::make_shared<ArcCost>((*vecGraph[d]), 0.0)); 
        mapCopy<ListDigraph,CombineMapCoeff,ArcCost>((*vecGraph[d]),combineMap,(*coeff[d])); // Copy to standard map to an easier utilization
    }
    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        setCoeffMapObj8();
    } 
    std::cout << "> Initial Coeffs were defined. " << std::endl;
}

void AbstractLagFormulation::setCoeffMapObj8() {
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        coeff8.emplace_back(std::make_shared<ArcCost>((*vecGraph[d]), 0.0)); 

        operatorSliceCoefficient operSliceCoeff; 
        operSliceCoeff.setSource(getFirstNodeFromLabel(d, getToBeRouted_k(d).getSource()));
        SourceMap<ListDigraph> sourceMap((*vecGraph[d])); 
        CombineMapCoeffSlice combineMapCS((*vecArcSlice[d]),sourceMap,operSliceCoeff);
        mapCopy<ListDigraph,CombineMapCoeffSlice,ArcCost>((*vecGraph[d]),combineMapCS,(*coeff8[d]));
    }
}

/* ****************************************************************************************************************
*                                                   UPDATE METHODS
**************************************************************************************************************** */

/*************************************************** MULTIPLIERS **************************************************/

/** Updates length multipliers **/
void AbstractLagFormulation::updateLengthMultiplier(double step){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double violation = -getLengthDirection_k(d); /* the constraint is <=, we have to pass to >= as it is a minimization problem*/
        double new_multipliplier = getLengthMultiplier_k(d) + (step*violation);
        setLengthMultiplier_k(d, std::max(new_multipliplier, 0.0)); /* the multiplier is >=0 */
    }
}

/** Updates source target multipliers **/
void AbstractLagFormulation::updateSourceTargetMultiplier(double step){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){ 
            if((v == getToBeRouted_k(d).getSource()) || v == getToBeRouted_k(d).getTarget()){
                double violation = -getSourceTargetDirection_k(d,v);
                double new_multipliplier = getSourceTargetMultiplier_k(d,v) + (step*violation); /* equality */
                setSourceTargetMultiplier_k(d,v,new_multipliplier); /* multiplier is a real number*/
            }else{
                double violation = -getSourceTargetDirection_k(d,v); /* the original constraints are <=, we have to change to >= because it is a min problem*/
                double new_multipliplier = getSourceTargetMultiplier_k(d,v) + (step*violation);
                setSourceTargetMultiplier_k(d,v,std::max(new_multipliplier, 0.0)); /* multiplier >= 0*/
            }
        }
    }
}

/** Updates flow multipliers **/
void AbstractLagFormulation::updateFlowMultiplier(double step){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < countNodes(*vecGraph[d]); v++){
            if(v != getTargetNodeIndex(d) && v != getSourceNodeIndex(d)){
                double violation = -getFlowDirection_k(d,v);
                double new_multiplier = getFlowMultiplier_k(d,v) + (step*violation); /* equality */
                setFlowMultiplier_k(d,v,new_multiplier); /* multiplier is a real number*/
                //std::cout << "Flow multiplier " << d << " " << v << " " <<  new_multipliplier << std::endl;
            }
        }
    }
}

/** update overlap multipliers **/
void AbstractLagFormulation::updateOverlapMultiplier(double step){
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            double violation = -getOverlapDirection_k( e, s);
            double new_multipliplier = getOverlapMultiplier_k( e, s) + (step*violation);
            setOverlapMultiplier_k( e, s, std::max(new_multipliplier, 0.0));
        }
    }  
}

/** update one slice per demand multipliers **/
void AbstractLagFormulation::updateOneSlicePerDemandMultiplier(double step){
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            double violation = -getOneSlicePerDemandDirection_k(e,d);
            double new_multipliplier = getOneSlicePerDemandMultiplier_k(e,d) +(step*violation);
            setOneSlicePerDemandMultiplier_k(e,d,std::max(new_multipliplier, 0.0));
        }
    }
}

/** update maximum used slice overall 1 multipliers **/
void AbstractLagFormulation::updateMaxUsedSliceOverallMultiplier(double step){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double violation = -getMaxUsedSliceOverallDirection_k(d);
        double new_multipliplier = getMaxUsedSliceOverallMultiplier_k(d) + (step*violation);
        setMaxUsedSliceOverallMultiplier_k(d,std::max(new_multipliplier, 0.0));
    }
}

/** update maximum used slice overall 1 (auxiliary) multipliers **/
void AbstractLagFormulation::updateMaxUsedSliceOverallAuxMultiplier(double step){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double violation = -getMaxUsedSliceOverallAuxDirection_k(d);
        double new_multipliplier = getMaxUsedSliceOverallAuxMultiplier_k(d) + (step*violation);
        setMaxUsedSliceOverallAuxMultiplier_k(d,std::max(new_multipliplier, 0.0));
    }
}

/** update maximum used slice overall 2 multipliers **/
void AbstractLagFormulation::updateMaxUsedSliceOverall2Multiplier(double step){
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < getNbSlicesLimitFromEdge(e); s++){
            double violation = -getMaxUsedSliceOverall2Direction_k(e,s);
            double new_multipliplier = getMaxUsedSliceOverall2Multiplier_k(e,s) + (step*violation);
            setMaxUsedSliceOverall2Multiplier_k(e,s,std::max(new_multipliplier, 0.0)); 
        }
    }
}

/** update maximum used slice overall 3 multipliers **/
void AbstractLagFormulation::updateMaxUsedSliceOverall3Multiplier(double step){
    for (int v = 0; v < instance.getNbNodes(); v++){
        for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
            double violation = -getMaxUsedSliceOverall3Direction_k(v,s);
            double new_multipliplier = getMaxUsedSliceOverall3Multiplier_k(v,s) + (step*violation);
            setMaxUsedSliceOverall3Multiplier_k(v,s,std::max(new_multipliplier, 0.0)); 
        }
    }
}

/********************************** MULTIPLIER CONSIDERING THE STABILITY CENTER *********************************/

/** update length multipliers considering stability center - volume**/
void AbstractLagFormulation::updateLengthMultiplier_v2(double step){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double violation = -getLengthSlack_v2_k(d); /* the constraint is <=, we have to pass to >= as it is a minimization problem*/
        double new_multipliplier = getLengthSC_k(d) + (step*violation);
        setLengthMultiplier_k(d, std::max(new_multipliplier, 0.0)); /* the multiplier is >=0 */
    }
}

/** update source/target multipliers considering stability center - volume**/
void AbstractLagFormulation::updateSourceTargetMultiplier_v2(double step){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){ 
            if((v == getToBeRouted_k(d).getSource()) || v == getToBeRouted_k(d).getTarget()){
                double violation = -getSourceTargetSlack_v2_k(d,v);
                double new_multipliplier = getSourceTargetSC_k(d,v) + (step*violation); /* equality */
                setSourceTargetMultiplier_k(d,v,new_multipliplier); /* multiplier is a real number*/

            }else{
                double violation = -getSourceTargetSlack_v2_k(d,v); /* the original constraints are <=, we have to change to >= because it is a min problem*/
                double new_multipliplier = getSourceTargetSC_k(d,v) + (step*violation);
                setSourceTargetMultiplier_k(d,v,std::max(new_multipliplier, 0.0)); /* multiplier >= 0*/
            }
        }
    }
}

/** update flow multipliers considering stability center - volume**/
void AbstractLagFormulation::updateFlowMultiplier_v2(double step){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < countNodes(*vecGraph[d]); v++){
            if(v != getTargetNodeIndex(d) && v != getSourceNodeIndex(d)){
                double violation = -getFlowSlack_v2_k(d,v);
                double new_multipliplier = getFlowSC_k(d,v) + (step*violation); /* equality */
                setFlowMultiplier_k(d,v,new_multipliplier); /* multiplier is a real number*/
            }
        }
    }
}

/** Updates overlap  multipliers considering stability center - volume**/
void AbstractLagFormulation::updateOverlapMultiplier_v2(double step){
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            double violation = -getOverlapSlack_v2_k( e, s);
            double new_multipliplier = getOverlapSC_k( e, s) + (step*violation);
            setOverlapMultiplier_k( e, s, std::max(new_multipliplier, 0.0));
        }
    }  
}

/** Updates one slice per demand  multipliers considering stability center - volume**/
void AbstractLagFormulation::updateOneSlicePerDemandMultiplier_v2(double step){
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            double violation = -getOneSlicePerDemandSlack_v2_k(e,d);
            double new_multipliplier = getOneSlicePerDemandSC_k(e,d) + (step*violation);
            setOneSlicePerDemandMultiplier_k(e,d,std::max(new_multipliplier,0.0));
        }
    }
}

/** Updates maximum used slice overall 1 multipliers considering stability center - volume**/
void AbstractLagFormulation::updateMaxUsedSliceOverallMultiplier_v2(double step){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double violation = -getMaxUsedSliceOverallSlack_v2_k(d);
        double new_multipliplier = getMaxUsedSliceOverallSC_k(d) + (step*violation);
        setMaxUsedSliceOverallMultiplier_k(d,std::max(new_multipliplier, 0.0));
    }
}

/** Updates maximum used slice overall 1 (auxiliary) multipliers considering stability center - volume**/
void AbstractLagFormulation::updateMaxUsedSliceOverallAuxMultiplier_v2(double step){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double violation = -getMaxUsedSliceOverallAuxSlack_v2_k(d);
        double new_multipliplier = getMaxUsedSliceOverallAuxSC_k(d) + (step*violation);
        setMaxUsedSliceOverallAuxMultiplier_k(d,std::max(new_multipliplier, 0.0));
    }
}

/** Updates maximum used slice overall 2 multipliers considering stability center - volume**/
void AbstractLagFormulation::updateMaxUsedSliceOverall2Multiplier_v2(double step){
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < getNbSlicesLimitFromEdge(e); s++){
            double violation = -getMaxUsedSliceOverall2Slack_v2_k(e,s);
            double new_multipliplier = getMaxUsedSliceOverall2SC_k(e,s) + (step*violation);
            setMaxUsedSliceOverall2Multiplier_k(e,s,std::max(new_multipliplier, 0.0)); 
        }
    }
}

/** Updates maximum used slice overall 3 multipliers considering stability center - volume**/
void AbstractLagFormulation::updateMaxUsedSliceOverall3Multiplier_v2(double step){
    for (int v = 0; v < instance.getNbNodes(); v++){
        for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
            double violation = -getMaxUsedSliceOverall3Slack_v2_k(v,s);
            double new_multipliplier = getMaxUsedSliceOverall3SC_k(v,s) + (step*violation);
            setMaxUsedSliceOverall3Multiplier_k(v,s,std::max(new_multipliplier, 0.0)); 
        }
    }
}

/************************************************ STABILITY CENTER ***************************************************/

/* Updates length stability center */
void AbstractLagFormulation::updateLengthSC(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        setLengthSC_k(d, getLengthMultiplier_k(d)); 
    }
}

/* Updates source/target stability center */
void AbstractLagFormulation::updateSourceTargetSC(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){ 
            setSourceTargetSC_k(d,v,getSourceTargetMultiplier_k(d,v)); 
        }
    }
}

/* Updates flow stability center */
void AbstractLagFormulation::updateFlowSC(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < countNodes(*vecGraph[d]); v++){
            setFlowSC_k(d,v, getFlowMultiplier_k(d,v)); 
        }
    }
}

/* Updates overlap stability center */
void AbstractLagFormulation::updateOverlapSC(){
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            setOverlapSC_k( e, s, getOverlapMultiplier_k( e, s));
        }
    }  
}

/* Updates one slice per demand stability center */
void AbstractLagFormulation::updateOneSlicePerDemandSC(){
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            setOneSlicePerDemandSC_k(e,d,getOneSlicePerDemandMultiplier_k(e,d));
        }
    }
}

/* Updates maximum used slice overall 1 stability center */
void AbstractLagFormulation::updateMaxUsedSliceOverallSC(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        setMaxUsedSliceOverallSC_k(d,getMaxUsedSliceOverallMultiplier_k(d));
    }
}

/* Updates maximum used slice overall 1 (auxiliary) stability center */
void AbstractLagFormulation::updateMaxUsedSliceOverallAuxSC(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        setMaxUsedSliceOverallAuxSC_k(d,getMaxUsedSliceOverallAuxMultiplier_k(d));
    }
}

/* Updates maximum used slice overall 2 stability center */
void AbstractLagFormulation::updateMaxUsedSliceOverall2SC(){
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < getNbSlicesLimitFromEdge(e); s++){
            setMaxUsedSliceOverall2SC_k(e,s,getMaxUsedSliceOverall2Multiplier_k(e,s));
        }
    }
}

/* Updates maximum used slice overall 3 stability center */
void AbstractLagFormulation::updateMaxUsedSliceOverall3SC(){
    for (int v = 0; v < instance.getNbNodes(); v++){
        for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
            setMaxUsedSliceOverall3SC_k(v,s,getMaxUsedSliceOverall3Multiplier_k(v,s));
        }
    }
}

/************************************************** SLACK ******************************************************/

/* Updates the slack of length constraints. */
void AbstractLagFormulation::updateLengthSlack(int d, const ListDigraph::Arc & arc){
    double maxReach = getToBeRouted_k(d).getMaxLength();
    lengthSlack[d] -= (getArcLength(arc,d)/maxReach);
}

/* Updates the slack of source target constraints. */
void AbstractLagFormulation::updateSourceTargetSlack(int d, const ListDigraph::Arc & arc){
    int nodeLabel = getNodeLabel((*vecGraph[d]).target(arc),d);
    if(nodeLabel== getToBeRouted_k(d).getTarget()){
        sourceTargetSlack[d][nodeLabel]     -= 1.0;
    }
    nodeLabel = getNodeLabel((*vecGraph[d]).source(arc),d);
    sourceTargetSlack[d][nodeLabel]     -= 1.0;
}

/* Updates the slack of source target constraints. */
void AbstractLagFormulation::updateFlowSlack(int d, const ListDigraph::Arc & arc){
    int nodeIndex = getNodeIndex((*vecGraph[d]).target(arc),d);
    flowSlack[d][nodeIndex] += 1.0;
    nodeIndex = getNodeIndex((*vecGraph[d]).source(arc),d);
    flowSlack[d][nodeIndex] -= 1.0;
}

/* Updates the slack of non-overlap constraints. */
void AbstractLagFormulation::updateOverlapSlack(int d,const ListDigraph::Arc & arc){
    int demandLoad = getToBeRouted_k(d).getLoad();
    int slice = getArcSlice(arc, d); 
    int label = getArcLabel(arc, d);
    for(int s = slice-demandLoad + 1;s<= slice;s++){
        overlapSlack[label][s] -= 1.0;
    }
}

/* Updates the slack of one slice per demand constraints. */
void AbstractLagFormulation::updateOneSlicePerDemandSlack(int d, const ListDigraph::Arc & arc){
    int label = getArcLabel(arc, d);
    oneSlicePerDemandSlack[label][d] -= 1.0;
}

/* Updates the slack of max used slice overall constraints. */
void AbstractLagFormulation::updateMaxUsedSliceOverallSlack(int d, const ListDigraph::Arc & arc){
    int sourceLabel = getToBeRouted_k(d).getSource();
    if(getNodeLabel((*vecGraph[d]).source(arc), d) == sourceLabel){
        maxUsedSliceOverallSlack[d]    -= getArcSlice(arc, d);
    }
}

void AbstractLagFormulation::updateMaxUsedSliceOverallSlack_aux(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        maxUsedSliceOverallSlack[d] += maxUsedSliceOverall;
    }
}


/* Updates the slack of max used slice overall (auxiliary) constraints. */
void AbstractLagFormulation::updateMaxUsedSliceOverallAuxSlack(int d, const ListDigraph::Arc & arc){
    int sourceLabel = getToBeRouted_k(d).getSource();
    if(getNodeLabel((*vecGraph[d]).source(arc), d) == sourceLabel){
        maxUsedSliceOverallAuxSlack[d] += getArcSlice(arc, d);
    }
}

void AbstractLagFormulation::updateMaxUsedSliceOverallAuxSlack_aux(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        maxUsedSliceOverallAuxSlack[d] = maxUsedSliceOverallAuxSlack[d] - maxUsedSliceOverall + (getInstance().getMaxSlice())*(1-varAuxZ[d]);      
    }
}

/* Updates the slack of max used slice overall 2 constraints. */
void AbstractLagFormulation::updateMaxUsedSliceOverall2Slack(int d, const ListDigraph::Arc & arc){
    int demandLoad = getToBeRouted_k(d).getLoad();
    int slice = getArcSlice(arc, d); 
    int label = getArcLabel(arc, d);
    for(int s = slice-demandLoad + 1;s<= slice;s++){
        maxUsedSliceOverallSlack2[label][s] -= getArcSlice(arc, d);
    }
}

/* Updates the slack of max used slice overall 3 constraints. */
void AbstractLagFormulation::updateMaxUsedSliceOverall3Slack(int d, const ListDigraph::Arc & arc){
    int demandLoad = getToBeRouted_k(d).getLoad();
    int slice = getArcSlice(arc, d); 
    int nodeLabel = getNodeLabel((*vecGraph[d]).source(arc),d);
    for(int s = slice-demandLoad + 1;s<= slice;s++){
        maxUsedSliceOverallSlack3[nodeLabel][s] -= getArcSlice(arc, d);
    }
}


/********************************************** SLACK PRIMAL APPROXIMATION *********************************************/

/*Updates the length slack considering the primal approximation. */
void AbstractLagFormulation::updateLengthPrimalSlack(double alpha){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        lengthSlack_v2[d] = alpha*lengthSlack[d] + (1-alpha)*lengthSlack_v2[d];
    }
}

/*Updates the source/target slack considering the primal approximation. */
void AbstractLagFormulation::updateSourceTargetPrimalSlack(double alpha){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){ 
            sourceTargetSlack_v2[d][v] = alpha*sourceTargetSlack[d][v] + (1-alpha)*sourceTargetSlack_v2[d][v];
        }
    }
}

/*Updates the flow slack considering the primal approximation. */
void AbstractLagFormulation::updateFlowPrimalSlack(double alpha){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < countNodes(*vecGraph[d]); v++){
            if(v != getTargetNodeIndex(d) && v != getSourceNodeIndex(d)){
                flowSlack_v2[d][v] = alpha*flowSlack[d][v] + (1-alpha)*flowSlack_v2[d][v];
            }
        }
    }
}

/*Updates the overlap slack considering the primal approximation. */
void AbstractLagFormulation::updateOverlapPrimalSlack(double alpha){
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            overlapSlack_v2[e][s] = alpha*overlapSlack[e][s] + (1-alpha)*overlapSlack_v2[e][s];
        }
    }
}

/*Updates the one slice per demand slack considering the primal variables. */
void AbstractLagFormulation::updateOneSlicePerDemandPrimalSlack(double alpha){
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            oneSlicePerDemandSlack_v2[e][d] = alpha*oneSlicePerDemandSlack[e][d] +(1-alpha)*oneSlicePerDemandSlack_v2[e][d];
        }
    }
}

/*Updates the max used slice overall slack considering the primal variables. */
void AbstractLagFormulation::updateMaxUsedSliceOverallPrimalSlack(double alpha){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        maxUsedSliceOverallSlack_v2[d] = alpha*maxUsedSliceOverallSlack[d] + (1-alpha)*maxUsedSliceOverallSlack_v2[d];
    }
}

/*Updates the max used slice overall (auxiliary) slack considering the primal variables. */
void AbstractLagFormulation::updateMaxUsedSliceOverallAuxPrimalSlack(double alpha){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        maxUsedSliceOverallAuxSlack_v2[d] = alpha*maxUsedSliceOverallAuxSlack[d] + (1-alpha)*maxUsedSliceOverallAuxSlack_v2[d];
    }
}

/*Updates the max used slice overall 2 slack considering the primal variables. */
void AbstractLagFormulation::updateMaxUsedSliceOverall2PrimalSlack(double alpha){
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < getNbSlicesLimitFromEdge(e); s++){
            maxUsedSliceOverallSlack2_v2[e][s] = alpha*maxUsedSliceOverallSlack2[e][s] + (1-alpha)*maxUsedSliceOverallSlack2_v2[e][s];
        }
    }
}

/*Updates the max used slice overall 3 slack considering the primal variables. */
void AbstractLagFormulation::updateMaxUsedSliceOverall3PrimalSlack(double alpha){
    for (int v = 0; v < instance.getNbNodes(); v++){
        for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
            maxUsedSliceOverallSlack3_v2[v][s] = alpha*maxUsedSliceOverallSlack3[v][s] + (1-alpha)*maxUsedSliceOverallSlack3_v2[v][s];
        }
    }
}

/******************************************************* DIRECTION *****************************************************/

/* Updates the slack of length constraints. */
void AbstractLagFormulation::updateLengthDirection(double theta){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double value = getLengthSlack_k(d) - theta*getLengthDirection_k(d);
        setLengthDirection_k(d,value);
    }
}

/* Updates the slack of source/target constraints. */
void AbstractLagFormulation::updateSourceTargetDirection(double theta){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){ 
            double value = getSourceTargetSlack_k(d,v) - theta*getSourceTargetDirection_k(d,v);
            setSourceTargetDirection_k(d,v,value);
        }
    }
}

/* Updates the slack of flow constraints. */
void AbstractLagFormulation::updateFlowDirection(double theta){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < countNodes(*vecGraph[d]); v++){
            if(v != getTargetNodeIndex(d) && v != getSourceNodeIndex(d)){
                double value = getFlowSlack_k(d,v) - theta*getFlowDirection_k(d,v);
                setFlowDirection_k(d,v,value);
            }
        }
    }
}

/* Updates the slack of non overlap constraints. */
void AbstractLagFormulation::updateOverlapDirection(double theta){
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            double value = getOverlapSlack_k(e,s) - theta*getOverlapDirection_k(e,s);
            setOverlapDirection_k(e,s,value);
        }
    }
}

/* Updates the slack of one slice per demand constraints. */
void AbstractLagFormulation::updateOneSlicePerDemandDirection(double theta){
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            double value = getOneSlicePerDemandSlack_k(e,d) - theta*getOneSlicePerDemandDirection_k(e,d);
            setOneSlicePerDemandDirection_k(e,d,value);
        }
    }
}

/* Updates the slack of max used slice overall 1 constraints. */
void AbstractLagFormulation::updateMaxUsedSliceOverallDirection(double theta){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double value = getMaxUsedSliceOverallSlack_k(d) - theta*getMaxUsedSliceOverallDirection_k(d);
        setMaxUsedSliceOverallDirection_k(d,value);
    }
}

/* Updates the slack of max used slice overall 1 (auxiliary) constraints. */
void AbstractLagFormulation::updateMaxUsedSliceOverallAuxDirection(double theta){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double value = getMaxUsedSliceOverallAuxSlack_k(d) - theta*getMaxUsedSliceOverallAuxDirection_k(d);
        setMaxUsedSliceOverallAuxDirection_k(d,value);
    }
}

/* Updates the slack of max used slice overall 2 constraints. */
void AbstractLagFormulation::updateMaxUsedSliceOverall2Direction(double theta){
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < getNbSlicesLimitFromEdge(e); s++){
            double value = getMaxUsedSliceOverall2Slack_k(e,s) - theta*getMaxUsedSliceOverall2Direction_k(e,s);
            setMaxUsedSliceOverall2Direction_k(e,s,value);
        }
    }
}

/* Updates the slack of max used slice overall 3 constraints. */
void AbstractLagFormulation::updateMaxUsedSliceOverall3Direction(double theta){
    for (int v = 0; v < instance.getNbNodes(); v++){
        for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
            double value = getMaxUsedSliceOverall3Slack_k(v,s) - theta*getMaxUsedSliceOverall3Direction_k(v,s);
            setMaxUsedSliceOverall3Direction_k(v,s,value);
        }
    }
}

/**************************************************** ASSIGNMENT MATRIX *****************************************************/

void AbstractLagFormulation::updateAssignment(){
    for(int d = 0; d < getNbDemandsToBeRouted(); d++){
        std::fill(assignmentMatrix_d[d].begin(), assignmentMatrix_d[d].end(), false);
    }
}

/**************************************************** PRIMAL SOLUTION *****************************************************/

void AbstractLagFormulation::updatePrimalSolution(double alpha){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for(int index = 0; index < primal_linear_solution[d].size(); index++){
            primal_linear_solution[d][index] = alpha*(double)assignmentMatrix_d[d][index] + (1.0-alpha)*primal_linear_solution[d][index];
        }
    }
    Input::ObjectiveMetric chosenMetric = getInstance().getInput().getChosenObj_k(0);
    if(chosenMetric == Input::OBJECTIVE_METRIC_8){
        primalMaxUsedSliceOverall = alpha*maxUsedSliceOverall + (1-alpha)*primalMaxUsedSliceOverall;
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            primalVarAuxZ[d] = alpha*varAuxZ[d] + (1-alpha)*primalVarAuxZ[d];
        }
    }
}

/************************************************** PRIMAL APPROXIMATION **********************************************/

void AbstractLagFormulation::updatePrimalApproximation(double alpha){
    updatePrimalSolution(alpha);
    updatePrimalSlack(alpha);

    double obj = alpha*getLagrCurrentCost() + (1.0-alpha)*getPrimalCurrentCost();
    setCurrentPrimalCost(obj);
}

/* ********************************************************************************************************************
*                                                   Check Feasibility METHODS
********************************************************************************************************************** */

/* Checks if all length slacks are non-negative. */
bool AbstractLagFormulation::checkLengthFeasibility(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        if (lengthSlack[d] < -DBL_EPSILON){
            return false;
        }
    }
    return true;
}

/* Checks if all source/target constraints are feasible. */
/*  For the source and destination : its an equality, it must be 0, for the rest it must be non negative */
bool AbstractLagFormulation::checkSourceTargetFeasibility(){
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

/* Checks if all flow constraints are feasible. */
/* Equalities: it most be 0.*/
bool AbstractLagFormulation::checkFlowFeasibility(){
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

/* Checks if all overlap slacks are non-negative. */
bool AbstractLagFormulation::checkOverlapFeasibility(){   
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            if (overlapSlack[e][s] < -DBL_EPSILON){
                return false;
            }
        }
    }
    return true;
}

/* Checks if all max used slice overall slacks are non-negative. */
bool AbstractLagFormulation::checkMaxUsedSliceOverallFeasibility(){
     for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        if (maxUsedSliceOverallSlack[d] < -DBL_EPSILON){
            return false;
        }
    }
    return true;
}

/* Checks if all length slacks (considering primal approximation) are non-negative. */
bool AbstractLagFormulation::checkLengthFeasibility_v2(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        if (lengthSlack_v2[d] < -DBL_EPSILON){
            return false;
        }
    }
    return true;
}

/* Checks if all source/target constraints are feasible. (slack considering primal approximation) */
/*  For the source and destination : its an equality, it must be 0, for the rest it must be non negative */
bool AbstractLagFormulation::checkSourceTargetFeasibility_v2(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){ 
            if((v == getToBeRouted_k(d).getSource()) ||  (v == getToBeRouted_k(d).getTarget())){
                if((sourceTargetSlack_v2[d][v] < -DBL_EPSILON) || (sourceTargetSlack_v2[d][v] > DBL_EPSILON)){
                    return false;
                }
            }else{
                if(sourceTargetSlack_v2[d][v]< -DBL_EPSILON){
                    return false;
                }
            }
        }
    }
    return true;
}

/* Checks if all flow constraints are feasible. (slack considering primal approximation)*/
/* Equalities: it most be 0.*/
bool AbstractLagFormulation::checkFlowFeasibility_v2(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){ 
            if((v != getToBeRouted_k(d).getSource()) &&  (v != getToBeRouted_k(d).getTarget())){
                if((flowSlack_v2[d][v] < -DBL_EPSILON) || (flowSlack_v2[d][v] > DBL_EPSILON)){
                    return false;
                }
            }
        }
    }
    return true;
}

/* Checks if all overlap slacks (considering primal approximation) are non-negative. */
bool AbstractLagFormulation::checkOverlapFeasibility_v2(){   
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            if (overlapSlack_v2[e][s] < -DBL_EPSILON){
                return false;
            }
        }
    }
    return true;
}

/* Checks if all max used slice overall slacks (considering primal approximation) are non-negative. */
bool AbstractLagFormulation::checkMaxUsedSliceOverallFeasibility_v2(){
     for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        if (maxUsedSliceOverallSlack_v2[d] < -DBL_EPSILON){
            return false;
        }
    }
    return true;
}

/* ********************************************************************************************************************
*                                            Slackness conditions methods
********************************************************************************************************************** */

/* Checks the slackness condition of length constraints. */
bool AbstractLagFormulation::checkLengthSlacknessCondition(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        if(!(lagrangianMultiplierLength[d] < DBL_EPSILON || (-lengthSlack[d] < DBL_EPSILON && -lengthSlack[d] > -DBL_EPSILON))){
            return false;
        }
    }
    return true;
}

/* Checks the slackness condition of source/target constraints. */
bool AbstractLagFormulation::checkSourceTargetSlacknessCondition(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){ 
            if((v == getToBeRouted_k(d).getSource()) ||  (v == getToBeRouted_k(d).getTarget())){
                if(!((lagrangianMultiplierSourceTarget[d][v] < DBL_EPSILON && lagrangianMultiplierSourceTarget[d][v] > -DBL_EPSILON) || (-sourceTargetSlack[d][v] < DBL_EPSILON && -sourceTargetSlack[d][v] > -DBL_EPSILON))){
                    return false;
                }
            }else{
                if(!(lagrangianMultiplierSourceTarget[d][v] < DBL_EPSILON || (-sourceTargetSlack[d][v] < DBL_EPSILON && -sourceTargetSlack[d][v] > -DBL_EPSILON))){
                    return false;
                }
            }
        }
    }
    return true;
}

/* Checks the slackness condition of flow constraints */
bool AbstractLagFormulation::checkFlowSlacknessCondition(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){ 
            if((v != getToBeRouted_k(d).getSource()) &&  (v != getToBeRouted_k(d).getTarget())){
                if(!((lagrangianMultiplierFlow[d][v] < DBL_EPSILON && lagrangianMultiplierFlow[d][v] > -DBL_EPSILON) || (-flowSlack[d][v] < DBL_EPSILON && -flowSlack[d][v] > -DBL_EPSILON))){
                    return false;
                }
            }
        }
    }
    return true;
}

/* Checks the slackness condition of overlap constraints. */
bool AbstractLagFormulation::checkOverlapSlacknessCondition(){   
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            if(!(lagrangianMultiplierOverlap[e][s] < DBL_EPSILON || (-overlapSlack[e][s] < DBL_EPSILON && -overlapSlack[e][s] > -DBL_EPSILON))){
                return false;
            }
        }
    }
    return true;
}

/* Checks the slackness condition of max used slice overall constraints. */
bool AbstractLagFormulation::checkMaxUsedSliceOverallSlacknessCondition(){
     for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        if(!(lagrangianMultiplierMaxUsedSliceOverall[d] < DBL_EPSILON || (-maxUsedSliceOverallAuxSlack[d] < DBL_EPSILON && -maxUsedSliceOverallAuxSlack[d] > -DBL_EPSILON))){
            return false;
        }
    }
    return true;
}

