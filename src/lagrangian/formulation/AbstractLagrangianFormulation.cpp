#include "AbstractLagrangianFormulation.h"

AbstractLagFormulation::AbstractLagFormulation(const Instance &instance): RSA(instance), time(ClockTime::getTimeNow()){}

void AbstractLagFormulation::updatePrimalSolution(double alpha){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            int index = getArcIndex(a,d);
            primal_linear_solution[d][index] = alpha*(double)assignmentMatrix_d[d][index] + (1-alpha)*primal_linear_solution[d][index];
        }
    }

}

/* Initializes the primal solution */
void AbstractLagFormulation::initPrimalSolution(){

    primal_linear_solution.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        primal_linear_solution[d].resize(countArcs(*vecGraph[d]));
        std::fill(primal_linear_solution[d].begin(), primal_linear_solution[d].end(), 0.0);
    } 

    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            int index = getArcIndex(a,d);
            primal_linear_solution[d][index] = (double)assignmentMatrix_d[d][index];
        }
    }
}

double AbstractLagFormulation::getPrimalObjective(){
    double total = 0.0;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            int index = getArcIndex(a,d);
            total += getCoeff(a, d)*primal_linear_solution[d][index];
        }
    }
    return total;
}

/* *******************************************************************************
*                             INITIALIZATION METHODS
******************************************************************************* */

/********************************* MULTIPLIERS ***********************************/

/* Sets the initial lagrangian multipliers associated with length constraints. */
void AbstractLagFormulation::initializeLengthMultipliers(double initialMultiplier){  
    lagrangianMultiplierLength.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        lagrangianMultiplierLength[d] = initialMultiplier;
    }
}

/* Sets the initial lagrangian multipliers associated with source/target constraints*/
void AbstractLagFormulation::initializeSourceTargetMultipliers(double initialMultiplier){
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
void AbstractLagFormulation::initializeFlowMultipliers(double initialMultiplier){
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

/* Sets the initial lagrangian multipliers associated with non-overlapping constraints. */
void AbstractLagFormulation::initializeOverlapMultipliers(double initialMultiplier){
    lagrangianMultiplierOverlap.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        lagrangianMultiplierOverlap[e].resize(instance.getPhysicalLinkFromIndex(e).getNbSlices());
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            lagrangianMultiplierOverlap[e][s] = initialMultiplier;
        }
    }
}

void AbstractLagFormulation::initializeMaxUsedSliceOverallMultipliers(double initialMultiplier){
    lagrangianMultiplierMaxUsedSliceOverall.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        lagrangianMultiplierMaxUsedSliceOverall[d] = initialMultiplier;
    }
}

void AbstractLagFormulation::initializeMaxUsedSliceOverall2Multipliers(double initialMultiplier){
    lagrangianMultiplierMaxUsedSliceOverall2.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        lagrangianMultiplierMaxUsedSliceOverall2[e].resize(getNbSlicesLimitFromEdge(e));
        for (int s = 0; s < getNbSlicesLimitFromEdge(e); s++){
            lagrangianMultiplierMaxUsedSliceOverall2[e][s] = initialMultiplier;
        }
    }
}

void AbstractLagFormulation::initializeMaxUsedSliceOverall3Multipliers(double initialMultiplier){
    lagrangianMultiplierMaxUsedSliceOverall3.resize(instance.getNbNodes());
    for (int v = 0; v < instance.getNbNodes(); v++){
        lagrangianMultiplierMaxUsedSliceOverall3[v].resize(getNbSlicesGlobalLimit());
        for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
            lagrangianMultiplierMaxUsedSliceOverall3[v][s] = initialMultiplier;
        }
    }
}

/********************************* STABILITY CENTER ***********************************/

/** Sets the initial lagrangian stability center associated with length constraints. **/
void AbstractLagFormulation::initializeLengthSC(){
    lagrangianSCLength.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        lagrangianSCLength[d] = lagrangianMultiplierLength[d];
    }
}

/** Sets the initial lagrangian stability center associated with Source/Target constraints **/
void AbstractLagFormulation::initializeSourceTargetSC(){
    lagrangianSCSourceTarget.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        lagrangianSCSourceTarget[d].resize(instance.getNbNodes());
        for (int v = 0; v < instance.getNbNodes(); v++){
            if((v == getToBeRouted_k(d).getSource()) || v == getToBeRouted_k(d).getTarget()) {
                //int load = getToBeRouted_k(d).getLoad();
                //lagrangianSCSourceTarget[d][v]= -load/2;
                lagrangianSCSourceTarget[d][v]= lagrangianMultiplierSourceTarget[d][v];
            }
            else{
                lagrangianMultiplierSourceTarget[d][v]= lagrangianMultiplierSourceTarget[d][v];
            }
        }
    }
}
        
/** Sets the initial lagrangian stability center associated with flow constraints **/
void AbstractLagFormulation::initializeFlowSC(){
    lagrangianSCFlow.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        lagrangianSCFlow[d].resize(instance.getNbNodes());
        for (int v = 0; v < instance.getNbNodes(); v++){
            /* The multiplier is not defined for the source and the target */
            if((v == getToBeRouted_k(d).getSource()) || v == getToBeRouted_k(d).getTarget()) {
                lagrangianSCFlow[d][v]= 0.0;
            }else{
                lagrangianSCFlow[d][v]= lagrangianMultiplierFlow[d][v];
            }
        }
    }

}

/* Sets the initial lagrangian stability center associated with non-overlapping constraints. */
void AbstractLagFormulation::initializeOverlapSC(){
    lagrangianSCOverlap.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        lagrangianSCOverlap[e].resize(instance.getPhysicalLinkFromIndex(e).getNbSlices());
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            lagrangianSCOverlap[e][s] = lagrangianMultiplierOverlap[e][s];
        }
    }
}

void AbstractLagFormulation::initializeMaxUsedSliceOverallSC(){
    lagrangianSCMaxUsedSliceOverall.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        lagrangianSCMaxUsedSliceOverall[d] = lagrangianMultiplierMaxUsedSliceOverall[d];
    }
}

void AbstractLagFormulation::initializeMaxUsedSliceOverall2SC(){
    lagrangianSCMaxUsedSliceOverall2.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        lagrangianSCMaxUsedSliceOverall2[e].resize(getNbSlicesLimitFromEdge(e));
        for (int s = 0; s < getNbSlicesLimitFromEdge(e); s++){
            lagrangianSCMaxUsedSliceOverall2[e][s] = lagrangianMultiplierMaxUsedSliceOverall2[e][s];
        }
    }
}

void AbstractLagFormulation::initializeMaxUsedSliceOverall3SC(){
    lagrangianSCMaxUsedSliceOverall3.resize(instance.getNbNodes());
    for (int v = 0; v < instance.getNbNodes(); v++){
        lagrangianSCMaxUsedSliceOverall3[v].resize(getNbSlicesGlobalLimit());
        for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
            lagrangianSCMaxUsedSliceOverall3[v][s] = lagrangianMultiplierMaxUsedSliceOverall3[v][s];
        }
    }
}

/********************************** SLACK ***************************************/

/* Initializes the slack of Length constraints. */ 
void AbstractLagFormulation::initializeLengthSlacks(){
    lengthSlack.resize(getNbDemandsToBeRouted(), 0.0);
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        //lengthSlack[d] = getToBeRouted_k(d).getMaxLength();
        //lengthSlack[d] = getToBeRouted_k(d).getMaxLength()/100;
        lengthSlack[d] = 1;
    }
}

/* Initializes the slack of Source/Target constraints. */
void AbstractLagFormulation::initializeSourceTargetSlacks(){
    sourceTargetSlack.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        sourceTargetSlack[d].resize(instance.getNbNodes());
        for(int v = 0; v < instance.getNbNodes();v++){
            sourceTargetSlack[d][v] = 1;
        }
    }
}

/* Initializes the slack of Flow constraints. */
void AbstractLagFormulation::initializeFlowSlacks(){
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

/* Initializes the slack of non-overlap constraints. */
void AbstractLagFormulation::initializeOverlapSlacks(){
    overlapSlack.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        overlapSlack[e].resize(instance.getPhysicalLinkFromIndex(e).getNbSlices());
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            overlapSlack[e][s] = 1;
        }
    }
}

void AbstractLagFormulation::initializeMaxUsedSliceOverallSlacks(){
    maxUsedSliceOverallSlack.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        maxUsedSliceOverallSlack[d] = 0.0;
    }
}

void AbstractLagFormulation::initializeMaxUsedSliceOverall2Slacks(){
    maxUsedSliceOverallSlack2.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        maxUsedSliceOverallSlack2[e].resize(getNbSlicesLimitFromEdge(e));
        for (int s = 0; s < getNbSlicesLimitFromEdge(e); s++){
            maxUsedSliceOverallSlack2[e][s] = 0.0;
        }
    }
}

void AbstractLagFormulation::initializeMaxUsedSliceOverall3Slacks(){
    maxUsedSliceOverallSlack3.resize(instance.getNbNodes());
    for (int v = 0; v < instance.getNbNodes(); v++){
        maxUsedSliceOverallSlack3[v].resize(getNbSlicesGlobalLimit());
        for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
            maxUsedSliceOverallSlack3[v][s] = 0.0;
        }
    }
}

/********************************** SLACK CONSIDERING PRIMAL VARIABLES ***************************************/

void AbstractLagFormulation::initializeLengthSlacks_v2(){
    lengthSlack_v2.resize(getNbDemandsToBeRouted(), 0.0);
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        //lengthSlack[d] = getToBeRouted_k(d).getMaxLength();
        //lengthSlack[d] = getToBeRouted_k(d).getMaxLength()/100;
        lengthSlack_v2[d] = 1;
    }
}
        
void AbstractLagFormulation::initializeSourceTargetSlacks_v2(){
    sourceTargetSlack_v2.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        sourceTargetSlack_v2[d].resize(instance.getNbNodes());
        for(int v = 0; v < instance.getNbNodes();v++){
            sourceTargetSlack_v2[d][v] = 1;
        }
    }
}

void AbstractLagFormulation::initializeFlowSlacks_v2(){
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

void AbstractLagFormulation::initializeOverlapSlacks_v2(){
    overlapSlack_v2.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        overlapSlack_v2[e].resize(instance.getPhysicalLinkFromIndex(e).getNbSlices());
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            overlapSlack_v2[e][s] = 1;
        }
    }
}

void AbstractLagFormulation::initializeMaxUsedSliceOverallSlacks_v2(){
    maxUsedSliceOverallSlack_v2.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        maxUsedSliceOverallSlack_v2[d] = 0.0;
    }
}

void AbstractLagFormulation::initializeMaxUsedSliceOverall2Slacks_v2(){
    maxUsedSliceOverallSlack2_v2.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        maxUsedSliceOverallSlack2_v2[e].resize(getNbSlicesLimitFromEdge(e));
        for (int s = 0; s < getNbSlicesLimitFromEdge(e); s++){
            maxUsedSliceOverallSlack2_v2[e][s] = 0.0;
        }
    }
}

void AbstractLagFormulation::initializeMaxUsedSliceOverall3Slacks_v2(){
    maxUsedSliceOverallSlack3_v2.resize(instance.getNbNodes());
    for (int v = 0; v < instance.getNbNodes(); v++){
        maxUsedSliceOverallSlack3_v2[v].resize(getNbSlicesGlobalLimit());
        for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
            maxUsedSliceOverallSlack3_v2[v][s] = 0.0;
        }
    }
}

/******************************** DIRECTION *********************************/

void AbstractLagFormulation::initializeLengthDirection(){
    lengthDirection.resize(getNbDemandsToBeRouted(), 0.0);
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        lengthDirection[d] = lengthSlack[d];
    }
}

void AbstractLagFormulation::initializeSourceTargetDirection(){
    sourceTargetDirection.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        sourceTargetDirection[d].resize(instance.getNbNodes());
        for(int v = 0; v < instance.getNbNodes();v++){
            sourceTargetDirection[d][v] = sourceTargetSlack[d][v];
        }
    }
}

void AbstractLagFormulation::initializeFlowDirection(){
    flowDirection.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        flowDirection[d].resize(instance.getNbNodes());
        for(int v = 0; v < instance.getNbNodes();v++){
            flowDirection[d][v] = flowSlack[d][v];
        }
    }
}

void AbstractLagFormulation::initializeOverlapDirection(){
    overlapDirection.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        overlapDirection[e].resize(instance.getPhysicalLinkFromIndex(e).getNbSlices());
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            overlapDirection[e][s] = overlapSlack[e][s];
        }
    }
}

void AbstractLagFormulation::initializeMaxUsedSliceOverallDirection(){
    maxUsedSliceOverallDirection.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        maxUsedSliceOverallDirection[d] = maxUsedSliceOverallSlack[d];
    }
}

void AbstractLagFormulation::initializeMaxUsedSliceOverall2Direction(){
    maxUsedSliceOverall2Direction.resize(instance.getNbEdges());
    for (int e = 0; e < instance.getNbEdges(); e++){
        maxUsedSliceOverall2Direction[e].resize(getNbSlicesLimitFromEdge(e));
        for (int s = 0; s < getNbSlicesLimitFromEdge(e); s++){
            maxUsedSliceOverall2Direction[e][s] = maxUsedSliceOverallSlack2[e][s];
        }
    }
}

void AbstractLagFormulation::initializeMaxUsedSliceOverall3Direction(){
    maxUsedSliceOverall3Direction.resize(instance.getNbNodes());
    for (int v = 0; v < instance.getNbNodes(); v++){
        maxUsedSliceOverall3Direction[v].resize(getNbSlicesGlobalLimit());
        for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
            maxUsedSliceOverall3Direction[v][s] = maxUsedSliceOverallSlack3[v][s];
        }
    }
}

/********************************** ASSIGNMENT MATRIX ***************************************/

void AbstractLagFormulation::initAssignmentMatrix(){
    assignmentMatrix_d.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        assignmentMatrix_d[d].resize(countArcs(*vecGraph[d]));
        std::fill(assignmentMatrix_d[d].begin(), assignmentMatrix_d[d].end(), false);
    }
    maxUsedSliceOverall = 0.0;
}

/* *******************************************************************************
*                             UPDATE METHODS
******************************************************************************* */

/********************************** MULTIPLIERS ***************************************/

/** Update length multipliers **/
void AbstractLagFormulation::updateLengthMultiplier(double step){
    /* length */
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double violation = -getLengthSlack_k(d); /* the constraint is <=, we have to pass to >= as it is a minimization problem*/
        double new_multipliplier = getLengthMultiplier_k(d) + (step*violation);
        setLengthMultiplier_k(d, std::max(new_multipliplier, 0.0)); /* the multiplier is >=0 */
    }
}

/** update source target multipliers **/
void AbstractLagFormulation::updateSourceTargetMultiplier(double step){
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

void AbstractLagFormulation::updateFlowMultiplier(double step){
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

void AbstractLagFormulation::updateOverlapMultiplier(double step){
    //overlap 
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            double violation = -getOverlapSlack_k( e, s);
            double new_multipliplier = getOverlapMultiplier_k( e, s) + (step*violation);
            setOverlapMultiplier_k( e, s, std::max(new_multipliplier, 0.0));
        }
    }  
}

void AbstractLagFormulation::updateMaxUsedSliceOverallMultiplier(double step){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double violation = -getMaxUsedSliceOverallSlack_k(d);
        double new_multipliplier = getMaxUsedSliceOverallMultiplier_k(d) + (step*violation);
        setMaxUsedSliceOverallMultiplier_k(d,std::max(new_multipliplier, 0.0));
    }
}

void AbstractLagFormulation::updateMaxUsedSliceOverall2Multiplier(double step){
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < getNbSlicesLimitFromEdge(e); s++){
            double violation = -getMaxUsedSliceOverall2Slack_k(e,s);
            double new_multipliplier = getMaxUsedSliceOverall2Multiplier_k(e,s) + (step*violation);
            setMaxUsedSliceOverall2Multiplier_k(e,s,std::max(new_multipliplier, 0.0)); 
        }
    }
}

void AbstractLagFormulation::updateMaxUsedSliceOverall3Multiplier(double step){
    for (int v = 0; v < instance.getNbNodes(); v++){
        for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
            double violation = -getMaxUsedSliceOverall3Slack_k(v,s);
            double new_multipliplier = getMaxUsedSliceOverall3Multiplier_k(v,s) + (step*violation);
            setMaxUsedSliceOverall3Multiplier_k(v,s,std::max(new_multipliplier, 0.0)); 
        }
    }
}

/****************** MULTIPLIER CONSIDERING THE STABILITY CENTER *********************/

void AbstractLagFormulation::updateLengthMultiplier_v2(double step){
    /* length */
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double violation = -getLengthSlack_v2_k(d); /* the constraint is <=, we have to pass to >= as it is a minimization problem*/
        double new_multipliplier = getLengthSC_k(d) + (step*violation);
        setLengthMultiplier_k(d, std::max(new_multipliplier, 0.0)); /* the multiplier is >=0 */
    }
}

void AbstractLagFormulation::updateSourceTargetMultiplier_v2(double step){
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

void AbstractLagFormulation::updateFlowMultiplier_v2(double step){
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

void AbstractLagFormulation::updateOverlapMultiplier_v2(double step){

    //overlap 
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            double violation = -getOverlapSlack_v2_k( e, s);
            double new_multipliplier = getOverlapSC_k( e, s) + (step*violation);
            setOverlapMultiplier_k( e, s, std::max(new_multipliplier, 0.0));
        }
    }  
}

void AbstractLagFormulation::updateMaxUsedSliceOverallMultiplier_v2(double step){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double violation = -getMaxUsedSliceOverallSlack_v2_k(d);
        double new_multipliplier = getMaxUsedSliceOverallSC_k(d) + (step*violation);
        setMaxUsedSliceOverallMultiplier_k(d,std::max(new_multipliplier, 0.0));
    }
}

void AbstractLagFormulation::updateMaxUsedSliceOverall2Multiplier_v2(double step){
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < getNbSlicesLimitFromEdge(e); s++){
            double violation = -getMaxUsedSliceOverall2Slack_v2_k(e,s);
            double new_multipliplier = getMaxUsedSliceOverall2SC_k(e,s) + (step*violation);
            setMaxUsedSliceOverall2Multiplier_k(e,s,std::max(new_multipliplier, 0.0)); 
        }
    }
}

void AbstractLagFormulation::updateMaxUsedSliceOverall3Multiplier_v2(double step){
    for (int v = 0; v < instance.getNbNodes(); v++){
        for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
            double violation = -getMaxUsedSliceOverall3Slack_v2_k(v,s);
            double new_multipliplier = getMaxUsedSliceOverall3SC_k(v,s) + (step*violation);
            setMaxUsedSliceOverall3Multiplier_k(v,s,std::max(new_multipliplier, 0.0)); 
        }
    }
}

/********************************** STABILITY CENTER ***************************************/

void AbstractLagFormulation::updateLengthSC(){
    /* length */
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        setLengthSC_k(d, getLengthMultiplier_k(d)); 
    }
}

void AbstractLagFormulation::updateSourceTargetSC(){
     /* source target */
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){ 
            setSourceTargetSC_k(d,v,getSourceTargetMultiplier_k(d,v)); 
        }
    }

}

void AbstractLagFormulation::updateFlowSC(){
    /* flow */
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){
            setFlowSC_k(d,v, getFlowMultiplier_k(d,v)); 
        }
    }
}

void AbstractLagFormulation::updateOverlapSC(){
    //overlap 
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            setOverlapSC_k( e, s, getOverlapMultiplier_k( e, s));
        }
    }  
}

void AbstractLagFormulation::updateMaxUsedSliceOverallSC(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        setMaxUsedSliceOverallSC_k(d,getMaxUsedSliceOverallMultiplier_k(d));
    }
}

void AbstractLagFormulation::updateMaxUsedSliceOverall2SC(){
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < getNbSlicesLimitFromEdge(e); s++){
            setMaxUsedSliceOverall2SC_k(e,s,getMaxUsedSliceOverall2Multiplier_k(e,s));
        }
    }
}

void AbstractLagFormulation::updateMaxUsedSliceOverall3SC(){
    for (int v = 0; v < instance.getNbNodes(); v++){
        for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
            setMaxUsedSliceOverall3SC_k(v,s,getMaxUsedSliceOverall3Multiplier_k(v,s));
        }
    }
}

/********************************** SLACK ***************************************/

/* Updates the slack of length constraints. */
void AbstractLagFormulation::updateLengthSlack(){
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

/* Updates the slack of source target constraints. */
void AbstractLagFormulation::updateSourceTargetSlack(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int label = 0; label < instance.getNbNodes(); label++){ 
            double exp = 0.0;
            for (ListDigraph::NodeIt v(*vecGraph[d]); v != INVALID; ++v){
                if (getNodeLabel(v, d) == label){
                    if(label == getToBeRouted_k(d).getTarget()){
                        for (ListDigraph::InArcIt a((*vecGraph[d]), v); a != INVALID; ++a){
                            int indexArc = getArcIndex(a,d);
                            exp += assignmentMatrix_d[d][indexArc];
                        }
                    }else{
                        for (ListDigraph::OutArcIt a((*vecGraph[d]), v); a != INVALID; ++a){
                            int indexArc = getArcIndex(a,d);
                            exp += assignmentMatrix_d[d][indexArc];
                        }
                    }
                }
            }
            sourceTargetSlack[d][label] = 1 - exp;
        }
    }
    //std::cout << "\t> Source/Target slack was updated. " << std::endl;
}

/* Updates the slack of source target constraints. */
void AbstractLagFormulation::updateFlowSlack(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::NodeIt v(*vecGraph[d]); v != INVALID; ++v){
            int label = getNodeLabel(v, d);
            /* The multiplier is not defined for the source and the target */
            if( (label != getToBeRouted_k(d).getSource()) && (label != getToBeRouted_k(d).getTarget()) ){
                double exp = 0.0;
                for (ListDigraph::OutArcIt a((*vecGraph[d]), v); a != INVALID; ++a){
                    int indexArc = getArcIndex(a,d);
                    exp -= assignmentMatrix_d[d][indexArc];
                }
                for (ListDigraph::InArcIt a((*vecGraph[d]), v); a != INVALID; ++a){
                    int indexArc = getArcIndex(a,d);
                    exp += assignmentMatrix_d[d][indexArc];
                }
                flowSlack[d][label] = 0 + exp;
            }
        }
    } 
    //std::cout << "\t> Flow slack was updated. " << std::endl; 
}

/* Updates the slack of non-overlap constraints. */
void AbstractLagFormulation::updateOverlapSlack(){
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
    //std::cout << "\t> Overlap slack was updated. " << std::endl;
}

/* Updates the slack of max used slice overall constraints. */
void AbstractLagFormulation::updateMaxUsedSliceOverallSlack(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double exp = 0.0;
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            if (getToBeRouted_k(d).getSource() == getNodeLabel((*vecGraph[d]).source(a), d)){
                int indexArc = getArcIndex(a, d);
                exp += assignmentMatrix_d[d][indexArc]*getArcSlice(a, d);
                //exp += assignmentMatrix_d[d][indexArc]*getArcSlice(a, d)/getInstance().getMaxSlice();
            }
        }
        maxUsedSliceOverallSlack[d] = 0 + maxUsedSliceOverall - exp;
    }
}

/* Updates the slack of max used slice overall 2 constraints. */
void AbstractLagFormulation::updateMaxUsedSliceOverall2Slack(){
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < getNbSlicesLimitFromEdge(e); s++){
            double exp = 0.0;
            for (int d = 0; d < getNbDemandsToBeRouted(); d++){
                int demandLoad = getToBeRouted_k(d).getLoad();
                for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
                    if ((getArcLabel(a, d) == e) && (getArcSlice(a, d) >= s) && (getArcSlice(a, d) <= s + demandLoad - 1)){
                        int indexArc = getArcIndex(a, d);
                        exp += assignmentMatrix_d[d][indexArc]*getArcSlice(a, d);
                        //exp += assignmentMatrix_d[d][indexArc]*getArcSlice(a, d)/getInstance().getMaxSlice();
                    }
                }
            }
            maxUsedSliceOverallSlack2[e][s] = 0 + maxUsedSliceOverall - exp;
        }
    }
}

/* Updates the slack of max used slice overall 3 constraints. */
void AbstractLagFormulation::updateMaxUsedSliceOverall3Slack(){
    for (int node = 0; node < instance.getNbNodes(); node++){
        for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
            double exp = 0.0;
            int degree = 0;
            for (ListGraph::NodeIt v(compactGraph); v != INVALID; ++v){
                if ((getCompactNodeLabel(v) == node)){
                    for (ListGraph::IncEdgeIt a(compactGraph, v); a != INVALID; ++a){
                        degree++;
                    }
                }
            }

            for (int d = 0; d < getNbDemandsToBeRouted(); d++){
                int demandLoad = getToBeRouted_k(d).getLoad();
                for (ListDigraph::NodeIt v(*vecGraph[d]); v != INVALID; ++v){
                    if ((getNodeLabel(v, d) == node)){
                        for (ListDigraph::OutArcIt a(*vecGraph[d], v); a != INVALID; ++a){
                            if ( (getArcSlice(a, d) >= s) && (getArcSlice(a, d) <= s + demandLoad - 1) ){
                                int indexArc = getArcIndex(a, d);
                                exp += assignmentMatrix_d[d][indexArc]*getArcSlice(a, d);
                                //exp += assignmentMatrix_d[d][indexArc]*getArcSlice(a, d)/getInstance().getMaxSlice();
                            }
                        }
                    }
                }
            }
            maxUsedSliceOverallSlack3[node][s] = 0 + maxUsedSliceOverall*degree - exp;
        }
    }
}

/************************ SLACK CONSIDERING THE PRIMAL VECTOR ********************/

/* Updates the slack (primal variables) of length constraints. */
void AbstractLagFormulation::updateLengthSlack_v2(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double exp = 0.0;

        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            int index = getArcIndex(a,d);
            //exp += assignmentMatrix_d[d][index]*getArcLength(a,d);
            exp += (primal_linear_solution[d][index]*getArcLength(a,d))/(getToBeRouted_k(d).getMaxLength());
        }
        //const int MAX_LENGTH = getToBeRouted_k(d).getMaxLength(); 
        const int MAX_LENGTH = 1;
        lengthSlack_v2[d] = (MAX_LENGTH - exp);
    }
}

/* Updates the slack (primal variables) of source target constraints. */
void AbstractLagFormulation::updateSourceTargetSlack_v2(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int label = 0; label < instance.getNbNodes(); label++){ 
            double exp = 0.0;
            for (ListDigraph::NodeIt v(*vecGraph[d]); v != INVALID; ++v){
                if (getNodeLabel(v, d) == label){
                    if(label == getToBeRouted_k(d).getTarget()){
                        for (ListDigraph::InArcIt a((*vecGraph[d]), v); a != INVALID; ++a){
                            int indexArc = getArcIndex(a,d);
                            exp += primal_linear_solution[d][indexArc];
                        }
                    }else{
                        for (ListDigraph::OutArcIt a((*vecGraph[d]), v); a != INVALID; ++a){
                            int indexArc = getArcIndex(a,d);
                            exp += primal_linear_solution[d][indexArc];
                        }
                    }
                }
            }
            sourceTargetSlack_v2[d][label] = 1 - exp;
        }
    }
    //std::cout << "\t> Source/Target slack was updated. " << std::endl;
}

/* Updates the slack (primal variables) of flow constraints. */
void AbstractLagFormulation::updateFlowSlack_v2(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::NodeIt v(*vecGraph[d]); v != INVALID; ++v){
            int label = getNodeLabel(v, d);
            /* The multiplier is not defined for the source and the target */
            if( (label != getToBeRouted_k(d).getSource()) && (label != getToBeRouted_k(d).getTarget()) ){
                double exp = 0.0;
                for (ListDigraph::OutArcIt a((*vecGraph[d]), v); a != INVALID; ++a){
                    int indexArc = getArcIndex(a,d);
                    exp -= primal_linear_solution[d][indexArc];
                }
                for (ListDigraph::InArcIt a((*vecGraph[d]), v); a != INVALID; ++a){
                    int indexArc = getArcIndex(a,d);
                    exp += primal_linear_solution[d][indexArc];
                }
                flowSlack_v2[d][label] = 0 + exp;
            }
        }
    } 
    //std::cout << "\t> Flow slack was updated. " << std::endl; 
}

/* Updates the slack (primal variables) of overlap constraints. */
void AbstractLagFormulation::updateOverlapSlack_v2(){
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


/* Updates the slack (primal variables) of max used slice overall constraints. */
void AbstractLagFormulation::updateMaxUsedSliceOverallSlack_v2(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double exp = 0.0;
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            if (getToBeRouted_k(d).getSource() == getNodeLabel((*vecGraph[d]).source(a), d)){
                int indexArc = getArcIndex(a, d);
                exp += primal_linear_solution[d][indexArc]*getArcSlice(a, d);
                //exp += primal_linear_solution[d][indexArc]*getArcSlice(a, d)/getInstance().getMaxSlice();
            }
        }
        maxUsedSliceOverallSlack_v2[d] = 0 + maxUsedSliceOverall - exp;
    }
}

/* Updates the slack (primal variables) of max used slice overall 2 constraints. */
void AbstractLagFormulation::updateMaxUsedSliceOverall2Slack_v2(){
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < getNbSlicesLimitFromEdge(e); s++){
            double exp = 0.0;
            for (int d = 0; d < getNbDemandsToBeRouted(); d++){
                int demandLoad = getToBeRouted_k(d).getLoad();
                for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
                    if ((getArcLabel(a, d) == e) && (getArcSlice(a, d) >= s) && (getArcSlice(a, d) <= s + demandLoad - 1)){
                        int indexArc = getArcIndex(a, d);
                        exp += primal_linear_solution[d][indexArc]*getArcSlice(a, d);
                        //exp += primal_linear_solution[d][indexArc]*getArcSlice(a, d)/getInstance().getMaxSlice();
                    }
                }
            }
            maxUsedSliceOverallSlack2_v2[e][s] = 0 + maxUsedSliceOverall - exp;
        }
    }
}

/* Updates the slack (primal variables) of max used slice overall 3 constraints. */
void AbstractLagFormulation::updateMaxUsedSliceOverall3Slack_v2(){
    for (int node = 0; node < instance.getNbNodes(); node++){
        for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
            double exp = 0.0;
            int degree = 0;
            for (ListGraph::NodeIt v(compactGraph); v != INVALID; ++v){
                if ((getCompactNodeLabel(v) == node)){
                    for (ListGraph::IncEdgeIt a(compactGraph, v); a != INVALID; ++a){
                        degree++;
                    }
                }
            }

            for (int d = 0; d < getNbDemandsToBeRouted(); d++){
                int demandLoad = getToBeRouted_k(d).getLoad();
                for (ListDigraph::NodeIt v(*vecGraph[d]); v != INVALID; ++v){
                    if ((getNodeLabel(v, d) == node)){
                        for (ListDigraph::OutArcIt a(*vecGraph[d], v); a != INVALID; ++a){
                            if ( (getArcSlice(a, d) >= s) && (getArcSlice(a, d) <= s + demandLoad - 1) ){
                                int indexArc = getArcIndex(a, d);
                                exp += primal_linear_solution[d][indexArc]*getArcSlice(a, d);
                                //exp += primal_linear_solution[d][indexArc]*getArcSlice(a, d)/getInstance().getMaxSlice();
                            }
                        }
                    }
                }
            }
            maxUsedSliceOverallSlack3_v2[node][s] = 0 + maxUsedSliceOverall*degree - exp;
        }
    }
}

/********************************* DIRECTION ***********************************/

void AbstractLagFormulation::updateLengthDirection(double theta){

    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double value = getLengthSlack_k(d) + theta*getLengthDirection_k(d);
        setLengthDirection_k(d,value);
    }

}

void AbstractLagFormulation::updateSourceTargetDirection(double theta){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (int v = 0; v < instance.getNbNodes(); v++){ 
            double value = getSourceTargetSlack_k(d,v) + theta*getSourceTargetDirection_k(d,v);
            setSourceTargetDirection_k(d,v,value);
        }
    }
}

void AbstractLagFormulation::updateFlowDirection(double theta){
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

void AbstractLagFormulation::updateOverlapDirection(double theta){

    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
            double value = getOverlapSlack_k(e,s) + theta*getOverlapDirection_k(e,s);
            setOverlapDirection_k(e,s,value);
        }
    }
    
}

void AbstractLagFormulation::updateMaxUsedSliceOverallDirection(double theta){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        double value = getMaxUsedSliceOverallSlack_k(d) + theta*getMaxUsedSliceOverallDirection_k(d);
        setMaxUsedSliceOverallDirection_k(d,value);
    }
}

void AbstractLagFormulation::updateMaxUsedSliceOverall2Direction(double theta){
    for (int e = 0; e < instance.getNbEdges(); e++){
        for (int s = 0; s < getNbSlicesLimitFromEdge(e); s++){
            double value = getMaxUsedSliceOverall2Slack_k(e,s) + theta*getMaxUsedSliceOverall2Direction_k(e,s);
            setMaxUsedSliceOverall2Direction_k(e,s,value);
        }
    }
}

void AbstractLagFormulation::updateMaxUsedSliceOverall3Direction(double theta){
    for (int v = 0; v < instance.getNbNodes(); v++){
        for (int s = 0; s < getNbSlicesGlobalLimit(); s++){
            double value = getMaxUsedSliceOverall3Slack_k(v,s) + theta*getMaxUsedSliceOverall3Direction_k(v,s);
            setMaxUsedSliceOverall3Direction_k(v,s,value);
        }
    }
}

/* *******************************************************************************
*                             Check Feasibility METHODS
******************************************************************************* */


/* Checks if all length slacks are non-negative. */
bool AbstractLagFormulation::checkLengthFeasibility(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        if (lengthSlack[d] < -DBL_EPSILON){
            return false;
        }
    }
    return true;
}

/*  For the source and destination : its an equality, it must be 0, for the rest it must be non negative*/
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

/* Checks if all overlap1 slacks are non-negative. */
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