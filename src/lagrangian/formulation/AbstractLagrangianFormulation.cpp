#include "AbstractLagrangianFormulation.h"

AbstractLagFormulation::AbstractLagFormulation(const Instance &instance): RSA(instance){

    primal_linear_solution.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        primal_linear_solution[d].resize(countArcs(*vecGraph[d]));
        std::fill(primal_linear_solution[d].begin(), primal_linear_solution[d].end(), 0.0);
    }

}


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


   