#include "AbstractHeuristic.h"

/* *******************************************************************************
*                              CONSTRUCTOR
******************************************************************************* */

AbstractHeuristic::AbstractHeuristic(AbstractLagFormulation* form, const Status &s):formulation(form),statusheuristic(s){
    heuristicSolution.resize(formulation->getNbDemandsToBeRouted());
    for (int d = 0; d < formulation->getNbDemandsToBeRouted(); d++){
        /* Initializing solution */
        heuristicSolution[d].resize(countArcs(*formulation->getVecGraphD(d)));
        std::fill(heuristicSolution[d].begin(), heuristicSolution[d].end(), false);
    } 
}

/* *******************************************************************************
*                              INITIALIZATION
******************************************************************************* */

/*Initializes the  heuristic Solution - to 0 for every element - to run the algorithm*/
void AbstractHeuristic::initSolution(){  
    for (int d = 0; d < formulation->getNbDemandsToBeRouted(); d++){
        std::fill(heuristicSolution[d].begin(), heuristicSolution[d].end(), false);
    }
}

void AbstractHeuristic::updateCostFromHeuristic(){
    setCurrentHeuristicCost(0.0);
    int somaval =0;
    for (int d = 0; d < formulation->getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(*formulation->getVecGraphD(d)); a != INVALID; ++a){
            int index = formulation->getArcIndex(a,d);
            if(heuristicSolution[d][index] == true){
                int val = formulation->getCoeff(a,d);
                somaval += val;
                //std::cout<< val << std::endl;
                incCurrentHeuristicCost(val);
            }
        }
    }
    //std::cout << getCurrentHeuristicCost() << std::endl;
    //printSolution();
}

void AbstractHeuristic::printSolution(){
    int somaslice = 0;
    for (int d = 0; d < formulation->getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(*formulation->getVecGraphD(d)); a != INVALID; ++a){
            int index = formulation->getArcIndex(a,d);
            if(heuristicSolution[d][index] == true){
                std::cout << "label: " << formulation->getArcLabel(a,d)+1 << " slice:" << formulation->getArcSlice(a,d) << " demand: " << d+1 << " load" << formulation->getToBeRouted_k(d).getLoad()<< std::endl;
                somaslice += formulation->getArcSlice(a,d);
            }
        }
    }
    std::cout<< "soma: " << somaslice << std::endl;

}
