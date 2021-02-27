#include "AbstractHeuristic.h"

/* *******************************************************************************
*                              CONSTRUCTOR
******************************************************************************* */

AbstractHeuristic::AbstractHeuristic(AbstractLagFormulation* form, const Status &s):formulation(form),statusheuristic(s),time(ClockTime::getTimeNow()){
    /* Initializing solution */
    heuristicSolution.resize(formulation->getNbDemandsToBeRouted());
    for (int d = 0; d < formulation->getNbDemandsToBeRouted(); d++){
        heuristicSolution[d].resize(countArcs(*formulation->getVecGraphD(d)),false);
        heuristicSolutionItBoolMap.emplace_back(std::make_shared<IterableBoolMap<ListDigraph, ListDigraph::Arc>>((*formulation->getVecGraphD(d)),false)); //mapFill<ListDigraph,IterableBoolMap<ListDigraph,ListDigraph::Arc>>((*formulation->getVecGraphD(d)),(*heuristicSolutionItBoolMap[d]),false);
    } 
    timeAux = 0.0;
}

/* *******************************************************************************
*                              INITIALIZATION
******************************************************************************* */

/*Initializes the  heuristic Solution - to 0 for every element - to run the algorithm*/
void AbstractHeuristic::initSolution(){  
    for (int d = 0; d < formulation->getNbDemandsToBeRouted(); d++){
        std::fill(heuristicSolution[d].begin(), heuristicSolution[d].end(), false);
        mapFill<ListDigraph,IterableBoolMap<ListDigraph,ListDigraph::Arc>>((*formulation->getVecGraphD(d)),(*heuristicSolutionItBoolMap[d]),false);
    }
}

void AbstractHeuristic::updateCostFromHeuristic(){
    setCurrentHeuristicCost(0.0);
    if(maxChangesPossible){
        setCurrentHeuristicCost(__DBL_MAX__);
    }
    else{
        Input::ObjectiveMetric chosenMetric = formulation->getInstance().getInput().getChosenObj_k(0);
        if(chosenMetric == Input::OBJECTIVE_METRIC_8){
            double val = 0.0;
            for (int d = 0; d < formulation->getNbDemandsToBeRouted(); d++){
                for (ListDigraph::ArcIt a(*formulation->getVecGraphD(d)); a != INVALID; ++a){
                    int index = formulation->getArcIndex(a,d);
                    if(heuristicSolution[d][index] == true && formulation->getCoeff(a,d) > val){
                        val = formulation->getCoeff(a,d);  
                    }
                }
            }
            incCurrentHeuristicCost(val);
        }else{
            for (int d = 0; d < formulation->getNbDemandsToBeRouted(); d++){
                for (ListDigraph::ArcIt a(*formulation->getVecGraphD(d)); a != INVALID; ++a){
                    int index = formulation->getArcIndex(a,d);
                    if(heuristicSolution[d][index] == true){
                        double val = formulation->getCoeff(a,d);
                        incCurrentHeuristicCost(val);
                    }
                }
            }  
        }
    }
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

AbstractHeuristic::~AbstractHeuristic(){
    heuristicSolution.clear();
}


