#include "shortestPathHeuristic.h"

/* *******************************************************************************
*                              CONSTRUCTOR
******************************************************************************* */

shortestPathHeuristic::shortestPathHeuristic(AbstractLagFormulation* form): AbstractHeuristic(form,STATUS_UNKNOWN){
    for (int d = 0; d < formulation->getNbDemandsToBeRouted(); d++){
        /* Initializing costs */
        heuristicCosts.emplace_back( std::make_shared<ArcCost>(*formulation->getVecGraphD(d), 0.0)); 
        freeArcs.emplace_back(std::make_shared<ListDigraph::ArcMap<bool>>(*formulation->getVecGraphD(d),true));
        mapItLabel.emplace_back(std::make_shared<IterableIntMap<ListDigraph, ListDigraph::Arc>>(*formulation->getVecGraphD(d)));
        mapCopy<ListDigraph,ArcMap,IterableIntMap<ListDigraph, ListDigraph::Arc>>((*formulation->getVecGraphD(d)),(*formulation->getArcLabelMap(d)),(*mapItLabel[d]));
    }  
}

/* *******************************************************************************
*                              INITIALIZATION
******************************************************************************* */

/* Initializes the Heuristics elements to run the algorithm */
void shortestPathHeuristic::init(){  
    initSolution();
    initDemandsSets();
    initCosts();
    statusheuristic = STATUS_UNKNOWN;
}

/* Initialize the sets of not analysed demands and analysed demands for the Heuristic - to run the algorithm*/
void shortestPathHeuristic::initDemandsSets(){
    notAnalysedDemands.clear();
    analysedDemands.clear();
    for (int d = 0; d < formulation->getNbDemandsToBeRouted(); d++){ notAnalysedDemands.insert(d);}
}

/** Initialize the considered cost for the heuristic using the values of the assgnmentMatrix. cost arc = 1-x - to run the algorithm **/
void shortestPathHeuristic::initCosts(){
    Input::LagMethod chosenMethod = formulation->getInstance().getInput().getChosenLagMethod();
    Input::ObjectiveMetric chosenMetric = formulation->getInstance().getInput().getChosenObj_k(0);
    switch (chosenMethod){
        case Input::VOLUME:{
            HeuristicCostVar operPrimal(formulation->getPrimalVariables()); // operator to combine the information
            if(chosenMetric == Input::OBJECTIVE_METRIC_8){
                for (int d = 0; d < formulation->getNbDemandsToBeRouted(); d++){
                    operPrimal.setDemand(d);
                    CombineHeuristicCostPrimal combPrimal((*formulation->getCoeffMapObj8(d)),(*formulation->getIndexMap(d)),operPrimal);
                    mapCopy<ListDigraph,CombineHeuristicCostPrimal,ArcCost>((*formulation->getVecGraphD(d)),combPrimal,(*heuristicCosts[d]));
                    mapFill<ListDigraph,ListDigraph::ArcMap<bool>>((*formulation->getVecGraphD(d)),(*freeArcs[d]),true);
                }
            }else{
                for (int d = 0; d < formulation->getNbDemandsToBeRouted(); d++){
                    operPrimal.setDemand(d);
                    CombineHeuristicCostPrimal combPrimal((*formulation->getCoeffMap(d)),(*formulation->getIndexMap(d)),operPrimal);
                    mapCopy<ListDigraph,CombineHeuristicCostPrimal,ArcCost>((*formulation->getVecGraphD(d)),combPrimal,(*heuristicCosts[d]));
                    mapFill<ListDigraph,ListDigraph::ArcMap<bool>>((*formulation->getVecGraphD(d)),(*freeArcs[d]),true);
                }
            }
            break;
        }
        case Input::SUBGRADIENT:{
            HeuristicCostAssign operVar(formulation->getVariables());
            if(chosenMetric == Input::OBJECTIVE_METRIC_8){
                for (int d = 0; d < formulation->getNbDemandsToBeRouted(); d++){
                    operVar.setDemand(d);
                    CombineHeuristicCostAssign combVar((*formulation->getCoeffMapObj8(d)),(*formulation->getIndexMap(d)),operVar);
                    mapCopy<ListDigraph,CombineHeuristicCostAssign,ArcCost>((*formulation->getVecGraphD(d)),combVar,(*heuristicCosts[d]));
                    mapFill<ListDigraph,ListDigraph::ArcMap<bool>>((*formulation->getVecGraphD(d)),(*freeArcs[d]),true);
                }
            }else{    
                for (int d = 0; d < formulation->getNbDemandsToBeRouted(); d++){
                    operVar.setDemand(d);
                    CombineHeuristicCostAssign combVar((*formulation->getCoeffMap(d)),(*formulation->getIndexMap(d)),operVar);
                    mapCopy<ListDigraph,CombineHeuristicCostAssign,ArcCost>((*formulation->getVecGraphD(d)),combVar,(*heuristicCosts[d]));
                    mapFill<ListDigraph,ListDigraph::ArcMap<bool>>((*formulation->getVecGraphD(d)),(*freeArcs[d]),true);
                }
            }
            break;
        }
        default:{
            HeuristicCostAssign operVar(formulation->getVariables());
            if(chosenMetric == Input::OBJECTIVE_METRIC_8){
                for (int d = 0; d < formulation->getNbDemandsToBeRouted(); d++){
                    operVar.setDemand(d);
                    CombineHeuristicCostAssign combVar((*formulation->getCoeffMapObj8(d)),(*formulation->getIndexMap(d)),operVar);
                    mapCopy<ListDigraph,CombineHeuristicCostAssign,ArcCost>((*formulation->getVecGraphD(d)),combVar,(*heuristicCosts[d]));
                    mapFill<ListDigraph,ListDigraph::ArcMap<bool>>((*formulation->getVecGraphD(d)),(*freeArcs[d]),true);
                }
            }else{    
                for (int d = 0; d < formulation->getNbDemandsToBeRouted(); d++){
                    operVar.setDemand(d);
                    CombineHeuristicCostAssign combVar((*formulation->getCoeffMap(d)),(*formulation->getIndexMap(d)),operVar);
                    mapCopy<ListDigraph,CombineHeuristicCostAssign,ArcCost>((*formulation->getVecGraphD(d)),combVar,(*heuristicCosts[d]));
                    mapFill<ListDigraph,ListDigraph::ArcMap<bool>>((*formulation->getVecGraphD(d)),(*freeArcs[d]),true);
                }
            }
            break;
        }
    }
    
}

/* *******************************************************************************
*                              RUNNING METHODS
******************************************************************************* */

void shortestPathHeuristic::run(){
    /* Initialization*/
    init();
    
    int maxChanges = 10*formulation->getNbDemandsToBeRouted();
    maxChangesPossible = false;
    
    int count = 0;
    while(!notAnalysedDemands.empty()){
        /* Choses a demand */
        int demand = choseDemand(notAnalysedDemands);
        notAnalysedDemands.erase(demand);
        std::pair<std::set<int>::iterator,bool> aux = analysedDemands.insert(demand);
        bool STOP = false;

        while(!STOP){
            /* Apply the shortest path */       
            STOP = heuristicRun(demand);  
            count ++;   
        }
        if(statusheuristic == STATUS_INFEASIBLE){
            std::cout <<"Infeasible\n";
            break;
        }
        if (count  > maxChanges){
            maxChangesPossible = true;
            break;
        }
    }
    updateCostFromHeuristic();
    
}

/* Find a feasible shortest path  for demand d - return true*/
/* If coul not find the path - return false */
bool shortestPathHeuristic::heuristicRun(int d){
    
    /** SOURCE and DEST from the graph d **/
    const ListDigraph::Node SOURCE = formulation->getFirstNodeFromLabel(d, formulation->getToBeRouted_k(d).getSource());
    const ListDigraph::Node TARGET = formulation->getFirstNodeFromLabel(d, formulation->getToBeRouted_k(d).getTarget());

    /** Shortest path with the heuristic cost**/
    FilterArcs<ListDigraph> subgraph(*formulation->getVecGraphD(d), (*freeArcs[d]));
    Dijkstra< FilterArcs<ListDigraph>, ListDigraph::ArcMap<double> > shortestPath(subgraph, (*heuristicCosts[d]));
    //Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > shortestPath(*formulation->getVecGraphD(d), (*heuristicCosts[d]));
    shortestPath.run(SOURCE, TARGET);
   
    time.setStart(ClockTime::getTimeNow());
    if (shortestPath.reached(TARGET) == false){   
        /** There is no path between the source and the destination for this demand **/
        /** It does not respect the non overlap constraints **/
        /** Remove a demand from the analysed Demands **/
        int demand2 = choseDemand(analysedDemands);
        analysedDemands.erase(demand2);
        notAnalysedDemands.insert(demand2);
        removePath_k(demand2);
        return false;

    }else if(getPathLength(d, shortestPath, SOURCE, TARGET) > formulation->getToBeRouted_k(d).getMaxLength()){
        /** The solution is infeasible because of the length constraints**/
        remove_Arc(d,shortestPath,SOURCE,TARGET);
        return false;
    }else{
        /** A feasible path for the demand was found **/
        insertPath_k(d,shortestPath,SOURCE,TARGET);
        timeAux += time.getTimeInSecFromStart();
        return true;
    }
}

/* Remove a found path for demand d. Then, we have to find another path for this demand*/
void shortestPathHeuristic::removePath_k(int d){
    for(ListDigraph::ArcIt a(*formulation->getVecGraphD(d)); a != INVALID; ++a){
        int index = formulation->getArcIndex(a,d);
        if(heuristicSolution[d][index] == true){
            /* the path uses this arc */
            heuristicSolution[d][index] = false; /* removing path from the solution */
            int slice = formulation->getArcSlice(a,d);
            int label = formulation->getArcLabel(a,d);
            int load = formulation->getToBeRouted_k(d).getLoad();
            /* include the arcs were removed when this path was chosen */
            /* they were removed from the posibilities in order to find feasible solutions considering non overlap constraints*/
            include_arcs(slice,label,load); 
        }
    }
}

/* Changes the cost of some arcs so they can be chosen - respect the non overlapping constraints*/
/* For all arcs, considering every demand, that have @param  slice-load+1 <= slice_k <= slice @param label_k = label */
void shortestPathHeuristic::include_arcs(int slice, int label, int load){
    for (int k = 0; k < formulation->getNbDemandsToBeRouted(); k++){
        int load_k = formulation->getToBeRouted_k(k).getLoad();
        for(IterableIntMap< ListDigraph, ListDigraph::Arc >::ItemIt arc((*mapItLabel[k]),label); arc != INVALID; ++arc){
            int slice_k = formulation->getArcSlice(arc,k);
            if((slice_k >= slice-load+1) && (slice_k <= slice + load_k-1)){
                (*freeArcs[k])[arc] = true;
            }
        }
    }
}

/* Changes the heuristic Solution considering the found path for demand d*/
void shortestPathHeuristic::insertPath_k(int d, Dijkstra< FilterArcs<ListDigraph>, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &SOURCE, const ListDigraph::Node &TARGET){
    ListDigraph::Node currentNode = TARGET;
    while (currentNode != SOURCE){
        ListDigraph::Arc arc = path.predArc(currentNode);
        int index = formulation->getArcIndex(arc, d);
        heuristicSolution[d][index] = true;
        int slice = formulation->getArcSlice(arc,d);
        int label = formulation->getArcLabel(arc,d);
        int load = formulation->getToBeRouted_k(d).getLoad();
        remove_arcs(slice,label,load);
        currentNode = path.predNode(currentNode);
    }
}

/* Changes the heuristic Solution considering the found path for demand d*/
void shortestPathHeuristic::insertPath_k(int d, Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &SOURCE, const ListDigraph::Node &TARGET){
    ListDigraph::Node currentNode = TARGET;
    while (currentNode != SOURCE){
        ListDigraph::Arc arc = path.predArc(currentNode);
        int index = formulation->getArcIndex(arc, d);
        heuristicSolution[d][index] = true;
        int slice = formulation->getArcSlice(arc,d);
        int label = formulation->getArcLabel(arc,d);
        int load = formulation->getToBeRouted_k(d).getLoad();
        remove_arcs(slice,label,load);
        currentNode = path.predNode(currentNode);
    }
}

/* Changes the cost (to infinite) of some arcs so they will not be chosen - respect the non overlapping constraints*/
/* For all arcs, considering every demand, that have @param  slice-load+1 <= slice_k <= slice @param label_k = label */
void shortestPathHeuristic::remove_arcs(int slice, int label, int load){    
    for (int k = 0; k < formulation->getNbDemandsToBeRouted(); k++){
        int load_k = formulation->getToBeRouted_k(k).getLoad();
        for(IterableIntMap< ListDigraph, ListDigraph::Arc >::ItemIt arc((*mapItLabel[k]),label); arc != INVALID; ++arc){
            int slice_k = formulation->getArcSlice(arc,k);
            if(!((slice_k < slice-load+1)||(slice_k-load_k+1 > slice))){
                (*freeArcs[k])[arc] = false;
            }
        }
    }
}

/* "Remove" (cost infinite) arc with highest length, so  it can not be selected -> respect length constraints  */
void shortestPathHeuristic::remove_Arc(int d, Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &SOURCE, const ListDigraph::Node &TARGET){
    ListDigraph::Node currentNode = TARGET;
    ListDigraph::Arc arcHighestLength = path.predArc(currentNode);
    int length = formulation->getArcLength(arcHighestLength,d);
    while (currentNode != SOURCE){
        ListDigraph::Arc arc = path.predArc(currentNode);
        if(formulation->getArcLength(arc,d) > length){
            arcHighestLength = arc;
        }
        currentNode = path.predNode(currentNode);
    }
    (*freeArcs[d])[arcHighestLength] = false;
}

/* "Remove" (cost infinite) arc with highest length, so  it can not be selected -> respect length constraints  */
void shortestPathHeuristic::remove_Arc(int d, Dijkstra< FilterArcs<ListDigraph>, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &SOURCE, const ListDigraph::Node &TARGET){
    ListDigraph::Node currentNode = TARGET;
    ListDigraph::Arc arcHighestLength = path.predArc(currentNode);
    int length = formulation->getArcLength(arcHighestLength,d);
    while (currentNode != SOURCE){
        ListDigraph::Arc arc = path.predArc(currentNode);
        if(formulation->getArcLength(arc,d) > length){
            arcHighestLength = arc;
        }
        currentNode = path.predNode(currentNode);
    }
    (*freeArcs[d])[arcHighestLength] = false;
}

/* Selects one demand to be analysed */
int shortestPathHeuristic::choseDemand(std::set<int> set){
    std::set<int>::iterator it;
    int load = 0;
    int demand = -1;
    for (it = set.begin(); it != set.end(); ++it){
        if(formulation->getToBeRouted_k(*it).getLoad() > load){
            load = formulation->getToBeRouted_k(*it).getLoad();
            demand = *it;
        }
    }
    return demand;
}

/* Returns the physical length of the path. */
double shortestPathHeuristic::getPathLength(int d, Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &s, const ListDigraph::Node &t){
    double pathLength = 0.0;
    ListDigraph::Node n = t;
    while (n != s){
        ListDigraph::Arc arc = path.predArc(n);
        n = path.predNode(n);
        pathLength += formulation->getArcLength(arc, d);
    }
    return pathLength;
}

/* Returns the physical length of the path. */
double shortestPathHeuristic::getPathLength(int d, Dijkstra< FilterArcs<ListDigraph>, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &s, const ListDigraph::Node &t){
    double pathLength = 0.0;
    ListDigraph::Node n = t;
    while (n != s){
        ListDigraph::Arc arc = path.predArc(n);
        n = path.predNode(n);
        pathLength += formulation->getArcLength(arc, d);
    }
    return pathLength;
}

/* *******************************************************************************
*                              DESTRUCTOR
******************************************************************************* */

shortestPathHeuristic::~shortestPathHeuristic(){
    heuristicCosts.clear();
    freeArcs.clear();
    mapItLabel.clear();

    while(!notAnalysedDemands.empty()){
        std::set<int>::iterator it;
        it = std::prev(notAnalysedDemands.end()); 
        
        notAnalysedDemands.erase(it);
    }
    notAnalysedDemands.clear();
    while(!analysedDemands.empty()){
        std::set<int>::iterator it;
        it = std::prev(analysedDemands.end()); 
        analysedDemands.erase(it);
    }
    analysedDemands.clear();
}