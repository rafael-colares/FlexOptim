#include "shortestPathHeuristic.h"

/* *******************************************************************************
*                              CONSTRUCTOR
******************************************************************************* */

shortestPathHeuristic::shortestPathHeuristic(AbstractLagFormulation* form): AbstractHeuristic(form,STATUS_UNKNOWN){
    for (int d = 0; d < formulation->getNbDemandsToBeRouted(); d++){
        /* Initializing costs */
        heuristicCosts.emplace_back( std::make_shared<ArcCost>(*formulation->getVecGraphD(d), 0.0)); 
        heuristicCostsAux.emplace_back( std::make_shared<ArcCost>(*formulation->getVecGraphD(d), 0.0)); 
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

    switch (chosenMethod){
        case Input::VOLUME:{
            std::vector< std::vector<double> > Variables = formulation->getPrimalVariables();
            auxCosts(Variables);
            break;
        }
        case Input::SUBGRADIENT:{
            std::vector<std::vector<bool> > Variables = formulation->getVariables();
            auxCosts(Variables);
            break;
        }
        default:{
            std::vector<std::vector<bool> > Variables = formulation->getVariables();
            auxCosts(Variables);
            break;
        }
    }
    
}

void shortestPathHeuristic::auxCosts(std::vector<std::vector<bool> > Variables){
    for(int d =0; d< formulation->getNbDemandsToBeRouted();d++){
        for(ListDigraph::ArcIt arc(*formulation->getVecGraphD(d)); arc != INVALID; ++arc){
            int indexArc = formulation->getArcIndex(arc,d);
            (*heuristicCosts[d])[arc] = (1 - Variables[d][indexArc])*formulation->getCoeff(arc,d); /* cost is (1 - x)c */
            (*heuristicCostsAux[d])[arc] = (1 - Variables[d][indexArc])*formulation->getCoeff(arc,d); /* cost is (1 - x)c */
           
        }
    }

}

void shortestPathHeuristic::auxCosts(std::vector<std::vector<double> > Variables){
    for(int d =0; d< formulation->getNbDemandsToBeRouted();d++){
        for(ListDigraph::ArcIt arc(*formulation->getVecGraphD(d)); arc != INVALID; ++arc){
            int indexArc = formulation->getArcIndex(arc,d);
            (*heuristicCosts[d])[arc] = (1 - Variables[d][indexArc])*formulation->getCoeff(arc,d); /* cost is (1 - x)c */
            (*heuristicCostsAux[d])[arc] = (1 - Variables[d][indexArc])*formulation->getCoeff(arc,d); /* cost is (1 - x)c */
           
        }
    } 
}

/* *******************************************************************************
*                              RUNNING METHODS
******************************************************************************* */

void shortestPathHeuristic::run(){
    /* Initialization*/
    init();

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
        if (count  > 10000){
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
    Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > shortestPath(*formulation->getVecGraphD(d), (*heuristicCosts[d]));
    shortestPath.run(SOURCE, TARGET);

    if (shortestPath.reached(TARGET) == false){   
        std::cout << "> RSA is infeasible because there is no path from : Heuristic" << formulation->getToBeRouted_k(d).getSource()+1 << " to " << formulation->getToBeRouted_k(d).getTarget()+1 << " required for routing demand " << formulation->getToBeRouted_k(d).getId()+1 << "." << std::endl;
        statusheuristic = STATUS_INFEASIBLE;
        return true;
    }

    const int PATH_LENGTH = getPathLength(d, shortestPath, SOURCE, TARGET);
    //std::cout << "path length:" << PATH_LENGTH << std::endl;
    //std::cout << shortestPath.dist(TARGET) << std::endl;
    if(shortestPath.dist(TARGET) >= __DBL_MAX__){
        /** There is no path between the source and the destination for this demand **/
        /** It does not respect the non overlap constraints **/
        /** Remove a demand from the analysed Demands **/
        int demand2 = choseDemand(analysedDemands);
        analysedDemands.erase(demand2);
        notAnalysedDemands.insert(demand2);
        removePath_k(demand2);
        //std::cout <<"opt1" << std::endl;
        return false;

    }else if(PATH_LENGTH > formulation->getToBeRouted_k(d).getMaxLength()){
        /** The solution is infeasible because of the length constraints**/
        remove_Arc(d,shortestPath,SOURCE,TARGET);
        //std::cout <<"opt2" << std::endl;
        return false;
    }else{
        /** A feasible path for the demand was found **/
        insertPath_k(d,shortestPath,SOURCE,TARGET);
        //std::cout <<"Feasible path" << std::endl;
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
        for(ListDigraph::ArcIt arc(*formulation->getVecGraphD(k)); arc != INVALID; ++arc){
            int slice_k = formulation->getArcSlice(arc,k);
            int label_k = formulation->getArcLabel(arc,k);
            int load_k = formulation->getToBeRouted_k(k).getLoad();
            if((slice_k >= slice-load+1) && (slice_k <= slice + load_k-1) && (label_k == label)){
                (*heuristicCosts[k])[arc] = (*heuristicCostsAux[k])[arc]; /* come back with arc - their cost will be no longer infinite*/
            }
        }
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
        //std::cout << slice << " " << label << " " << load << std::endl;
        remove_arcs(slice,label,load);
        currentNode = path.predNode(currentNode);
    }
}

/* Changes the cost (to infinite) of some arcs so they will not be chosen - respect the non overlapping constraints*/
/* For all arcs, considering every demand, that have @param  slice-load+1 <= slice_k <= slice @param label_k = label */
void shortestPathHeuristic::remove_arcs(int slice, int label, int load){
    for (int k = 0; k < formulation->getNbDemandsToBeRouted(); k++){
        for(ListDigraph::ArcIt arc(*formulation->getVecGraphD(k)); arc != INVALID; ++arc){
            int slice_k = formulation->getArcSlice(arc,k);
            int label_k = formulation->getArcLabel(arc,k);
            int load_k = formulation->getToBeRouted_k(k).getLoad();
            if(!((slice_k < slice-load+1)||(slice_k-load_k+1 > slice))){
                if(label_k == label){
                    (*heuristicCosts[k])[arc] = __DBL_MAX__; /* "remove "arc - their cost will be infinite*/
                    //std::cout << (*heuristicCosts[k])[arc] << std::endl;
                }
            }
            //if((slice_k >= slice-load+1) && (slice_k <= slice + load_k-1) && (label_k == label)){
            //    (*heuristicCosts[k])[arc] = __DBL_MAX__; /* "remove "arc - their cost will be infinite*/
            //}
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
    (*heuristicCosts[d])[arcHighestLength] = __DBL_MAX__;
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

/* *******************************************************************************
*                              DESTRUCTOR
******************************************************************************* */

shortestPathHeuristic::~shortestPathHeuristic(){
    heuristicCosts.clear();
    heuristicCostsAux.clear();

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