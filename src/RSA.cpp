#include "RSA.h"

/* Constructor. A graph associated to the initial mapping (instance) is built as well as an extended graph for each demand to be routed. */
RSA::RSA(const Instance &inst) : instance(inst), compactArcId(compactGraph), compactArcLabel(compactGraph), 
                                compactArcLength(compactGraph), compactNodeId(compactGraph), 
                                compactNodeLabel(compactGraph){
    /* Creates compact graph. */
    buildCompactGraph();

    /* Set demands to be routed. */
    this->setToBeRouted(instance.getNextDemands());
    displayToBeRouted();

    /* Creates an extended graph for each one of the demands to be routed. */
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        vecGraph.emplace_back(new ListDigraph);
        vecArcId.emplace_back(new ArcMap(*vecGraph[d]));
        vecArcLabel.emplace_back(new ArcMap(*vecGraph[d]));
        vecArcSlice.emplace_back(new ArcMap(*vecGraph[d]));
        vecArcLength.emplace_back(new ArcCost(*vecGraph[d]));
        vecNodeId.emplace_back(new NodeMap(*vecGraph[d]));
        vecNodeLabel.emplace_back(new NodeMap(*vecGraph[d]));
        vecNodeSlice.emplace_back(new NodeMap(*vecGraph[d]));
        vecOnPath.emplace_back(new ArcMap(*vecGraph[d]));
    
        for (int i = 0; i < instance.getNbEdges(); i++){
            int linkSourceLabel = instance.getPhysicalLinkFromId(i).getSource();
            int linkTargetLabel = instance.getPhysicalLinkFromId(i).getTarget();
            for (int s = 0; s < instance.getPhysicalLinkFromId(i).getNbSlices(); s++){
                /* IF SLICE s IS NOT USED */
                if (instance.getPhysicalLinkFromId(i).getSlice_i(s).isUsed() == false){
                    /* CREATE NODES (u, s) AND (v, s) IF THEY DO NOT ALREADY EXIST AND ADD AN ARC BETWEEN THEM */
                    addArcs(d, linkSourceLabel, linkTargetLabel, i, s, instance.getPhysicalLinkFromId(i).getLength());
                    addArcs(d, linkTargetLabel, linkSourceLabel, i, s, instance.getPhysicalLinkFromId(i).getLength());
                }
            }
        }
    }

    /* Calls preprocessing. */
    preprocessing();
}

/* Builds the simple graph associated with the initial mapping. */
void RSA::buildCompactGraph(){
    for (int i = 0; i < instance.getNbNodes(); i++){
        ListDigraph::Node n = compactGraph.addNode();
        compactNodeLabel[n] = i;
        compactNodeId[n] = compactGraph.id(n);
    }
    for (int i = 0; i < instance.getNbEdges(); i++){
        PhysicalLink edge = instance.getPhysicalLinkFromId(i);
        int sourceLabel = edge.getSource();
        int targetLabel = edge.getTarget();
        ListDigraph::Node sourceNode = INVALID;
        ListDigraph::Node targetNode = INVALID;
        for (ListDigraph::NodeIt v(compactGraph); v != INVALID; ++v){
            if(compactNodeLabel[v] == sourceLabel){
                sourceNode = v;
            }
            if(compactNodeLabel[v] == targetLabel){
                targetNode = v;
            }
        }
        if (targetNode != INVALID && sourceNode != INVALID){
            ListDigraph::Arc a = compactGraph.addArc(sourceNode, targetNode);
            compactArcId[a] = compactGraph.id(a);
            compactArcLabel[a] = edge.getId();
            compactArcLength[a] = edge.getLength();
        }
    }
}

/* Creates an arc -- and its nodes if necessary -- between nodes (source,slice) and (target,slice) on a graph. */
void RSA::addArcs(int d, int linkSourceLabel, int linkTargetLabel, int linkLabel, int slice, double l){
    ListDigraph::Node arcSource = getNode(d, linkSourceLabel, slice);
    ListDigraph::Node arcTarget = getNode(d, linkTargetLabel, slice);

    if (arcSource == INVALID){
        arcSource = vecGraph[d]->addNode();
        int id = vecGraph[d]->id(arcSource);
        setNodeId(arcSource, d, id);
        setNodeLabel(arcSource, d, linkSourceLabel);
        setNodeSlice(arcSource, d, slice);
        //displayNode(edgeSource);
    }
    if (arcTarget == INVALID){
        arcTarget = vecGraph[d]->addNode();
        int id = vecGraph[d]->id(arcTarget);
        setNodeId(arcTarget, d, id);
        setNodeLabel(arcTarget, d, linkTargetLabel);
        setNodeSlice(arcTarget, d, slice);
        //displayNode(edgeSource);
    }
    
    // CREATE ARC BETWEEN NODES arcSource AND arcTarget
    ListDigraph::Arc a = vecGraph[d]->addArc(arcSource, arcTarget);
    int id = vecGraph[d]->id(a);
    setArcId(a, d, id);
    setArcLabel(a, d, linkLabel);
    setArcSlice(a, d, slice);
    setArcLength(a, d, l);
    (*vecOnPath[d])[a] = -1;
    //displayEdge(a);
}

/* Returns the first node identified by (label, slice) on graph #d. */
ListDigraph::Node RSA::getNode(int d, int label, int slice){
    for (ListDigraph::NodeIt n(*vecGraph[d]); n != INVALID; ++n){
        if(getNodeLabel(n, d) == label && getNodeSlice(n, d) == slice){
            return n;
        }
    }
    return INVALID;
}

/* Updates the mapping stored in the given instance with the results obtained from RSA solution (i.e., vecOnPath).*/
void RSA::updateInstance(Instance &i){
    //instance.displaySlices();
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            if ((*vecOnPath[d])[a] != -1){
                int id = (*vecOnPath[d])[a];
                Demand demand = i.getTabDemand()[id];
                i.assignSlicesOfLink(getArcLabel(a, d), getArcSlice(a, d), demand);
            }
        }
    }
    instance.displaySlices();
}

/* Returns the first node with a given label from the graph associated with the d-th demand to be routed. If such node does not exist, return INVALID. */
ListDigraph::Node RSA::getFirstNodeFromLabel(int d, int label){
    for (ListDigraph::NodeIt v(*vecGraph[d]); v != INVALID; ++v){
        if (getNodeLabel(v, d) == label){
            return v;
        }
    }
    return INVALID;
}

/* Contract nodes with the same given label from the graph associated with the d-th demand to be routed. */
void RSA::contractNodesFromLabel(int d, int label){
    int nb = 0;
    ListDigraph::NodeIt previousNode(*vecGraph[d]);
    ListDigraph::Node n = getFirstNodeFromLabel(d, label);
    ListDigraph::NodeIt v(*vecGraph[d]);
    ListDigraph::NodeIt currentNode(*vecGraph[d], v);
    if (n != INVALID){
        (*vecNodeSlice[d])[n] = -1;
        while (v != INVALID){
            currentNode = v;
            ListDigraph::NodeIt nextNode(*vecGraph[0], ++currentNode);
            currentNode = v;
            if ( (getNodeLabel(v, d) == label) && ((*vecGraph[d]).id(n) != (*vecGraph[d]).id(v)) ){
                (*vecGraph[0]).contract(n, v);
                nb++;
            }
            v = nextNode;
        }
    }
    std::cout << "> Number of nodes with label " << label << " contracted: " << nb << std::endl; 
}

/* Delete arcs that are known 'a priori' to be unable to route on graph #d. */
void RSA::eraseNonRoutableArcs(int d){
    int nb = 0;
    ListDigraph::ArcIt previousArc(*vecGraph[d]);
    ListDigraph::ArcIt a(*vecGraph[d]);
    ListDigraph::ArcIt currentArc(*vecGraph[d], a);
    int demandSource = getToBeRouted_k(d).getSource();
    int demandTarget = getToBeRouted_k(d).getTarget();
    while (a != INVALID){
        currentArc = a;
        ListDigraph::ArcIt nextArc(*vecGraph[d], ++currentArc);
        currentArc = a;
        int label = getArcLabel(a, d);
        int slice = getArcSlice(a, d);
        int uLabel = getNodeLabel((*vecGraph[d]).source(a), d);
        int vLabel = getNodeLabel((*vecGraph[d]).target(a), d);
        if ( (instance.hasEnoughSpace(label, slice, getToBeRouted_k(d)) == false) || (uLabel == demandTarget) || (vLabel == demandSource) ){
            (*vecGraph[d]).erase(a);
            nb++;
        }
        a = nextArc;
    }
    std::cout << "> Number of non-routable arcs erased on graph #" << d << ": " << nb << std::endl; 
}

/* Runs preprocessing on every extended graph. */
void RSA::preprocessing(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        std::cout << "> Number of arcs in graph #" << d << ": " << countArcs((*vecGraph[d])) << std::endl;
    }
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        eraseNonRoutableArcs(d);
    }
    if (getInstance().getInput().getChosenPreprLvl() >= Input::PREPROCESSING_LVL_PARTIAL){
        // do partial preprocessing;
        bool keepPreprocessing = lengthPreprocessing();
        
        if (getInstance().getInput().getChosenPreprLvl() >= Input::PREPROCESSING_LVL_FULL){
            // do full preprocessing;
            while (keepPreprocessing){
                keepPreprocessing = lengthPreprocessing();
            }
        }
    }
}


/* Erases every arc from graph #d having the given slice and returns the number of arcs removed. */
int RSA::eraseAllArcsFromSlice(int d, int slice){
    int nb = 0;
    ListDigraph::ArcIt previousArc(*vecGraph[d]);
    ListDigraph::ArcIt a(*vecGraph[d]);
    ListDigraph::ArcIt currentArc(*vecGraph[d], a);
    while (a != INVALID){
        currentArc = a;
        ListDigraph::ArcIt nextArc(*vecGraph[0], ++currentArc);
        currentArc = a;
        if (getArcSlice(a, d) == slice){
            //std::cout << "Erase arc ";
            //displayArc(d, a);
            (*vecGraph[d]).erase(a);
            nb++;
        }
        a = nextArc;
    }
    return nb;
}
/* Performs preprocessing based on the arc lengths and returns true if at least one arc is erased. */
bool RSA::lengthPreprocessing(){
    int totalNb = 0;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        //displayGraph(d);
        int nb = 0;
        ListDigraph::ArcIt previousArc(*vecGraph[d]);
        ListDigraph::ArcIt a(*vecGraph[d]);
        ListDigraph::ArcIt currentArc(*vecGraph[d], a);
        while (a != INVALID){
            currentArc = a;
            ListDigraph::ArcIt nextArc(*vecGraph[0], ++currentArc);
            currentArc = a;
            int slice = getArcSlice(a, d);
            ListDigraph::Node source = getNode(d, getToBeRouted_k(d).getSource(), slice);
            ListDigraph::Node target = getNode(d, getToBeRouted_k(d).getTarget(), slice);
            if (source != INVALID && target != INVALID){
                if (shortestDistance(d, source, a, target) >= getToBeRouted_k(d).getMaxLength() + DBL_EPSILON){
                    (*vecGraph[d]).erase(a);
                    nb++;
                }
            }
            else{
                nb += eraseAllArcsFromSlice(d, slice);
            }
            a = nextArc;
        }
        std::cout << "> Number of erased arcs due to length in graph #" << d << ": " << nb << std::endl;
        totalNb += nb;
    }
    if (totalNb >= 1){
        return true;
    }
    return false;
}

/* Returns the distance of the shortest path from source to target passing through arc a. */
double RSA::shortestDistance(int d, ListDigraph::Node &s, ListDigraph::Arc &a, ListDigraph::Node &t){
    double distance = 0.0;
    Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > s_u_path((*vecGraph[d]), (*vecArcLength[d]));

    s_u_path.run(s,(*vecGraph[d]).source(a));
    if (s_u_path.reached((*vecGraph[d]).source(a))){
        distance += s_u_path.dist((*vecGraph[d]).source(a));
    }
    else{
        return DBL_MAX;
    }
    
    distance += getArcLength(a, d);

    Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > v_t_path((*vecGraph[d]), (*vecArcLength[d]));
    v_t_path.run((*vecGraph[d]).target(a), t);
    if (v_t_path.reached(t)){
        distance += v_t_path.dist(t);
    }
    else{
        return DBL_MAX;
    }

    return distance;
}

/* Returns the coefficient of an arc according to metric 1 on graph #d. */
double RSA::getCoeffObj1(const ListDigraph::Arc &a, int d){
    double coeff = 0.0;
    ListDigraph::Node u = (*vecGraph[d]).source(a);
    int uLabel = getNodeLabel(u, d);
    int arcSlice = getArcSlice(a, d);
    if(uLabel == getToBeRouted_k(d).getSource()){
        coeff = arcSlice + 1; 
    }
    else{
        coeff = 1; 
    }
    return coeff;
}

/* Returns the coefficient of an arc according to metric 1p on graph #d. */
double RSA::getCoeffObj1p(const ListDigraph::Arc &a, int d){
    double coeff = 0.0;
    int arcLabel = getArcLabel(a, d);
    int arcSlice = getArcSlice(a, d);
    int maxSliceUsedOnLink = instance.getPhysicalLinkFromId(arcLabel).getMaxUsedSlicePosition();
    if(arcSlice <= maxSliceUsedOnLink){
        coeff = maxSliceUsedOnLink; 
    }
    else{
        coeff = arcSlice; 
    }
    return coeff;
}

/* Returns the coefficient of an arc according to metric 2 on graph #d. */
double RSA::getCoeffObj2(const ListDigraph::Arc &a, int d){
    double coeff = 1.0;
    return coeff;
}

/* Returns the coefficient of an arc according to metric 4 on graph #d. */
double RSA::getCoeffObj4(const ListDigraph::Arc &a, int d){
    return getArcLength(a, d);
}

/* Returns the coefficient of an arc according to metric 8 on graph #d. */
double RSA::getCoeffObj8(const ListDigraph::Arc &a, int d){
    double coeff = 0.0;
    int maxSliceUsed = 0;
    for (int i = 0; i < instance.getNbEdges(); i++){
        int maxSliceUsedOnLink = instance.getPhysicalLinkFromId(i).getMaxUsedSlicePosition();
        if (maxSliceUsedOnLink >= maxSliceUsed){
            maxSliceUsed = maxSliceUsedOnLink;
        }
    }
    int arcSlice = getArcSlice(a, d);
    ListDigraph::Node u = (*vecGraph[d]).source(a);
    int uLabel = getNodeLabel(u, d);
    if(uLabel == getToBeRouted_k(d).getSource()){
        if(arcSlice <= maxSliceUsed){
            coeff = maxSliceUsed + 1; 
        }
        else{
            coeff = arcSlice + 1; 
        }
    }
    else{
        coeff = 1;
    }
    return coeff;
}

/* Returns the coefficient of an arc (according to the chosen metric) on graph #d. */
double RSA::getCoeff(const ListDigraph::Arc &a, int d){
    double coeff = 0.0;
    switch (getInstance().getInput().getChosenObj()){
        case Input::OBJECTIVE_METRIC_1:
        {
            coeff = getCoeffObj1(a, d);
            break;
        }
        case Input::OBJECTIVE_METRIC_1p:
        {
            coeff = getCoeffObj1p(a, d);
            break;
        }
        case Input::OBJECTIVE_METRIC_2:
        {
            coeff = getCoeffObj2(a, d);
            break;
        }
        case Input::OBJECTIVE_METRIC_4:
        {
            coeff = getCoeffObj4(a, d);
            break;
        }
        case Input::OBJECTIVE_METRIC_8:
        {
            coeff = getCoeffObj8(a, d);
            break;
        }
        default:
        {
            std::cerr << "Objective metric out of range.\n";
            exit(0);
            break;
        }
    }
    
    return coeff;
}

/* Displays the demands to be routed in the next optimization. */
void RSA::displayToBeRouted(){
    std::cout << "--- ROUTING DEMANDS ";
    for (int i = 0; i < getNbDemandsToBeRouted(); i++){
        std::cout << "#" << toBeRouted[i].getId()+1 << " (" << toBeRouted[i].getSource()+1 << ", " << toBeRouted[i].getTarget()+1 << "), ";
    }
    std::cout << " --- " << std::endl;
	
}

/* Displays the nodes of graph #d that are incident to the node identified by (label,slice). */
void RSA::displayNodesIncidentTo(int d, int label, int slice){
    ListDigraph::Node target = getNode(d, label, slice);
    std::cout << "Nodes incident to (" << label << ", " << slice << "): " << std::endl;
    for (ListDigraph::InArcIt a(*vecGraph[d], target); a != INVALID; ++a){
        ListDigraph::Node incident = (*vecGraph[d]).source(a);
        std::cout << "(" << getNodeLabel(incident, d) << ", " << getNodeSlice(incident, d) << ")" << std::endl;
    }
    for (ListDigraph::OutArcIt a(*vecGraph[d], target); a != INVALID; ++a){
        ListDigraph::Node incident = (*vecGraph[d]).source(a);
        std::cout << "(" <<  getNodeLabel(incident, d) << ", " << getNodeSlice(incident, d) << ")" << std::endl;
    }
}

/* Displays the nodes of graph #d that have a given label. */
void RSA::displayNodesFromLabel(int d, int label){
    std::cout << "Nodes with label " << label << ": " << std::endl;
    for (ListDigraph::NodeIt n(*vecGraph[d]); n != INVALID; ++n){
        if(getNodeLabel(n, d) == label){
            std::cout << "(" << getNodeLabel(n, d) << ", " << getNodeSlice(n, d) << ")" << std::endl;
        }
    }
}

/* Displays the paths found for each of the new routed demands. */
void RSA::displayPaths(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            if ((*vecOnPath[d])[a] != -1){
                displayArc(d, a);
            }
        }
        std::cout << std::endl;
    }
}

/* Displays an arc from the graph #d. */
void RSA::displayArc(int d, const ListDigraph::Arc &a){
    std::cout << "(" << getNodeLabel((*vecGraph[d]).source(a), d) + 1 << ", " <<  getNodeSlice((*vecGraph[d]).source(a), d) + 1 << ")";
    std::cout << "--";
    std::cout << "(" <<  getNodeLabel((*vecGraph[d]).target(a), d) + 1 << ", " << getNodeSlice((*vecGraph[d]).target(a), d) + 1 << ")" << std::endl;
}


/* Display all arcs from the graph #d. */
void RSA::displayGraph(int d){
    for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
        displayArc(d, a);
    }
}

/* Displays a node from the graph #d. */
void RSA::displayNode(int d, const ListDigraph::Node &n){
    std::cout << "(" << getNodeLabel(n, d)+1 << "," << getNodeSlice(n, d)+1 << ")";
}
