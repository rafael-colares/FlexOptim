#include "RSA.h"

RSA::RSA(const Instance &inst) : instance(inst), compactArcId(compactGraph), compactArcLabel(compactGraph), 
                                compactArcLength(compactGraph), compactNodeId(compactGraph), 
                                compactNodeLabel(compactGraph){
    // Creates compact graph.
    for (int i = 0; i < instance.getNbNodes(); i++){
        ListDigraph::Node n = compactGraph.addNode();
        compactNodeLabel[n] = i;
        compactNodeId[n] = compactGraph.id(n);
    }
    for (int i = 0; i < instance.getNbEdges(); i++){
        PhysicalLink edge = instance.getPhysicalLinkFromId(i);
        int sourceLabel = edge.getSource();
        ListDigraph::Node sourceNode = INVALID;
        for (ListDigraph::NodeIt v(compactGraph); v != INVALID; ++v){
            if(compactNodeLabel[v] == sourceLabel){
                sourceNode = v;
            }
        }
        int targetLabel = edge.getTarget();
        ListDigraph::Node targetNode = INVALID;
        for (ListDigraph::NodeIt v(compactGraph); v != INVALID; ++v){
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
    // for (ListDigraph::NodeIt u(compactGraph); u != INVALID; ++u){
    //     for (ListDigraph::NodeIt v(compactGraph); v != INVALID; ++v){
    //         if (instance.hasLink(compactNodeLabel[u], compactNodeLabel[v])){
    //             PhysicalLink link = instance.getPhysicalLinkBetween(compactNodeLabel[u], compactNodeLabel[v]);
    //             ListDigraph::Arc a = compactGraph.addArc(u, v);
    //             compactArcId[a] = compactGraph.id(a);
    //             compactArcLabel[a] = link.getId();
    //             compactArcLength[a] = link.getLength();
    //         }
    //     }
    // }

    // Set demands to be routed.
    this->setToBeRouted(instance.getNextDemands());
    displayToBeRouted();

    // Creates an extended graph for each one of the demands to be routed.
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
                // IF SLICE s IS NOT USED
                if (instance.getPhysicalLinkFromId(i).getSlice_i(s).isUsed() == false){
                    // CREATE NODES (u, s) AND (v, s) IF THEY DO NOT ALREADY EXIST AND ADD AN ARC BETWEEN THEM   
                    addArcs(d, linkSourceLabel, linkTargetLabel, i, s, instance.getPhysicalLinkFromId(i).getLength());
                    addArcs(d, linkTargetLabel, linkSourceLabel, i, s, instance.getPhysicalLinkFromId(i).getLength());
                }
            }
        }
    }
}

/* On graph #d, creates an arc -- and its nodes if necessary -- between nodes (linkSourceLabel,slice) and (linkTargetLabel,slice) */
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

/* For graph #d, returns the node with id (label, slice). If it does not exist, returns INVALID. */
ListDigraph::Node RSA::getNode(int d, int label, int slice){
    for (ListDigraph::NodeIt n(*vecGraph[d]); n != INVALID; ++n){
        if(getNodeLabel(n, d) == label && getNodeSlice(n, d) == slice){
            return n;
        }
    }
    return INVALID;
}

/* Updates the mapping stored in instance i with the results obtained from RSA solution.*/
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

/* Returns the first node with a given label from the graph associated with the d-th demand. If such node does not exist, return INVALID. */
ListDigraph::Node RSA::getFirstNodeFromLabel(int d, int label){
    for (ListDigraph::NodeIt v(*vecGraph[d]); v != INVALID; ++v){
        if (getNodeLabel(v, d) == label){
            return v;
        }
    }
    return INVALID;
}

/* Contract nodes with the same given label from the graph associated with the d-th demand. */
void RSA::contractNodesFromLabel(int d, int label){
    int nb = 0;
    ListDigraph::NodeIt previousNode(*vecGraph[d]);
    ListDigraph::Node n = getFirstNodeFromLabel(d, label);
    (*vecNodeSlice[d])[n] = -1;
    if (n != INVALID){
        for (ListDigraph::NodeIt v(*vecGraph[d]); v != INVALID; ++v){
            if ( (getNodeLabel(v, d) == label) && ((*vecGraph[d]).id(n) != (*vecGraph[d]).id(v)) ){
                (*vecGraph[d]).contract(n,v);
                v = previousNode;
                nb++;
            }
            previousNode = v;
        }
    }
    std::cout << "> Number of nodes with label " << label << " contracted: " << nb << std::endl; 
}

/* Delete arcs that are known a priori to be unable to route on graph #d. */
void RSA::eraseNonRoutableArcs(int d){
    int nb = 0;
    ListDigraph::ArcIt previousArc(*vecGraph[d]);
    for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
        int label = getArcLabel(a, d);
        int slice = getArcSlice(a, d);
        if (instance.hasEnoughSpace(label, slice, getToBeRouted_k(d)) == false){
            (*vecGraph[d]).erase(a);
            a = previousArc;
            nb++;
        }
        previousArc = a;
    }
    std::cout << "> Number of non-routable arcs erased on graph #" << d << ": " << nb << std::endl; 
}

/* Returns the coefficient of arc a (according to the chosen metric) on graph #d. */
double RSA::getCoeff(const ListDigraph::Arc &a, int d){
    double coeff = 0.0;
    ListDigraph::Node s = (*vecGraph[d]).source(a);
    int label = getNodeLabel(s, d);
    int slice = getArcSlice(a, d);
    if(label == getToBeRouted_k(d).getSource()){
        coeff = slice + 1; 
    }
    else{
        coeff = 1; 
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

/* For graph #d, display the nodes that are incident to node (label,slice). */
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

/* For graph #d, display the nodes that have a given label. */
void RSA::displayNodesFromLabel(int d, int label){
    std::cout << "Nodes with label " << label << ": " << std::endl;
    for (ListDigraph::NodeIt n(*vecGraph[d]); n != INVALID; ++n){
        if(getNodeLabel(n, d) == label){
            std::cout << "(" << getNodeLabel(n, d) << ", " << getNodeSlice(n, d) << ")" << std::endl;
        }
    }
}

/* Display the paths found for each of the new routed demands. */
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

/* Display arc a from the graph associated with the d-th demand. */
void RSA::displayArc(int d, const ListDigraph::Arc &a){
    std::cout << "(" << getNodeLabel((*vecGraph[d]).source(a), d) + 1 << ", " <<  getNodeSlice((*vecGraph[d]).source(a), d) + 1 << ")";
    std::cout << "--";
    std::cout << "(" <<  getNodeLabel((*vecGraph[d]).target(a), d) + 1 << ", " << getNodeSlice((*vecGraph[d]).target(a), d) + 1 << ")" << std::endl;
}


/* Display node n from the graph associated with the d-th demand. */
void RSA::displayNode(int d, const ListDigraph::Node &n){
    std::cout << "(" << getNodeLabel(n, d)+1 << "," << getNodeSlice(n, d)+1 << ")";
}

