#include "ShortestPath.h"

ShortestPath::ShortestPath(Instance &instance, const Demand &demand) : arcId(g), arcLabel(g), arcSlice(g), length(g), nodeId(g), nodeLabel(g), nodeSlice(g), onPath(g), SOURCE(demand.getSource()), TARGET(demand.getTarget()) {
    std::cout << "--- CREATING INITIAL MAPPING LEMON GRAPH... --- " << std::endl;
    // SET SOURCE NODE
    ListDigraph::Node demandSource = g.addNode();
    nodeId[demandSource] = g.id(demandSource);
    nodeLabel[demandSource] = SOURCE;
    nodeSlice[demandSource] = -1;

    // SET TARGET NODE
    ListDigraph::Node demandTarget = g.addNode();
    nodeId[demandTarget] = g.id(demandTarget);
    nodeLabel[demandTarget] = TARGET;
    nodeSlice[demandTarget] = -1;

    // FOR EACH PHYSICAL LINK (u,v) AND EACH SLICE s    
    for (int i = 0; i < instance.getNbEdges(); i++){
        int linkSourceLabel = instance.getPhysicalLinkFromId(i).getSource();
        int linkTargetLabel = instance.getPhysicalLinkFromId(i).getTarget();
        for (int s = 0; s < instance.getPhysicalLinkFromId(i).getNbSlices(); s++){
            // IF SLICE EVERY SLICE FROM s - load + 1 TO s IS NOT USED
            if (instance.hasEnoughSpace(i, s, demand)){
                // CREATE NODES (u, s) AND (v, s) IF THEY DO NOT ALREADY EXIST AND ADD AN ARC BETWEEN THEM    
                if ((linkSourceLabel != TARGET) && (linkTargetLabel != SOURCE)){
                    addArcs(linkSourceLabel, linkTargetLabel, i, s, SOURCE, TARGET, instance.getPhysicalLinkFromId(i).getLength());
                }
                if ((linkTargetLabel != TARGET) && (linkSourceLabel != SOURCE)){
                    addArcs(linkTargetLabel, linkSourceLabel, i, s, SOURCE, TARGET, instance.getPhysicalLinkFromId(i).getLength());
                }
            }
        }
    }

    //displayNodesIncidentToTarget();
    //displayNodesFromLabel(SOURCE);
}

// Creates an arc -- and nodes if necessary -- between nodes (linkSourceLabel,slice) and (linkTargetLabel,slice)
void ShortestPath::addArcs(int linkSourceLabel, int linkTargetLabel, int i, int slice, const int SOURCE, const int TARGET, double l){
    ListDigraph::Node arcSource = getNode(linkSourceLabel, slice);
    ListDigraph::Node arcTarget = getNode(linkTargetLabel, slice);
    if(linkSourceLabel == SOURCE){
        arcSource = getNode(SOURCE, -1);
    }
    if(linkTargetLabel == TARGET){
        arcTarget = getNode(TARGET,-1);
    }

    if (arcSource == INVALID){
        arcSource = g.addNode();
        nodeId[arcSource] = g.id(arcSource);
        nodeLabel[arcSource] = linkSourceLabel;
        nodeSlice[arcSource] = slice;
        //displayNode(edgeSource);
    }
    if (arcTarget == INVALID){
        arcTarget = g.addNode();
        nodeId[arcTarget] = g.id(arcTarget);
        nodeLabel[arcTarget] = linkTargetLabel;
        nodeSlice[arcTarget] = slice;
        //displayNode(edgeSource);
    }
    
    // CREATE ARC BETWEEN NODES (u, s) AND (v, s)
    if((linkSourceLabel != TARGET) && (linkTargetLabel != SOURCE)){
        ListDigraph::Arc a = g.addArc(arcSource, arcTarget);
        arcId[a] = g.id(a);
        arcLabel[a] = i;
        arcSlice[a] = slice;
        onPath[a] = false;
        length[a] = l;
        //displayEdge(a);
    }
}
ListDigraph::Node ShortestPath::getNode(int label, int slice){
    for (ListDigraph::NodeIt n(g); n != INVALID; ++n){
        if(nodeLabel[n] == label && nodeSlice[n] == slice){
            return n;
        }
    }
    return INVALID;
}


void ShortestPath::displayNodesIncidentToTarget(){
    ListDigraph::Node target = getNode(TARGET,-1);
    std::cout << "Nodes incident to target: " << std::endl;
    for (ListDigraph::InArcIt a(g, target); a != INVALID; ++a){
        ListDigraph::Node incident = g.source(a);
        std::cout << "(" << nodeLabel[incident] << ", " << nodeSlice[incident] << ")" << std::endl;
    }
}

void ShortestPath::displayNodesFromLabel(int label){
    std::cout << "Nodes with label " << label << ": " << std::endl;
    for (ListDigraph::NodeIt n(g); n != INVALID; ++n){
        if(nodeLabel[n] == label){
            std::cout << "(" << nodeLabel[n] << ", " << nodeSlice[n] << ")" << std::endl;
        }
    }
}


void ShortestPath::displayPath(){
    ListDigraph::Node source = getNode(SOURCE,-1);
    ListDigraph::Arc current = getNextOnPath(source);
    std::cout << "(" << nodeLabel[g.source(current)] + 1 << ", " << nodeSlice[g.source(current)] << ")";
    while(current != INVALID){
        std::cout << " --> (" << nodeLabel[g.target(current)] + 1 << ", " << nodeSlice[g.target(current)] << ")";
        current = getNextOnPath(g.target(current));
    }
    std::cout << std::endl;
}

ListDigraph::Arc ShortestPath::getNextOnPath(const ListDigraph::Node &n){
    for (ListDigraph::OutArcIt a(g, n); a != INVALID; ++a){
        if (onPath[a] == true){
            return a;
        }
    }
    return INVALID;
}

void ShortestPath::updateInstance(Instance &instance, const Demand &demand){
    //instance.displaySlices();
    ListDigraph::Node source = getNode(SOURCE,-1);
    ListDigraph::Arc currentArc = getNextOnPath(source);
    while(currentArc != INVALID){
        int currentEdgeLabel = arcLabel[currentArc];
        int currentEdgeSlice = arcSlice[currentArc];
        instance.assignSlicesOfLink(currentEdgeLabel, currentEdgeSlice, demand);
        //instance.getPhysicalLinkFromId(currentEdgeLabel).assignSlices(demand, currentEdgeSlice);
        currentArc = getNextOnPath(g.target(currentArc));
    }
    instance.displaySlices();
}


void ShortestPath::displayArc(const ListDigraph::Arc &a){
    std::cout << "(" << nodeLabel[g.source(a)] + 1 << ", " << nodeSlice[g.source(a)] << ")";
    std::cout << "--";
    std::cout << "(" << nodeLabel[g.target(a)] + 1 << ", " << nodeSlice[g.target(a)] << ")" << std::endl;
}
