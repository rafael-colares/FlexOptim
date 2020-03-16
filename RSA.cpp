#include "RSA.h"

RSA::RSA(Instance &instance, const std::vector<Demand> &d) : arcId(g), arcLabel(g), arcSlice(g), length(g), nodeId(g), nodeLabel(g), nodeSlice(g), onPath(g), demands(d) {
    std::cout << "--- CREATING INITIAL MAPPING LEMON GRAPH... --- " << std::endl;
    // FOR EACH PHYSICAL LINK (u,v) AND EACH SLICE s    
    for (int i = 0; i < instance.getNbEdges(); i++){
        int linkSourceLabel = instance.getPhysicalLinkFromId(i).getSource();
        int linkTargetLabel = instance.getPhysicalLinkFromId(i).getTarget();
        for (int s = 0; s < instance.getPhysicalLinkFromId(i).getNbSlices(); s++){
            // IF SLICE s IS NOT USED
            if (instance.getPhysicalLinkFromId(i).getSlice_i(s).isUsed() == false){
                // CREATE NODES (u, s) AND (v, s) IF THEY DO NOT ALREADY EXIST AND ADD AN ARC BETWEEN THEM   
                addArcs(linkSourceLabel, linkTargetLabel, i, s, instance.getPhysicalLinkFromId(i).getLength());
                addArcs(linkTargetLabel, linkSourceLabel, i, s, instance.getPhysicalLinkFromId(i).getLength());
            }
        }
    }
    //displayNodesIncidentToTarget();
    //displayNodesFromLabel(SOURCE);
}

// Creates an arc -- and nodes if necessary -- between nodes (linkSourceLabel,slice) and (linkTargetLabel,slice)
void RSA::addArcs(int linkSourceLabel, int linkTargetLabel, int linkLabel, int slice, double l){
    ListDigraph::Node arcSource = getNode(linkSourceLabel, slice);
    ListDigraph::Node arcTarget = getNode(linkTargetLabel, slice);

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
    
    // CREATE ARC BETWEEN NODES arcSource and arcTarget
    ListDigraph::Arc a = g.addArc(arcSource, arcTarget);
    arcId[a] = g.id(a);
    arcLabel[a] = linkLabel;
    arcSlice[a] = slice;
    onPath[a] = false;
    length[a] = l;
    //displayEdge(a);
}

ListDigraph::Node RSA::getNode(int label, int slice){
    for (ListDigraph::NodeIt n(g); n != INVALID; ++n){
        if(nodeLabel[n] == label && nodeSlice[n] == slice){
            return n;
        }
    }
    return INVALID;
}


void RSA::displayNodesIncidentTo(int label, int slice){
    ListDigraph::Node target = getNode(label, slice);
    std::cout << "Nodes incident to (" << label << ", " << slice << "): " << std::endl;
    for (ListDigraph::InArcIt a(g, target); a != INVALID; ++a){
        ListDigraph::Node incident = g.source(a);
        std::cout << "(" << nodeLabel[incident] << ", " << nodeSlice[incident] << ")" << std::endl;
    }
    for (ListDigraph::OutArcIt a(g, target); a != INVALID; ++a){
        ListDigraph::Node incident = g.source(a);
        std::cout << "(" << nodeLabel[incident] << ", " << nodeSlice[incident] << ")" << std::endl;
    }
}

void RSA::displayNodesFromLabel(int label){
    std::cout << "Nodes with label " << label << ": " << std::endl;
    for (ListDigraph::NodeIt n(g); n != INVALID; ++n){
        if(nodeLabel[n] == label){
            std::cout << "(" << nodeLabel[n] << ", " << nodeSlice[n] << ")" << std::endl;
        }
    }
}


void RSA::displayPaths(){
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a){
        if (onPath[a] != -1){
           displayArc(a);
        }
    }
    std::cout << std::endl;
}

void RSA::updateInstance(Instance &instance, const std::vector<Demand> &demand){
    //instance.displaySlices();
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a){
        if (onPath[a] != -1){
            int demandOnPath = onPath[a];
            instance.assignSlicesOfLink(arcLabel[a], arcSlice[a], demand[demandOnPath]);
        }
    }
    instance.displaySlices();
}


void RSA::displayArc(const ListDigraph::Arc &a){
    std::cout << "(" << nodeLabel[g.source(a)] + 1 << ", " << nodeSlice[g.source(a)] << ")";
    std::cout << "--";
    std::cout << "(" << nodeLabel[g.target(a)] + 1 << ", " << nodeSlice[g.target(a)] << ")" << std::endl;
}
