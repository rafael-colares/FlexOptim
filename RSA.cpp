#include "RSA.h"

RSA::RSA(const Instance &inst) : instance(inst), arcId(g), arcLabel(g), arcSlice(g), length(g), nodeId(g), nodeLabel(g), nodeSlice(g), onPath(g) {
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
    
    this->setToBeRouted(instance.getNextDemands());
    displayToBeRouted();
    
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a){
        onPath[a] = -1;
    }
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


void RSA::displayArc(const ListDigraph::Arc &a){
    std::cout << "(" << nodeLabel[g.source(a)] + 1 << ", " << nodeSlice[g.source(a)] + 1 << ")";
    std::cout << "--";
    std::cout << "(" << nodeLabel[g.target(a)] + 1 << ", " << nodeSlice[g.target(a)] + 1 << ")" << std::endl;
}

void RSA::updateInstance(Instance &i){
    //instance.displaySlices();
    
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a){
        if (onPath[a] != -1){
            int id = onPath[a];
            Demand demand = i.getTabDemand()[id];
            i.assignSlicesOfLink(arcLabel[a], arcSlice[a], demand);
        }
    }
    instance.displaySlices();
}
ListDigraph::Node RSA::getFirstNodeFromLabel(int label){
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v){
        if (nodeLabel[v] == label){
            return v;
        }
    }
    return INVALID;
}

/* Contract nodes with the same given label. */
void RSA::contractNodesFromLabel(int label){
    int nb = 0;
    ListDigraph::NodeIt lastNode(g);
    ListDigraph::Node n = getFirstNodeFromLabel(label);
    nodeSlice[n] = -1;
    if (n != INVALID){
        for (ListDigraph::NodeIt v(g); v != INVALID; ++v){
            if ( (nodeLabel[v] == label) && (g.id(n) != g.id(v)) ){
                g.contract(n,v);
                v = lastNode;
                nb++;
            }
            lastNode = v;
        }
    }
    std::cout << "> Number of nodes with label " << label << " contracted: " << nb << std::endl; 
}


/* Delete arcs that are known to be unable to root demands a priori. */
void RSA::eraseNonRoutableArcs(){
    int nb = 0;
    ListDigraph::ArcIt lastArc(g);
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a){
        if (instance.isRoutable(arcLabel[a], arcSlice[a], getToBeRouted()[0]) == false){
            g.erase(a);
            a = lastArc;
            nb++;
        }
        lastArc = a;
    }
    std::cout << "> Number of non-routable arcs erased: " << nb << std::endl; 
}

double RSA::getCoeff(const ListDigraph::Arc &a, const Demand &demand){
    double coeff = 0.0;
    ListDigraph::Node s = g.source(a);
    if(nodeLabel[s] == demand.getSource()){
        coeff = arcSlice[a]+1; 
    }
    else{
        coeff = 1; 
    }
    return coeff;
}

void RSA::displayToBeRouted(){
    std::cout << "--- ROUTING DEMANDS ";
    for (int i = 0; i < getNbDemandsToBeRouted(); i++){
        std::cout << "#" << toBeRouted[i].getId()+1 << " (" << toBeRouted[i].getSource()+1 << ", " << toBeRouted[i].getTarget()+1 << "), ";
    }
    std::cout << " --- " << std::endl;
	
}

void RSA::displayNode(const ListDigraph::Node &n){
    std::cout << "(" << nodeLabel[n]+1 << "," << nodeSlice[n]+1 << ")";
}