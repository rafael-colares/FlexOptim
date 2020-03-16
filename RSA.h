#ifndef __RSA__h
#define __RSA__h

#include "Instance.h"

#include<ilcplex/ilocplex.h>
#include <lemon/list_graph.h>
#include <lemon/concepts/graph.h>

using namespace lemon;

typedef ListDigraph::NodeMap<int> NodeMap;
typedef ListDigraph::ArcMap<int> ArcMap;
typedef ListDigraph::ArcMap<double> ArcCost;

class RSA{

protected:
    ListDigraph g;
    ArcMap arcId;        // provides the lemonId of an edge
    ArcMap arcLabel;     // provides the id e of an edge (e,s)
    ArcMap arcSlice;     // provides the slice s of an edge (e,s)
    ArcCost length;      // provides the length of an arc
    NodeMap nodeId;       // provides the lemonId of a node
    NodeMap nodeLabel;    // provides the id n of a node (n,s)
    NodeMap nodeSlice;    // provides the slice s of a node (n,s)
    ArcMap onPath;        // provides information on which demand is assigned to the arc. 
    std::vector<Demand> demands;

public:
    // Build the extended graph from an initial mapping (instance).
    RSA(Instance &instance, const std::vector<Demand> & d);

    // Add an arc going from linkSourceLabel to linkTargetLabel through slice based on the physicalLink linkLabel.
    void addArcs(int linkSourceLabel, int linkTargetLabel, int linkLabel, int slice, double l);
    
    // Returns the node with id (label, slice). If it does not exist, returns INVALID.
    ListDigraph::Node getNode(int label, int slice);

    void displayNodesIncidentTo(int label, int slice);
    void displayNodesFromLabel(int label);
    void displayPaths();
    void displayArc(const ListDigraph::Arc &a);
    
    void updateInstance(Instance &instance, const std::vector<Demand> &demand);
    
};
#endif