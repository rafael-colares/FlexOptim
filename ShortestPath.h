#ifndef __ShortestPath__h
#define __ShortestPath__h

#include "Instance.h"

#include<ilcplex/ilocplex.h>
#include <lemon/list_graph.h>
#include <lemon/concepts/graph.h>

using namespace lemon;

typedef ListDigraph::NodeMap<int> NodeMap;
typedef ListDigraph::ArcMap<int> ArcMap;
typedef ListDigraph::ArcMap<double> ArcCost;

class ShortestPath{

protected:
    ListDigraph g;
    ArcMap arcId;        // provides the lemonId of an edge
    ArcMap arcLabel;     // provides the id e of an edge (e,s)
    ArcMap arcSlice;     // provides the slice s of an edge (e,s)
    ArcCost length;      // provides the length of an arc
    NodeMap nodeId;       // provides the lemonId of a node
    NodeMap nodeLabel;    // provides the id n of a node (n,s)
    NodeMap nodeSlice;    // provides the slice s of a node (n,s)
    ListDigraph::ArcMap<bool> onPath;   // provides information on whether or not an arc belong to the computed path
    const int SOURCE;
    const int TARGET;
    

public:
    // Build the extended graph from an initial mapping (instance).
    ShortestPath(Instance &instance, const Demand &demand);

    // Add an arc going from linkSourceId to linkTargetId through slice s based on the physicalLink i.
    void addArcs(int linkSourceId, int linkTargetId, int i, int s, const int SOURCE, const int TARGET, double l);
    
    // Returns the node with id (label, slice). If it does not exist, returns INVALID.
    ListDigraph::Node getNode(int label, int slice);

    void displayNodesIncidentToTarget();
    void displayNodesFromLabel(int label);
    void displayPath();
    void displayArc(const ListDigraph::Arc &a);
    
    ListDigraph::Arc getNextOnPath(const ListDigraph::Node &n);
    void updateInstance(Instance &instance, const Demand &demand);
    
};
#endif