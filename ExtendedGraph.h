#ifndef __ExtendedGraph__h
#define __ExtendedGraph__h
#include <vector>
#include <string>
#include "Instance.h"
#include "PhysicalLink.h"
#include <lemon/list_graph.h>

using namespace lemon;

typedef ListGraph::NodeMap<int> IdNodeMap;
typedef ListGraph::EdgeMap<int> LinkMap;

/************************************************************************************/
/*  ExtendedGraph models the initial mapping as a graph implemented via lemon 	    */
/*	For each non-used Slice s of each PhysicalLink (u,v),                           */
/*	we create an edge going from node (u,s) to node (v,s).							*/												
/************************************************************************************/
class ExtendedGraph{

public: 
    ListGraph g;
private:
    LinkMap edgeMap;        // provides the lemonId of an edge
    LinkMap edgeLabel;      // provides the id e of an edge (e,s)
    LinkMap edgeSlice;      // provides the slice s of an edge (e,s)
    IdNodeMap nodeId;       // provides the lemonId of a node
    IdNodeMap nodeLabel;    // provides the id n of a node (n,s)
    IdNodeMap nodeSlice;    // provides the slice s of a node (n,s)

public:
    ExtendedGraph(const Instance &instance);

	/************************************************/
	/*					Getters						*/
	/************************************************/
    int getNbNodes(){ return countNodes(g); }
    int getNbEdges(){ return countEdges(g); }
    int getLemonIdEdge(ListGraph::Edge e){ return edgeMap[e]; }
    std::string getEdgeName(ListGraph::Edge e);
    std::string getNodeName(ListGraph::Node v);

    bool checkNodeExistence(int label, int slice);
    ListGraph::Node getNode(int label, int slice);
    void displayNode (ListGraph::Node n);
    void displayEdge (ListGraph::Edge e);
};

#endif