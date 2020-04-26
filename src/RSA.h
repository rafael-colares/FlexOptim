#ifndef __RSA__h
#define __RSA__h

#include "Instance.h"

#include<ilcplex/ilocplex.h>
#include <lemon/list_graph.h>
#include <lemon/concepts/graph.h>
#include <lemon/dijkstra.h>
#include <lemon/list_graph.h>
#include <lemon/concepts/graph.h>

using namespace lemon;

typedef ListDigraph::NodeMap<int> NodeMap;
typedef ListDigraph::ArcMap<int> ArcMap;
typedef ListDigraph::ArcMap<double> ArcCost;

/**********************************************************************************************
 * This class stores the input needed for solving the Routing and Spectrum Allocation problem.
 * This consists of an initial mapping and a set of demands to be routed. A graph associated 
 * to this initial mapping is built as well as an extended graph for each demand to be routed. 
 * \note It uses the LEMON library to build the associated graphs. 
 * *******************************************************************************************/
class RSA{

protected:
    Instance instance;                  /**< An instance describing the initial mapping. **/

    std::vector<Demand> toBeRouted;     /**< The list of demands to be routed in the next optimization. **/

    /** A list of pointers to the extended graph associated with each demand to be routed. 
        \note (*vecGraph[i]) is the graph associated with the i-th demand to be routed. **/
    std::vector< std::shared_ptr<ListDigraph> > vecGraph;

    /** A list of pointers to the ArcMap storing the arc ids of the graph associated with each demand to be routed. 
        \note (*vecArcId[i])[a] is the id of arc a in the graph associated with the i-th demand to be routed. **/
    std::vector< std::shared_ptr<ArcMap> > vecArcId;
    
    /** A list of pointers to the ArcMap storing the arc labels of the graph associated with each demand to be routed. 
        \note (*vecArcLabel[i])[a] is the label of arc a in the graph associated with the i-th demand to be routed. **/
    std::vector< std::shared_ptr<ArcMap> > vecArcLabel;

    /** A list of pointers to the ArcMap storing the arc slices of the graph associated with each demand to be routed. 
        \note (*vecArcSlice[i])[a] is the slice of arc a in the graph associated with the i-th demand to be routed. **/
    std::vector< std::shared_ptr<ArcMap> > vecArcSlice;

    /** A list of pointers to the map storing the arc lengths of the graph associated with each demand to be routed. 
        \note (*vecArcLength[i])[a] is the length of arc a in the graph associated with the i-th demand to be routed. **/
    std::vector< std::shared_ptr<ArcCost> > vecArcLength;

    /** A list of pointers to the NodeMap storing the node ids of the graph associated with each demand to be routed. 
        \note (*vecNodeId[i])[v] is the id of node v in the graph associated with the i-th demand to be routed. **/  
    std::vector< std::shared_ptr<NodeMap> > vecNodeId;

    /** A list of pointers to the NodeMap storing the node labels of the graph associated with each demand to be routed. 
        \note (*vecNodeLabel[i])[v] is the label of node v in the graph associated with the i-th demand to be routed. **/
    std::vector< std::shared_ptr<NodeMap> > vecNodeLabel;
    
    /** A list of pointers to the NodeMap storing the node slices of the graph associated with each demand to be routed. 
        \note (*vecNodeSlice[i])[v] is the slice of node v in the graph associated with the i-th demand to be routed. **/
    std::vector< std::shared_ptr<NodeMap> > vecNodeSlice;

    /** A list of pointers to the ArcMap storing the RSA solution. It stores the id of the demand that is routed through each arc of the graph #i. 
        \note (*vecOnPath[i])[a] is the id the demand routed through arc a in the graph associated with the i-th demand to be routed. **/
    std::vector< std::shared_ptr<ArcMap> > vecOnPath;


    ListDigraph compactGraph;   /**< The simple graph associated with the initial mapping. **/
    ArcMap compactArcId;        /**< ArcMap storing the arc ids of the simple graph associated with the initial mapping. **/
    ArcMap compactArcLabel;     /**< ArcMap storing the arc labels of the simple graph associated with the initial mapping. **/
    ArcCost compactArcLength;   /**< ArcMap storing the arc lengths of the simple graph associated with the initial mapping. **/
    NodeMap compactNodeId;      /**< NodeMap storing the LEMON node ids of the simple graph associated with the initial mapping. **/
    NodeMap compactNodeLabel;   /**< NodeMap storing the node labels of the simple graph associated with the initial mapping. **/

public:
	/****************************************************************************************/
	/*										Constructor										*/
	/****************************************************************************************/

    /** Constructor. A graph associated with the initial mapping (instance) is built as well as an extended graph for each demand to be routed. **/
    RSA(const Instance &instance);

	/****************************************************************************************/
	/*										Getters 										*/
	/****************************************************************************************/
    
    /** Returns the input instance. **/
    Instance getInstance() const{ return instance; }

    /** Returns the vector of demands to be routed. **/
    std::vector<Demand> getToBeRouted() { return toBeRouted; } 
    
    /** Returns the i-th demand to be routed. @param k The index of the required demand. **/
    Demand getToBeRouted_k(int k){ return toBeRouted[k]; }

    /** Returns the number of demands to be routed. **/
    int getNbDemandsToBeRouted() { return toBeRouted.size(); }

    /** Returns the id of a node in a graph. @param n The node. @param d The graph #d. **/
    int getNodeId(const ListDigraph::Node &n, int d) const { return (*vecNodeId[d])[n]; }
    
    /** Returns the label of a node in a graph. @param n The node. @param d The graph #d. **/
    int getNodeLabel(const ListDigraph::Node &n, int d) const { return (*vecNodeLabel[d])[n]; }
    
    /** Returns the slice of a node in a graph. @param n The node. @param d The graph #d. **/
    int getNodeSlice(const ListDigraph::Node &n, int d) const { return (*vecNodeSlice[d])[n]; }
    
    /** Returns the lemon id of an arc in a graph. @param a The arc. @param d The graph #d. **/
    int getArcId(const ListDigraph::Arc &a, int d) const { return (*vecArcId[d])[a]; }
    
    /** Returns the label of an arc in a graph. @param a The arc. @param d The graph #d. **/
    int getArcLabel(const ListDigraph::Arc &a, int d) const { return (*vecArcLabel[d])[a]; }
    
    /** Returns the slice of an arc in a graph. @param a The arc. @param d The graph #d. **/
    int getArcSlice(const ListDigraph::Arc &a, int d) const { return (*vecArcSlice[d])[a]; }
    
    /** Returns the length of an arc in a graph. @param a The arc. @param d The graph #d. **/
    double getArcLength(const ListDigraph::Arc &a, int d) const  {return (*vecArcLength[d])[a]; }

    /** Returns the first node identified by (label, slice) on graph #d. @param d The graph #d. @param label The node's label. @param slice The node's slice. \warning If it does not exist, returns INVALID. **/
    ListDigraph::Node getNode(int d, int label, int slice);

    /** Returns the length of an arc on the compact graph. @param a The arc. */
    double getCompactLength(const ListDigraph::Arc &a) { return compactArcLength[a]; }
    
    /** Returns the coefficient of an arc (according to the chosen metric) on graph #d. @param a The arc. @param d The graph #d. @todo Implement other metrics. For the moment, only 1p is used. **/
    double getCoeff(const ListDigraph::Arc &a, int d);
    
    /** Returns the coefficient of an arc according to metric 1 on graph #d. @param a The arc. @param d The graph #d. **/
    double getCoeffObj1(const ListDigraph::Arc &a, int d);

    /** Returns the coefficient of an arc according to metric 1p on graph #d. @param a The arc. @param d The graph #d. @warning Only adapted for the case of treating one demand at a time. **/
    double getCoeffObj1p(const ListDigraph::Arc &a, int d);

    /** Returns the coefficient of an arc according to metric 2 on graph #d. @param a The arc. @param d The graph #d. @todo I do not really understand this metric -> to be implemented. **/
    double getCoeffObj2(const ListDigraph::Arc &a, int d);

    /** Returns the coefficient of an arc according to metric 4 on graph #d. @param a The arc. @param d The graph #d. **/
    double getCoeffObj4(const ListDigraph::Arc &a, int d);

    /** Returns the coefficient of an arc according to metric 8 on graph #d. @param a The arc. @param d The graph #d. @warning Only adapted for the case of treating one demand at a time. **/
    double getCoeffObj8(const ListDigraph::Arc &a, int d);

	/****************************************************************************************/
	/*										Setters											*/
	/****************************************************************************************/

    /** Changes the vector of demands to be routed. @param vec The new vector of demands. **/
    void setToBeRouted(const std::vector<Demand> &vec){this->toBeRouted = vec;}

    /** Changes the id of a node in a graph. @param n The node. @param d The graph #d. @param val The new id. **/
    void setNodeId(const ListDigraph::Node &n, int d, int val) { (*vecNodeId[d])[n] = val; }
    
    /** Changes the label of a node in a graph. @param n The node. @param d The graph #d. @param val The new label. **/
    void setNodeLabel(const ListDigraph::Node &n, int d, int val) { (*vecNodeLabel[d])[n] = val; }
    
    /** Changes the slice of a node in a graph. @param n The node. @param d The graph #d. @param val The new slice position. **/
    void setNodeSlice(const ListDigraph::Node &n, int d, int val) { (*vecNodeSlice[d])[n] = val; }
    
    /** Changes the id of an arc in a graph. @param a The arc. @param d The graph #d. @param val The new id. **/
    void setArcId(const ListDigraph::Arc &a, int d, int val) { (*vecArcId[d])[a] = val; }
    
    /** Changes the label of an arc in a graph. @param a The arc. @param d The graph #d. @param val The new label. **/
    void setArcLabel(const ListDigraph::Arc &a, int d, int val) { (*vecArcLabel[d])[a] = val; }
    
    /** Changes the slice of an arc in a graph. @param a The arc. @param d The graph #d. @param val The new slice position. **/
    void setArcSlice(const ListDigraph::Arc &a, int d, int val) { (*vecArcSlice[d])[a] = val; }
    
    /** Changes the length of an arc in a graph. @param a The arc. @param d The graph #d. @param val The new length. **/
    void setArcLength(const ListDigraph::Arc &a, int d, double val) { (*vecArcLength[d])[a] = val; }

	/****************************************************************************************/
	/*										Methods											*/
	/****************************************************************************************/

    /** Builds the simple graph associated with the initial mapping. **/
    void buildCompactGraph();

    /** Creates an arc -- and its nodes if necessary -- between nodes (source,slice) and (target,slice) on a graph. @param d The graph #d. @param source The source node's id. @param target The target node's id. @param linkLabel The arc's label. @param slice The arc's slice position. @param l The arc's length. **/
    void addArcs(int d, int source, int target, int linkLabel, int slice, double l);    
    
    /** Updates the mapping stored in the given instance with the results obtained from RSA solution (i.e., vecOnPath). @param i The instance to be updated.*/
    void updateInstance(Instance &i);

    /** Returns the first node with a given label from the graph associated with the d-th demand to be routed. @note If such node does not exist, returns INVALID. @param d The graph #d. @param label The node's label. **/
    ListDigraph::Node getFirstNodeFromLabel(int d, int label);
    
    /** Contract nodes with the same given label from the graph associated with the d-th demand to be routed. @param d The graph #d. @param label The node's label. **/
    void contractNodesFromLabel(int d, int label);

    /** Delete arcs that are known 'a priori' to be unable to route on graph #d. Erase arcs that do not support the demand's load. @param d The graph #d to be inspected. **/
    void eraseNonRoutableArcs(int d);
    
    /** Erases every arc from graph #d having the given slice and returns the number of arcs removed. @param d The graph #d. @param slice The slice to be removed. **/
    int eraseAllArcsFromSlice(int d, int slice);

    /** Runs preprocessing on every extended graph. **/
    void preprocessing();
    
    /** Performs preprocessing based on the arc lengths and returns true if at least one arc is erased. An arc (u,v) can only be part of a solution if the distance from demand source to u, plus the distance from v to demand target plus the arc length is less than or equal to the demand's maximum length. **/
    bool lengthPreprocessing();

    /** Returns the distance of the shortest path from source to target passing through arc a. \note If there exists no st-path, returns +Infinity. @param d The graph #d. @param source The source node.  @param a The arc required to be present. @param target The target node.  **/
    double shortestDistance(int d, ListDigraph::Node &source, ListDigraph::Arc &a, ListDigraph::Node &target);

	/****************************************************************************************/
	/*										Display											*/
	/****************************************************************************************/

    /** Displays the nodes of graph #d that are incident to the node identified by (label,slice).  @param d The graph #d. @param label The node's label. @param slice The node's slice position. **/
    void displayNodesIncidentTo(int d, int label, int slice);

    /** Displays the nodes of graph #d that have a given label. @param d The graph #d. @param label The node's label. **/
    void displayNodesFromLabel(int d, int label);

    /** Displays the paths found for each of the new routed demands. **/
    void displayPaths();
    
    /** Displays an arc from the graph #d. @param d The graph #d. @param a The arc to be displayed. **/
    void displayArc(int d, const ListDigraph::Arc &a);

    /** Displays a node from the graph #d. @param d The graph #d. @param n The node to be displayed. */
    void displayNode(int d, const ListDigraph::Node &n);
    
    /** Display all arcs from the graph #d. @param d The graph #d. **/
    void displayGraph(int d);

    /** Displays the demands to be routed in the next optimization. **/
    void displayToBeRouted();

};
#endif