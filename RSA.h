#ifndef __RSA__h
#define __RSA__h

#include "Instance.h"
#include "ExtendedGraph.h"

#include<ilcplex/ilocplex.h>
#include <lemon/list_graph.h>
#include <lemon/concepts/graph.h>

using namespace lemon;

typedef ListDigraph::NodeMap<int> NodeMap;
typedef ListDigraph::ArcMap<int> ArcMap;
typedef ListDigraph::ArcMap<double> ArcCost;

class RSA{

protected:
    /* An instance describing the initial mapping. */
    Instance instance;

    /* The list of demands to be routed in the next optimization. */
    std::vector<Demand> toBeRouted; 

    /* A list of pointers to the extended graph associated with each demand to be routed. 
        Ex.: (*vecGraph[i]) is the graph associated with the i-th demand to be routed. */
    std::vector< std::shared_ptr<ListDigraph> > vecGraph;

    /* A list of pointers to the ArcMap storing the arc ids of the graph associated with each demand to be routed. 
        Ex.: (*vecArcId[i])[a] is the id of arc a in the graph associated with the i-th demand to be routed. */
    std::vector< std::shared_ptr<ArcMap> > vecArcId;
    
    /* A list of pointers to the ArcMap storing the arc labels of the graph associated with each demand to be routed. 
        Ex.: (*vecArcLabel[i])[a] is the label of arc a in the graph associated with the i-th demand to be routed. */
    std::vector< std::shared_ptr<ArcMap> > vecArcLabel;

    /* A list of pointers to the ArcMap storing the arc slices of the graph associated with each demand to be routed. 
        Ex.: (*vecArcSlice[i])[a] is the slice of arc a in the graph associated with the i-th demand to be routed. */
    std::vector< std::shared_ptr<ArcMap> > vecArcSlice;

    /* A list of pointers to the map storing the arc lengths of the graph associated with each demand to be routed. 
        Ex.: (*vecArcLength[i])[a] is the length of arc a in the graph associated with the i-th demand to be routed. */
    std::vector< std::shared_ptr<ArcCost> > vecArcLength;

    /* A list of pointers to the NodeMap storing the node ids of the graph associated with each demand to be routed. 
        Ex.: (*vecNodeId[i])[v] is the id of node v in the graph associated with the i-th demand to be routed. */  
    std::vector< std::shared_ptr<NodeMap> > vecNodeId;

    /* A list of pointers to the NodeMap storing the node labels of the graph associated with each demand to be routed. 
        Ex.: (*vecNodeLabel[i])[v] is the label of node v in the graph associated with the i-th demand to be routed. */
    std::vector< std::shared_ptr<NodeMap> > vecNodeLabel;
    
    /* A list of pointers to the NodeMap storing the node slices of the graph associated with each demand to be routed. 
        Ex.: (*vecNodeSlice[i])[v] is the slice of node v in the graph associated with the i-th demand to be routed. */
    std::vector< std::shared_ptr<NodeMap> > vecNodeSlice;

    /* A list of pointers to the ArcMap storing the demand ids that are routed through each arc of the graph associated with each demand to be routed. 
        Ex.: (*vecOnPath[i])[a] is the id the demand routed through arc a in the graph associated with the i-th demand to be routed. */
    std::vector< std::shared_ptr<ArcMap> > vecOnPath;


    ListDigraph compactGraph;
    ArcMap compactArcId;        // provides the lemonId of an arc
    ArcMap compactArcLabel;     // provides the id e of an arc (e,s)
    ArcCost compactArcLength;      // provides the length of an arc
    NodeMap compactNodeId;       // provides the lemonId of a node
    NodeMap compactNodeLabel;    // provides the id n of a node (n,s)

public:
	/************************************************/
	/*				Constructor						*/
	/************************************************/

    /* Build the extended graphs from an initial mapping (instance). */
    RSA(const Instance &instance);

	/************************************************/
	/*				    Getters						*/
	/************************************************/

    Instance getInstance() const{ return instance; }
    std::vector<Demand> getToBeRouted() { return toBeRouted; } 
    Demand getToBeRouted_k(int k){ return toBeRouted[k]; }
    int getNbDemandsToBeRouted() { return toBeRouted.size(); }

    int getNodeId(const ListDigraph::Node &n, int d) const { return (*vecNodeId[d])[n]; }
    int getNodeLabel(const ListDigraph::Node &n, int d) const { return (*vecNodeLabel[d])[n]; }
    int getNodeSlice(const ListDigraph::Node &n, int d) const { return (*vecNodeSlice[d])[n]; }
    int getArcId(const ListDigraph::Arc &a, int d) const { return (*vecArcId[d])[a]; }
    int getArcLabel(const ListDigraph::Arc &a, int d) const { return (*vecArcLabel[d])[a]; }
    int getArcSlice(const ListDigraph::Arc &a, int d) const { return (*vecArcSlice[d])[a]; }
    double getArcLength(const ListDigraph::Arc &a, int d) const  {return (*vecArcLength[d])[a]; }

    /* For graph #d, returns the node with id (label, slice). If it does not exist, returns INVALID. */
    ListDigraph::Node getNode(int d, int label, int slice);

    /* Returns the length of arc a on the compact graph. */
    double getCompactLength(const ListDigraph::Arc &a) { return compactArcLength[a]; }
    
    /* Returns the coefficient of arc a (according to the chosen metric) on graph #d. */
    double getCoeff(const ListDigraph::Arc &a, int d);
    
	/************************************************/
	/*				    Setters						*/
	/************************************************/
    void setToBeRouted(const std::vector<Demand> &vec){this->toBeRouted = vec;}

    void setNodeId(const ListDigraph::Node &n, int d, int val) { (*vecNodeId[d])[n] = val; }
    void setNodeLabel(const ListDigraph::Node &n, int d, int val) { (*vecNodeLabel[d])[n] = val; }
    void setNodeSlice(const ListDigraph::Node &n, int d, int val) { (*vecNodeSlice[d])[n] = val; }
    void setArcId(const ListDigraph::Arc &a, int d, int val) { (*vecArcId[d])[a] = val; }
    void setArcLabel(const ListDigraph::Arc &a, int d, int val) { (*vecArcLabel[d])[a] = val; }
    void setArcSlice(const ListDigraph::Arc &a, int d, int val) { (*vecArcSlice[d])[a] = val; }
    void setArcLength(const ListDigraph::Arc &a, int d, double val) { (*vecArcLength[d])[a] = val; }

	/************************************************/
	/*				    Display						*/
	/************************************************/

    /* For graph #d, display the nodes that are incident to node (label,slice). */
    void displayNodesIncidentTo(int d, int label, int slice);

    /* For graph #d, display the nodes that have a given label. */
    void displayNodesFromLabel(int d, int label);

    /* Display the paths found for each of the new routed demands. */
    void displayPaths();
    
    /* Display arc a from the graph associated with the d-th demand. */
    void displayArc(int d, const ListDigraph::Arc &a);

    /* Display node n from the graph associated with the d-th demand. */
    void displayNode(int d, const ListDigraph::Node &n);

    /* Displays the demands to be routed in the next optimization. */
    void displayToBeRouted();

	/************************************************/
	/*				    Methods						*/
	/************************************************/

    /* On graph #d, creates an arc -- and its nodes if necessary -- between nodes (linkSourceLabel,slice) and (linkTargetLabel,slice) */
    void addArcs(int d, int linkSourceLabel, int linkTargetLabel, int linkLabel, int slice, double l);    
    
    /* Updates the mapping stored in instance i with the results obtained from RSA solution.*/
    void updateInstance(Instance &i);

    /* Returns the first node with a given label from the graph associated with the d-th demand. If such node does not exist, return INVALID. */
    ListDigraph::Node getFirstNodeFromLabel(int d, int label);
    
    /* Contract nodes with the same given label from the graph associated with the d-th demand. */
    void contractNodesFromLabel(int d, int label);

    /* Delete arcs that are known a priori to be unable to route on graph #d. */
    void eraseNonRoutableArcs(int d);
    
};
#endif