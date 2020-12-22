#ifndef LAG_NON_OVERLAPPING_H
#define LAG_NON_OVERLAPPING_H

#include "AbstractLagrangianFormulation.h"
#include <set>

using namespace lemon;

class lagNonOverlapping: public AbstractLagFormulation{

    private:

        /********************************* MULTIPLIERS ***********************************/

        /** A vector storing the value of the Lagrangian multipliers associated with Length Constraints. **/
        std::vector<double> lagrangianMultiplierLength;

        /** A 2-dimensional vector storing the value of the Lagrangian multipliers associated with each Source/Target constraint. **/
        std::vector<std::vector<double>> lagrangianMultiplierSourceTarget;

        /** A 2-dimensional vector storing the value of the Lagrangian multipliers associated with each Flow constraint. **/
        std::vector<std::vector<double>> lagrangianMultiplierFlow;

        /********************************** SLACK ***************************************/

        /** Stores the value of the slack of lengths constraints (i.e., b - Dx). */
        std::vector<double> lengthSlack;

        /** Stores the value of the slack of Source/Target constraints (i.e., b - Dx). */
        std::vector<std::vector<double>> sourceTargetSlack;

        /** Stores the value of the slack of flow constraints (i.e., b - Dx). */
        std::vector<std::vector<double>> flowSlack;

        /***************************** ASSIGNMENT MATRIX *******************************/

        std::vector< std::vector<bool> > assignmentMatrix;

        /***************************** BUILDING AUXILIARY GRAPH *******************************/

        /** A list of pointers to the extended graph associated with each EDGE to be routed. 
        \note (*vecGraph[i]) is the graph associated with the i-th EDGE to be routed. **/
        std::vector<std::shared_ptr<ListDigraph>> vecEGraph;

        /** A list of pointers to the NodeMap storing the node LEMON IDs of the graph associated with each EDGE to be routed. 
        \note (*vecENodeID[i])[v] is the LEMON ID of node v in the graph associated with the i-th EDGE to be routed. **/
        /* using in the assigment matrix */
        std::vector<std::shared_ptr<NodeMap>> vecENodeID;

        /** A list of pointers to the NodeMap storing the node INDEXs of the graph associated with each EDGE to be routed. 
        \note (*vecENodeID[i])[v] is the LEMON ID of node v in the graph associated with the i-th EDGE to be routed. **/
        std::vector<std::shared_ptr<NodeMap>> vecENodeIndex;

        /** A list of pointers to the NodeMap storing the node DEMANDS of the graph associated with EDGE to be routed
        \note (*vecENodeDemand[i])[v] is the DEMAND of node v in the graph associated with the i-th EDGE to be routed. **/
        /** Indeed, here each node refers to an arc in the graph for each d (G(D)), so, its the DEMAND of the original arc **/
        std::vector<std::shared_ptr<NodeMap>> vecENodeDemand;

        /** A list of pointers to the NodeMap storing the node SLICES of the graph associated with EDGE to be routed
        \note (*vecENodeSlice[i])[v] is the SLICE of node v in the graph associated with the i-th EDGE to be routed. **/
        /** Indeed, here each node refers to an arc in the graph for each d (G(D)), so, its the SLICE of the original arc **/
        std::vector<std::shared_ptr<NodeMap>> vecENodeSlice;

        /** A list of pointers to the NodeMap storing the node FIRST CONSTRAINT of the graph associated with EDGE to be routed
        \note (*vecENodeFirstConst[i])[v] is the INDEX OF THE FIRST CONSTRAINT of node v in the graph associated with the i-th 
        EDGE to be routed. **/
        /** Indeed, here each node refers to an arc in the graph for each d (G(D)), so, its the FIRST CONSTRAINT of the original arc **/
        std::vector<std::shared_ptr<NodeMap>> vecENodeFirstConst;

        /** A list of pointers to the ArcMap storing the arc ids of the graph associated with each demand to be routed. 
        \note (*vecArcId[i])[a] is the id of arc a in the graph associated with the i-th demand to be routed. **/
        std::vector<std::shared_ptr<ArcMap>> vecEArcId;

        /** A list of pointers to the map storing the arc COSTS of the graph associated with each EDGE to be routed. 
        \note (*vecECost[i])[a] is the COST of arc a in the graph associated with the i-th EDGE to be routed. **/
        /** The cost is associated with the destination node of a, this node represents an arc in the  graph for each d (G(D)) */
        std::vector<std::shared_ptr<ArcCost>> vecECost;

        /* Refers to the cost of an arc during iteration k. cost = c_{ij} + u_k*length_{ij} */
        /* This is the cost of the original variables, corresponding to the original graph G(D)*/
        std::vector<std::shared_ptr<ArcCost>> cost; 

        /** A list of pointers to the map storing the node ORIGINAL ARCS of the graph associated with each EDGE to be routed. 
        \note (*vecENodeArc[i])[a] is the ARC of nde a in the graph associated with the i-th EDGE to be routed. **/
        /** The ARC is in the  graph for each d (G(D)) */
        std::vector<std::shared_ptr<ListDigraph::NodeMap<std::shared_ptr<ListDigraph::Arc>>>> vecENodeArc;

        /* Vector of the id of the node source for each edge*/
        std::vector<int> vecESourceIndex;

        /* Vector of the id of the destination node for each edge*/
        std::vector<int> vecEDestinationIndex;

    public:
        /* *******************************************************************************
        *                             INITIALIZATION METHODS
        ******************************************************************************* */
        lagNonOverlapping(const Instance &instance):AbstractLagFormulation(instance){}

        /** Sets all initial parameters **/
        void init();

        /********************************* MULTIPLIERS ***********************************/

        /** Sets the initial lagrangian multipliers values for the subgradient to run. **/
        void initMultipliers();

        /** Sets the initial lagrangian multipliers associated with length constraints. **/
        void initializeLengthMultipliers();

        /** Sets the initial lagrangian multipliers associated with Source/Target constraints **/
        void initializeSourceTargetMultipliers();
        
        /** Sets the initial lagrangian multipliers associated with flow constraints **/
        void initializeFlowMultipliers();

        /********************************** SLACK ***************************************/

        /** Initializes the slack of relaxed constraints. **/
        void initSlacks();

        /** Initializes the slack of length constraints. **/
        void initializeLengthSlacks();
        
        /** Initializes the slack of Source/Target constraints. **/
        void initializeSourceTargetSlacks();

        /** Initializes the slack of Flow constraints. **/
        void initializeFlowSlacks();

        /******************************** COSTS *********************************/

        /** Initializes the costs in the objective function. **/
        void initCosts();

        /******************************** ASSIGMENT MATRIX *********************************/

        /** Initializes the assignement matrix. **/
        void initAssignmentMatrix();

        /***************************** BUILDING AUXILIARY GRAPH *******************************/

        /** Constructs the auxiliary graphs for each edge e **/
        void build_Graphs_e();

        /** Add a node in the auxiliary graph **/
        void addENode(int,int,int,int, const ListDigraph::Arc &);

        /* *******************************************************************************
        *                             RUNING METHODS
        ******************************************************************************* */
        void run();

        double getRealCostFromPath(int, BellmanFord< ListDigraph, ListDigraph::ArcMap<double> > &, const ListDigraph::Node &, const ListDigraph::Node &);

        void updateAssignment_e(int, BellmanFord< ListDigraph, ListDigraph::ArcMap<double> > &, const ListDigraph::Node &, const ListDigraph::Node &);

        void subtractConstantValuesFromLagrCost();

        /** Checks with the slacks if the solution is feasible. **/
        bool checkFeasibility();

        bool checkLengthFeasibility();

        bool checkSourceTargetFeasibility();

        bool checkFlowFeasibility();
 
        /****************************************************************************************/
        /*										Getters 										*/
        /****************************************************************************************/

        /** Returns the multiplier for the length constraint k **/
        double getLengthMultiplier_k(int k) const { return lagrangianMultiplierLength[k]; }

        /** Returns the multiplier for the source target constraint k,v **/
        double getSourceTargetMultiplier_k(int k, int v) const { return lagrangianMultiplierSourceTarget[k][v]; }

        /** Returns the multiplier for the flow constraint k,v **/
        double getFlowMultiplier_k(int k, int v) const { return lagrangianMultiplierFlow[k][v]; }

        /** Return the slack for the the length constraint k **/
        double getLengthSlack_k(int k) const { return lengthSlack[k];}

        /** Return the slack for the the length constraint k **/
        double getSourceTargetSlack_k(int k, int v) const { return sourceTargetSlack[k][v];}

        /** Return the slack for the the length constraint k **/
        double getFlowSlack_k(int k,int v) const { return flowSlack[k][v];}
        
        /** Returns the LEMON id of a node in a graph. @param n The node. @param label The graph index. **/
        int getNodeEId(const ListDigraph::Node &n, int label) const { return (*vecENodeID[label])[n]; }
        
        /** Returns the demand of a node in a graph. @param n The node. @param label The graph index. **/
        int getNodeEDemand(const ListDigraph::Node &n, int label) const { return (*vecENodeDemand[label])[n]; }

        /** Returns the slice of a node in a graph. @param n The node. @param label The graph index. **/
        int getNodeESlice(const ListDigraph::Node &n, int label) const { return (*vecENodeSlice[label])[n]; }

        /** Returns the firstconst of a node in a graph. @param n The node. @param label The graph index. **/
        int getNodeEFirstConst(const ListDigraph::Node &n, int label) const{ return (*vecENodeFirstConst[label])[n]; }

        /** Returns the index of a node in a graph. @param n The node. @param label The graph index. **/
        int getNodeEIndex(const ListDigraph::Node &n, int label) const { return (*vecENodeIndex[label])[n]; }

        /** Returns the LEMON id of an arc in a graph. @param a The arc. @param label The graph index. **/
        int getArcEId(const ListDigraph::Arc &a, int label) const { return (*vecEArcId[label])[a]; }

        /** Returns the cost of an arc in a graph. @param a The arc. @param label The graph index. **/
        double getArcECost(const ListDigraph::Arc &a, int label) const  {return (*vecECost[label])[a]; } 

        /** Returns the original ARC of an node in a graph. @param a The arc. @param label The graph index. **/
        const ListDigraph::Arc & getNodeEArc(const ListDigraph::Node &n, int label) const  {return (*(*vecENodeArc[label])[n]); }

        /** Changes the cost of an arc in a graph. @param a The arc. @param d The graph index. **/
	    double getArcCost(const ListDigraph::Arc &a, int d) const { return (*cost[d])[a]; }

        /** Returns the source index of the graph corresponding to label**/
        int getIndexSource(int label){ return vecESourceIndex[label]; }

        /** Returns the destination of the graph corresponding to label**/
        int getIndexDestination(int label){ return vecEDestinationIndex[label]; }

        /** Returns the corresponding node for the heuristic auxiliary graph. @ param v The node from original graph d. @param d The graph index**/
        //ListDigraph::Node getNodeRef(const ListDigraph::Node &v, int d) const { return (*vecNodeRef[d])[v]; }

        ListDigraph::Node getNodeFromIndex(int, int);

        /** Returns the constraints slack module **/
        double getSlackModule();


        /****************************************************************************************/
        /*										Setters											*/
        /****************************************************************************************/

        void setLengthMultiplier_k (int k, double val) { lagrangianMultiplierLength[k] = val; }
        void setSourceTargetMultiplier_k (int k, int v, double val) { lagrangianMultiplierSourceTarget[k][v] = val; }
        void setFlowMultiplier_k (int k, int v, double val) { lagrangianMultiplierFlow[k][v] = val; }

        void setLengthSlack_k (int k, double val) { lengthSlack[k] = val; }
        void setSourceTargetSlack_k (int k, int v, double val) { sourceTargetSlack[k][v] = val; }
        void setFlowSlack_k (int k, int v, double val) { flowSlack[k][v] = val; }
    

        /** Changes the id of a node in a graph. @param n The node. @param label The graph index. @param val The new id. **/
        void setNodeEId(const ListDigraph::Node &n, int label, int val){ (*vecENodeID[label])[n] = val; }

        /** Changes the demand of a node in a graph. @param n The node. @param label The graph index. @param val The new demand. **/
        void setNodeEDemand(const ListDigraph::Node &n, int label, int val){ (*vecENodeDemand[label])[n] = val; }

        /** Changes the slice of a node in a graph. @param n The node. @param label The graph index. @param val The new slice. **/
        void setNodeESlice(const ListDigraph::Node &n, int label, int val){ (*vecENodeSlice[label])[n] = val; }

        /** Changes the firstconst of a node in a graph. @param n The node. @param label The graph index. @param val The new firstconst. **/
        void setNodeEFirstConst(const ListDigraph::Node &n, int label, int val){ (*vecENodeFirstConst[label])[n] = val; }

        /** Changes the id of an arc in a graph. @param n The node. @param label The graph index. @param val The new id. **/
        void setArcEId(const ListDigraph::Arc &a, int label, int val) { (*vecEArcId[label])[a] = val; }

        /** Changes the cost of an arc in a graph. @param a The arc. @param label The graph index. @param val The new length. **/
        void setArcECost(const ListDigraph::Arc &a, int label, double val) { (*vecECost[label])[a] = val; }

        /** Changes the cost of an arc in a graph. @param a The arc. @param d The graph index. @param val The new cost value. **/
	    void setArcCost(const ListDigraph::Arc &a, int d, double val) { (*cost[d])[a] = val; }

	    /** Increments the cost of an arc in a graph. @param a The arc. @param d The graph index. @param val The value to be added to cost. **/
	    void incArcCost(const ListDigraph::Arc &a, int d, double val) { (*cost[d])[a] += val; }

        /** Changes the firstconst of a node in a graph. @param n The node. @param label The graph index. @param val The new firstconst. **/
        void setNodeEArc(const ListDigraph::Node &n, int label, const ListDigraph::Arc &a){ (*vecENodeArc[label])[n] =  std::make_shared<ListDigraph::Arc>(a); }

        /** Changes the Index of a node in a graph. @param n The node. @param label The graph index. @param val The new index. **/
        void setNodeEIndex(const ListDigraph::Node &n, int label, int val){ (*vecENodeIndex[label])[n] = val; }

        /** Changes the source index of the graph corresponding to label**/
        void setIndexSource(int id, int label){ vecESourceIndex[label] = id; }

        /** Changes the destination index of the graph corresponding to label**/
        void setIndexDestination(int id, int label){ vecEDestinationIndex[label] = id; }

        /** Updates the slack of a lehgth constraint using the assigment matrix **/
        void updateLengthSlack();

        /** Updates the slack of a source/target constraint using the assigment matrix **/
        void updateSourceTargetSlack();

        /** Updates the slack of a flow constraint using the assigment matrix **/
        void updateFlowSlack();

        /** Updates the arc costs. @note cost = coeff + u_d*length **/
        void updateCosts();

        /* Updates lagrangian multipliers with the rule: u[k+1] = u[k] + t[k]*violation */
        void updateMultiplier(double);

        void updateLengthMultiplier(double);

        void updateSourceTargetMultiplier(double);

        void updateFlowMultiplier(double);

        void displayEArc(int, const ListDigraph::Arc &);

        void displayEGraph(int);

        void displayAssignmentMatrix(int);

        void displayENode(const ListDigraph::Node &, int);
        void displaySlack(std::ostream & = std::cout);

        void displayMultiplier(std::ostream & = std::cout);
        void createGraphFile(int);

        /* *******************************************************************************
        *                             DESTRUCTOR
        ******************************************************************************* */

        ~lagNonOverlapping();
};




#endif