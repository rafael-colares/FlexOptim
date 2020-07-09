#ifndef __YoussoufForm__h
#define __YoussoufForm__h

#include "solver.h"
#include "genericCallback.h"

typedef IloArray<IloBoolVarArray> IloBoolVarMatrix;
typedef IloArray<IloBoolVarMatrix> IloBoolVarMatrix3D;


/*********************************************************************************************
* This class implements the Online Routing and Spectrum Allocation through the Youssouf's 
* Edge-Node MIP formulation using CPLEX.	
*********************************************************************************************/
class YoussoufForm : public Solver{

private:
    IloBoolVarMatrix x;				/**< Demand-Edge variables. x[e][k] equals 1 if demand k is routed through edge e. **/
    IloBoolVarMatrix z;				/**< Demand-Slot variables. z[s][k] equals 1 if slot s is the last slot allocated for demand k. **/
    IloBoolVarMatrix3D t;			/**< Demand-Edge-Slot variables. t[e][s][k] equals 1 if slot s is assigned to demand k on edge e. **/
    IloIntVarArray maxSlicePerLink;	/**< The array of variables used in the MIP for verifying the max used slice position for each link in the topology network. maxSlicePerLink[i]=p if p is the max used slice position from the link with id i. **/
	IloIntVar maxSliceOverall;		/**< The max used slice position throughout all the network. **/

public:

	/****************************************************************************************/
	/*										Constructors									*/
	/****************************************************************************************/

	/** Constructor. Builds the Online RSA mixed-integer program and solves it using CPLEX.  @param instance The instance to be solved. **/
    YoussoufForm(const Instance &instance);

	/****************************************************************************************/
	/*										Getters											*/
	/****************************************************************************************/

	/** Returns the objective function expression. **/
    IloExpr getObjFunction();

	/** Returns the source constraint associated with a demand and a node. @param k The demand index. **/
    IloRange getOriginConstraint_k(int k);
	
	/** Returns the flow conservation constraint associated with a demand and a node. @param var The variable matrix x[d][a]. @param mod The CPLEX model. @param v The node. @param demand The demand. @param d The demand index. **
    IloRange getFlowConservationConstraint_i_d(IloBoolVarMatrix &var, IloModel &mod, ListDigraph::Node &v, const Demand & demand, int d);

	/** Returns the target constraint associated with a demand. @param k The demand index. **/
    IloRange getDestinationConstraint_k(int k);
	
	/** Returns the degree constraint associated with a demand k and a node v. @param k The demand index. @param v The node. **/
	IloRange getDegreeConstraint_k(int k, ListGraph::Node &v);
	
    /** Returns the transmission reach constraint associated with a demand. @param k The demand index. **/
    IloRange getTransmissionReachConstraint_k(int k);
	
	/** Returns the channel-selection constraint associated with a demand. @param k The demand index. **/
	IloRange getChannelSelectionConstraint_k(int k);
	
	/** Returns the Forbidden-Slot constraint associated with a demand. @param k The demand index. **/
	IloRange getForbiddenSlotConstraint_k(int k);

	/** Returns the Edge-Slot constraint associated with a demand and an edge. @param k The demand index. @param e The edge index. **/
	IloRange getEdgeSlotConstraint_k_e(int k, int e);
	
	/** Returns the Demand-Edge-Slot constraint associated with a demand, an edge and a slice. @param k The demand index. @param e The edge index. @param s The slice index. **/
	IloRange getDemandEdgeSlotConstraint_k_e_s(int k, int e, int s);

    /** Returns the non-overlapping constraint associated with an arc and a pair of demands. @param e The edge index. @param s The slice index. **/
    IloRange getNonOverlappingConstraint_e_s(int e, int s);
    
	/** Returns the Link's Max Used Slice Position constraint associated with a link and a demand. @param var The assignment variable matrix. @param maxSlicePerLink The array of variables refering to the max used slice position of each link. @param linkIndex The link's index. @param d The demand's index. @param mod The CPLEX model. **
	IloRange getMaxUsedSlicePerLinkConstraints(IloBoolVarMatrix &var, IloIntVarArray &maxSlicePerLink, int linkIndex, int d, IloModel &mod);

	/** Returns the Overall Max Used Slice Position constraints associated with a link and a demand. @param var The assignment variable matrix. @param maxSliceOverall The variable refering to the max used slice position throughout the whole network. @param linkIndex The link's index. @param d The demand's index. @param mod The CPLEX model.**
	IloRange getMaxUsedSliceOverallConstraints(IloBoolVarMatrix &var, IloIntVar &maxSliceOverall, int linkIndex, int d, IloModel &mod);
	
	/****************************************************************************************/
	/*										Setters											*/
	/****************************************************************************************/

	/** Defines the decision variables needed in the MIP formulation. **/
    void setVariables();

	/** Defines the objective function. **/
    void setObjective();

	/** Defines Origin constraints. At most one arc leaves each node and exactly one arc leaves the source. @param var The variable matrix. @param mod The CPLEX model. **/
    void setOriginConstraints();

	/** Defines Flow Conservation constraints. If an arc enters a node, then an arc must leave it. @param var The variable matrix. @param mod The CPLEX model. **
    void setFlowConservationConstraints(IloBoolVarMatrix &var, IloModel &mod);
	
	/** Defines Destination constraints. **/
    void setDestinationConstraints();
	
	/** Defines Degree constraints. At most two edges are adjacent to any node. **/
	void setDegreeConstraints();

	/** Defines Transmission-Reach constraints. Demands must be routed within a length limit. **/
    void setTransmissionReachConstraints();

	/** Defines Channel-Selection constraints. Only one slot can be the last assigned to a demand. **/
	void setChannelSelectionConstraints();

	/** Defines Forbidden-Slot constraints. **/
	void setForbiddenSlotConstraints();
	
	/** Defines Edge-Slot constraints. **/
	void setEdgeSlotConstraints();
	
	/** Defines Demand-Edge-Slot constraints. **/
	void setDemandEdgeSlotConstraints();

	/** Defines Non-Overlapping constraints. Demands must not overlap eachother's slices. **/
    void setNonOverlappingConstraints();
	
	/** Defines the Link's Max Used Slice Position constraints. The max used slice position on each link must be greater than every slice position used in the link. @param var The assignment variable matrix. @param maxSlicePerLink The array of variables refering to the max used slice position of each link. @param mod The CPLEX model.**
	void setMaxUsedSlicePerLinkConstraints(IloBoolVarMatrix &var, IloIntVarArray &maxSlicePerLink, IloModel &mod);

	/** Defines the Overall Max Used Slice Position constraints. The max used slice position overall must be greater than every other slice position used in the network. @param var The assignment variable matrix. @param maxSliceOverall The variable refering to the max used slice position throughout the whole network. @param mod The CPLEX model.**
	void setMaxUsedSliceOverallConstraints(IloBoolVarMatrix &var, IloIntVar maxSliceOverall, IloModel &mod);

	/****************************************************************************************/
	/*										Methods											*/
	/****************************************************************************************/

	/** Recovers the obtained MIP solution and builds a path for each demand on its associated graph from RSA. **/
    void updatePath();

	/****************************************************************************************/
	/*										Display											*/
	/****************************************************************************************/

	/** Displays the obtained paths. **
    void displayOnPath();

	/** Displays the value of each variable in the obtained solution. **
    void displayVariableValues();

	/****************************************************************************************/
	/*										Destructor										*/
	/****************************************************************************************/

	/** Destructor. Clears the variable matrices, cplex model and environment. **/
	~YoussoufForm();

};


#endif