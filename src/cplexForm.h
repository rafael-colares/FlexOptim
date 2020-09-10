#ifndef __cplexForm__h
#define __cplexForm__h

#include "solver.h"
#include "FormulationComponents.h"

typedef IloArray<IloBoolVarArray> IloBoolVarMatrix;
typedef IloArray<IloNumVarArray> IloNumVarMatrix;
typedef std::vector<std::vector<Variable> > VarMatrix;


/*********************************************************************************************
* This class implements the Online Routing and Spectrum Allocation through a Flow MIP 
* formulation using CPLEX.	
*********************************************************************************************/
class FlowForm : public Solver{

private:
    IloBoolVarMatrix x;				/**< The matrix of assignement variables used in the MIP. x[d][a]=1 if the d-th demand is routed through the arc from index a. **/
    IloIntVarArray maxSlicePerLink;	/**< The array of variables used in the MIP for verifying the max used slice position for each link in the topology network. maxSlicePerLink[i]=p if p is the max used slice position from the link with id i. **/
	IloIntVar maxSliceOverall;		/**< The max used slice position throughout all the network. **/

public:

	/****************************************************************************************/
	/*										Constructors									*/
	/****************************************************************************************/

	/** Constructor. Builds the Online RSA mixed-integer program and solves it using CPLEX.  @param instance The instance to be solved. **/
    FlowForm(const Instance &instance);

	/****************************************************************************************/
	/*										Getters											*/
	/****************************************************************************************/

	/** Returns the objective function expression. **/
    IloExpr getObjFunction(Input::ObjectiveMetric chosenObjective);

	/** Returns the source constraint associated with a demand and a node. @param demand The demand. @param d The demand index. @param nodeLabel The node label. **/
    IloRange getSourceConstraint_d_n(const Demand & demand, int d, int nodeLabel);
	
	/** Returns the flow conservation constraint associated with a demand and a node. @param v The node. @param demand The demand. @param d The demand index. **/
    IloRange getFlowConservationConstraint_i_d(ListDigraph::Node &v, const Demand & demand, int d);

	/** Returns the target constraint associated with a demand. @param demand The demand. @param d The demand index. **/
    IloRange getTargetConstraint_d(const Demand & demand, int d);
	
    /** Returns the length constraint associated with a demand. @param demand The demand. @param d The demand index. **/
    IloRange getLengthConstraint(const Demand &demand, int d);

    /** Returns the non-overlapping constraint associated with an arc and a pair of demands. @param linkLabel The arc label. @param slice The arc slice. @param demand1 The first demand. @param d1 The first demand index. @param demand2 The second demand. @param d2 The second demand index. **/
    IloRange getNonOverlappingConstraint(int linkLabel, int slice, const Demand & demand1, int d1, const Demand & demand2, int d2);
    
	/** Returns the Link's Max Used Slice Position constraint associated with a link and a demand. @param linkIndex The link's index. @param d The demand's index. **/
	IloRange getMaxUsedSlicePerLinkConstraints(int linkIndex, int d);

	/** Returns the Overall Max Used Slice Position constraints associated with a link and a demand. @param linkIndex The link's index. @param d The demand's index. **/
	IloRange getMaxUsedSliceOverallConstraints(int d);
	/** Returns the Overall Max Used Slice Position constraints associated with a link and a demand. @param linkIndex The link's index. **/
	IloRange getMaxUsedSliceOverallConstraints2(int linkIndex, int s);
	/** Returns the Overall Max Used Slice Position constraints associated with a link and a demand. @param linkIndex The link's index. **/
	IloRange getMaxUsedSliceOverallConstraints3(int linkIndex, int s);
	
	/** Returns the constraint from the first set of improved non-overlapping associated with an arc, a demand and a load. @param linkLabel The arc label. @param slice The arc slice. @param min_load The minimum load required for a demand to be part of D[k]. @param demand2 The second demand. @param d2 The second demand index.**/
	IloRange getImprovedNonOverlappingConstraint_1(int linkLabel, int slice, int min_load, const Demand & demand2, int d2);
	
	/** Returns the constraint from the second set of improved non-overlapping associated with an arc and a pair of loads. @param linkLabel The arc label. @param slice The arc slice. **/
	IloRange getImprovedNonOverlappingConstraint_2(int linkLabel, int slice);
	
	/****************************************************************************************/
	/*										Setters											*/
	/****************************************************************************************/

	/** Defines the decision variables need in the MIP formulation. **/
    void setVariables();

	/** Defines the objective function. @param chosenObjective The chosen objective metric. **/
    void setObjective(Input::ObjectiveMetric chosenObjective);

	/** Defines Source constraints. At most one arc leaves each node and exactly one arc leaves the source. **/
    void setSourceConstraints();

	/** Defines Flow Conservation constraints. If an arc enters a node, then an arc must leave it. **/
    void setFlowConservationConstraints();
	
	/** Defines Target constraints. Exactly one arc enters the demand's target. **/
    void setTargetConstraints();
	
	/** Defines Length constraints. Demands must be routed within a length limit. **/
    void setLengthConstraints();

	/** Defines Non-Overlapping constraints. Demands must not overlap eachother's slices. **/
    void setNonOverlappingConstraints();
	
	/** Defines the Link's Max Used Slice Position constraints. The max used slice position on each link must be greater than every slice position used in the link. **/
	void setMaxUsedSlicePerLinkConstraints();

	/** Defines the Overall Max Used Slice Position constraints. The max used slice position overall must be greater than every other slice position used in the network. **/
	void setMaxUsedSliceOverallConstraints();
	/** Defines the Overall Max Used Slice Position constraints. The max used slice position overall must be greater than every other slice position used in the network. **/
	void setMaxUsedSliceOverallConstraints2();
	/** Defines the Overall Max Used Slice Position constraints. The max used slice position overall must be greater than every other slice position used in the network. **/
	void setMaxUsedSliceOverallConstraints3();

	/** Defines the first set of Improved Non-Overlapping constraints. **/
	void setImprovedNonOverlappingConstraints_1();
		
	/** Defines the second set of Improved Non-Overlapping constraints. **/
	void setImprovedNonOverlappingConstraints_2();

	/****************************************************************************************/
	/*										Methods											*/
	/****************************************************************************************/

	/** Recovers the obtained MIP solution and builds a path for each demand on its associated graph from RSA. **/
    void updatePath();

	/****************************************************************************************/
	/*										Display											*/
	/****************************************************************************************/

	/** Displays the obtained paths. **/
    void displayOnPath();

	/** Displays the value of each variable in the obtained solution. **/
    void displayVariableValues();

	/****************************************************************************************/
	/*										Destructor										*/
	/****************************************************************************************/

	/** Destructor. Clears the variable matrices, cplex model and environment. **/
	~FlowForm();

};


#endif