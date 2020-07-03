#ifndef __cplexForm__h
#define __cplexForm__h

#include "solver.h"

typedef IloArray<IloBoolVarArray> IloBoolVarMatrix;


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

	/** Returns the objective function expression. @param var The variable matrix x[d][a]. @param maxSliceFromLink The array of variables refering to the max used slice position of each link. @param maxSlice The variable storing the max used slice position in the whole network. @param mod The CPLEX model. **/
    IloExpr getObjFunction(IloBoolVarMatrix &var, IloIntVarArray &maxSliceFromLink, IloIntVar &maxSlice, IloModel &mod);

	/** Returns the source constraint associated with a demand and a node. @param var The variable matrix x[d][a]. @param mod The CPLEX model. @param demand The demand. @param d The demand index. @param label The node label. **/
    IloRange getSourceConstraint_d_n(IloBoolVarMatrix &var, IloModel &mod, const Demand & demand, int d, int label);
	
	/** Returns the flow conservation constraint associated with a demand and a node. @param var The variable matrix x[d][a]. @param mod The CPLEX model. @param v The node. @param demand The demand. @param d The demand index. **/
    IloRange getFlowConservationConstraint_i_d(IloBoolVarMatrix &var, IloModel &mod, ListDigraph::Node &v, const Demand & demand, int d);

	/** Returns the target constraint associated with a demand. @param var The variable matrix x[d][a]. @param mod The CPLEX model. @param demand The demand. @param d The demand index. **/
    IloRange getTargetConstraint_d(IloBoolVarMatrix &var, IloModel &mod, const Demand & demand, int d);
	
    /** Returns the length constraint associated with a demand. @param var The variable matrix x[d][a]. @param mod The CPLEX model. @param demand The demand. @param d The demand index. **/
    IloRange getLengthConstraint(IloBoolVarMatrix &var, IloModel &mod, const Demand &demand, int d);

    /** Returns the non-overlapping constraint associated with an arc and a pair of demands. @param var The variable matrix x[d][a]. @param mod The CPLEX model. @param linkLabel The arc label. @param slice The arc slice. @param demand1 The first demand. @param d1 The first demand index. @param demand2 The second demand. @param d2 The second demand index. **/
    IloRange getNonOverlappingConstraint(IloBoolVarMatrix &var, IloModel &mod, int linkLabel, int slice, const Demand & demand1, int d1, const Demand & demand2, int d2);
    
	/** Returns the Link's Max Used Slice Position constraint associated with a link and a demand. @param var The assignment variable matrix. @param maxSlicePerLink The array of variables refering to the max used slice position of each link. @param linkIndex The link's index. @param d The demand's index. @param mod The CPLEX model. **/
	IloRange getMaxUsedSlicePerLinkConstraints(IloBoolVarMatrix &var, IloIntVarArray &maxSlicePerLink, int linkIndex, int d, IloModel &mod);

	/** Returns the Overall Max Used Slice Position constraints associated with a link and a demand. @param var The assignment variable matrix. @param maxSliceOverall The variable refering to the max used slice position throughout the whole network. @param linkIndex The link's index. @param d The demand's index. @param mod The CPLEX model.**/
	IloRange getMaxUsedSliceOverallConstraints(IloBoolVarMatrix &var, IloIntVar &maxSliceOverall, int linkIndex, int d, IloModel &mod);
	
	/** Returns the constraint from the first set of improved non-overlapping associated with an arc, a demand and a load.  @param var The variable matrix x[d][a]. @param mod The CPLEX model. @param linkLabel The arc label. @param slice The arc slice. @param min_load The minimum load required for a demand to be part of D[k]. @param demand2 The second demand. @param d2 The second demand index.**/
	IloRange getImprovedNonOverlappingConstraint_1(IloBoolVarMatrix &var, IloModel &mod, int linkLabel, int slice, int min_load, const Demand & demand2, int d2);
	
	/** Returns the constraint from the second set of improved non-overlapping associated with an arc and a pair of loads.  @param var The variable matrix x[d][a]. @param mod The CPLEX model. @param linkLabel The arc label. @param slice The arc slice. @param min_load1 The minimum load required for a demand to be part of D[k1]). @param min_load2 The minimum load required for a demand to be part of D[k2]). **/
	IloRange getImprovedNonOverlappingConstraint_2(IloBoolVarMatrix &var, IloModel &mod, int linkLabel, int slice, int min_load1, int min_load2);
	
	/****************************************************************************************/
	/*										Setters											*/
	/****************************************************************************************/

	/** Defines the decision variables need in the MIP formulation. @param var The assignement variable matrix. @param maxSliceFromLink The array of variables refering to the max used slice position of each link. @param maxSlice The variable storing the max used slice position in the whole network. @param mod The CPLEX model. **/
    void setVariables(IloBoolVarMatrix &var, IloIntVarArray &maxSliceFromLink, IloIntVar &maxSlice, IloModel &mod);

	/** Defines the objective function. @param var The variable matrix. @param maxSliceFromLink The array of variables refering to the max used slice position of each link. @param maxSlice The variable storing the max used slice position in the whole network. @param mod The CPLEX model. **/
    void setObjective(IloBoolVarMatrix &var, IloIntVarArray &maxSliceFromLink, IloIntVar &maxSlice, IloModel &mod);

	/** Defines Source constraints. At most one arc leaves each node and exactly one arc leaves the source. @param var The variable matrix. @param mod The CPLEX model. **/
    void setSourceConstraints(IloBoolVarMatrix &var, IloModel &mod);

	/** Defines Flow Conservation constraints. If an arc enters a node, then an arc must leave it. @param var The variable matrix. @param mod The CPLEX model. **/
    void setFlowConservationConstraints(IloBoolVarMatrix &var, IloModel &mod);
	
	/** Defines Target constraints. Exactly one arc enters the demand's target. @param var The variable matrix. @param mod The CPLEX model. **/
    void setTargetConstraints(IloBoolVarMatrix &var, IloModel &mod);
	
	/** Defines Length constraints. Demands must be routed within a length limit. @param var The variable matrix. @param mod The CPLEX model. **/
    void setLengthConstraints(IloBoolVarMatrix &var, IloModel &mod);

	/** Defines Non-Overlapping constraints. Demands must not overlap eachother's slices. @param var The variable matrix. @param mod The CPLEX model. **/
    void setNonOverlappingConstraints(IloBoolVarMatrix &var, IloModel &mod);
	
	/** Defines the Link's Max Used Slice Position constraints. The max used slice position on each link must be greater than every slice position used in the link. @param var The assignment variable matrix. @param maxSlicePerLink The array of variables refering to the max used slice position of each link. @param mod The CPLEX model.**/
	void setMaxUsedSlicePerLinkConstraints(IloBoolVarMatrix &var, IloIntVarArray &maxSlicePerLink, IloModel &mod);

	/** Defines the Overall Max Used Slice Position constraints. The max used slice position overall must be greater than every other slice position used in the network. @param var The assignment variable matrix. @param maxSliceOverall The variable refering to the max used slice position throughout the whole network. @param mod The CPLEX model.**/
	void setMaxUsedSliceOverallConstraints(IloBoolVarMatrix &var, IloIntVar maxSliceOverall, IloModel &mod);

	/** Defines the first set of Improved Non-Overlapping constraints. @param var The variable matrix. @param mod The CPLEX model. **/
	void setImprovedNonOverlappingConstraints_1(IloBoolVarMatrix &var, IloModel &mod);
		
	/** Defines the second set of Improved Non-Overlapping constraints. @param var The variable matrix. @param mod The CPLEX model. **/
	void setImprovedNonOverlappingConstraints_2(IloBoolVarMatrix &var, IloModel &mod);

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