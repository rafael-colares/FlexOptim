#ifndef __solver__h
#define __solver__h

#include "RSA.h"


//typedef IloArray<IloNumVarArray> IloNumVarMatrix;
typedef IloArray<IloBoolVarArray> IloBoolVarMatrix;


/**********************************************************************************************
 * This class specifies which MIP solver to be used and stores the new arc maps used for 
 * identifying arcs after preprocessing has erased some arcs. The methods needed for solving 
 * the Routing and Spectrum Allocation problem via CPLEX. 
 * \note It uses the LEMON library to build the arc maps and CPLEX concert library to build 
 * expressions and constraints. 
 * *******************************************************************************************/
class Solver : public RSA{
protected:

public:
	/****************************************************************************************/
	/*										Constructor										*/
	/****************************************************************************************/

	/** Constructor. The RSA constructor is called and the arc map storing the index of the preprocessed graphs associated is built. @param inst The instance to be solved. **/
    Solver(const Instance &inst);

	/****************************************************************************************/
	/*										 Getters										*/
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

	/****************************************************************************************/
	/*										 Setters										*/
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

};    
#endif