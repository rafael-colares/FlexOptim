#ifndef __solver__h
#define __solver__h

#include "RSA.h"


typedef IloArray<IloNumVarArray> IloNumVarMatrix;


/**********************************************************************************************
 * This class specifies which MIP solver to be used and stores the new arc maps used for 
 * identifying arcs after preprocessing has erased some arcs. The methods needed for solving 
 * the Routing and Spectrum Allocation problem via CPLEX. 
 * \note It uses the LEMON library to build the arc maps and CPLEX concert library to build 
 * expressions and constraints. 
 * *******************************************************************************************/
class Solver : public RSA{
protected:
	/** A list of pointers to the ArcMap storing the arc index of the preprocessed graph associated with each demand to be routed. 
        \note (*vecArcIndex[i])[a] is the index of the arc a in the preprocessed graph associated with the i-th demand to be routed. **/
	std::vector< std::shared_ptr<ArcMap> > vecArcIndex; 

public:
	/****************************************************************************************/
	/*										Constructor										*/
	/****************************************************************************************/

	/** Constructor. The RSA constructor is called and the arc map storing the index of the preprocessed graphs associated is built. @param inst The instance to be solved. **/
    Solver(const Instance &inst);

	/****************************************************************************************/
	/*										 Getters										*/
	/****************************************************************************************/

	/** Returns the index of an arc in a preprocessed graph. @param a The arc. @param d The graph #d. **/
    int getArcIndex(const ListDigraph::Arc &a, int d) const { return (*vecArcIndex[d])[a]; }

	/** Returns the objective function expression. @param var The variable matrix x[d][a]. @param mod The CPLEX model. **/
    IloExpr getObjFunction(IloNumVarMatrix &var, IloModel &mod);

	/** Returns the source constraint associated with a demand and a node. @param var The variable matrix x[d][a]. @param mod The CPLEX model. @param demand The demand. @param d The demand index. @param label The node label. **/
    IloRange getSourceConstraint_d_n(IloNumVarMatrix &var, IloModel &mod, const Demand & demand, int d, int label);
	
	/** Returns the flow conservation constraint associated with a demand and a node. @param var The variable matrix x[d][a]. @param mod The CPLEX model. @param v The node. @param demand The demand. @param d The demand index. **/
    IloRange getFlowConservationConstraint_i_d(IloNumVarMatrix &var, IloModel &mod, ListDigraph::Node &v, const Demand & demand, int d);

	/** Returns the target constraint associated with a demand. @param var The variable matrix x[d][a]. @param mod The CPLEX model. @param demand The demand. @param d The demand index. **/
    IloRange getTargetConstraint_d(IloNumVarMatrix &var, IloModel &mod, const Demand & demand, int d);
	
    /** Returns the length constraint associated with a demand. @param var The variable matrix x[d][a]. @param mod The CPLEX model. @param demand The demand. @param d The demand index. **/
    IloRange getLengthConstraint(IloNumVarMatrix &var, IloModel &mod, const Demand &demand, int d);

    /** Returns the non-overlapping constraint associated with an arc and a pair of demands. @param var The variable matrix x[d][a]. @param mod The CPLEX model. @param linkLabel The arc label. @param slice The arc slice. @param demand1 The first demand. @param d1 The first demand index. @param demand2 The second demand. @param d2 The second demand index. **/
    IloRange getNonOverlappingConstraint(IloNumVarMatrix &var, IloModel &mod, int linkLabel, int slice, const Demand & demand1, int d1, const Demand & demand2, int d2);
    
	/****************************************************************************************/
	/*										 Setters										*/
	/****************************************************************************************/

	/** Changes the index of an arc in a graph. @param a The arc. @param d The graph #d. @param val The new index value. **/
	void setArcIndex(const ListDigraph::Arc &a, int d, int val) { (*vecArcIndex[d])[a] = val; }

	/** Defines variables x[d][a] for every arc a in the extedend graph #d. @param var The variable matrix. @param mod The CPLEX model. **/
    void setVariables(IloNumVarMatrix &var, IloModel &mod);

	/** Defines the objective function. @param var The variable matrix. @param mod The CPLEX model. **/
    void setObjective(IloNumVarMatrix &var, IloModel &mod);

	/** Defines Source constraints. At most one arc leaves each node and exactly one arc leaves the source. @param var The variable matrix. @param mod The CPLEX model. **/
    void setSourceConstraints(IloNumVarMatrix &var, IloModel &mod);

	/** Defines Flow Conservation constraints. If an arc enters a node, then an arc must leave it. @param var The variable matrix. @param mod The CPLEX model. **/
    void setFlowConservationConstraints(IloNumVarMatrix &var, IloModel &mod);
	
	/** Defines Target constraints. Exactly one arc enters the demand's target. @param var The variable matrix. @param mod The CPLEX model. **/
    void setTargetConstraints(IloNumVarMatrix &var, IloModel &mod);
	
	/** Defines Length constraints. Demands must be routed within a length limit. @param var The variable matrix. @param mod The CPLEX model. **/
    void setLengthConstraints(IloNumVarMatrix &var, IloModel &mod);

	/** Defines Non-Overlapping constraints. Demands must not overlap eachother's slices. @param var The variable matrix. @param mod The CPLEX model. **/
    void setNonOverlappingConstraints(IloNumVarMatrix &var, IloModel &mod);
};    
#endif