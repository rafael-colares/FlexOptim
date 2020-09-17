#ifndef __solverCplex__h
#define __solverCplex__h

#include <ilcplex/ilocplex.h>
#include "FlowForm.h"
#include "solver.h"


/*********************************************************************************************
* This class implements the Online Routing and Spectrum Allocation through a Flow MIP 
* formulation using CPLEX.	
*********************************************************************************************/
class SolverCplex : public Solver{

private:
    IloEnv env;						/**< The CPLEX environment. **/
    IloModel model;					/**< The CPLEX model. **/
    IloCplex cplex;					/**< The CPLEX engine. **/
	static int count;				/**< Counts how many times the solver is called. **/
    IloNumVarArray var;				/**< The array of variables used in the MIP. **/
	IloObjective obj;

public:
	/****************************************************************************************/
	/*										Constructors									*/
	/****************************************************************************************/

	/** Constructor. Builds the Online RSA mixed-integer program and solves it using CPLEX.  @param instance The instance to be solved. **/
    SolverCplex(Instance &instance);

	/****************************************************************************************/
	/*										Getters											*/
	/****************************************************************************************/
	IloCplex getCplex(){ return cplex; }

	Solver::Status getStatus();

	std::vector<double> getSolution();
	/****************************************************************************************/
	/*										Setters											*/
	/****************************************************************************************/

	/** Defines the decision variables needed in the MIP formulation. **/
    void setVariables(const std::vector<Variable> &myVars);

	/** Defines the constraints needed in the MIP formulation. **/
    void setConstraints(const std::vector<Constraint> &myConstraints);

	/** Defines the objective function. **/
    void setObjective(const ObjectiveFunction &myObjective);

	/** Defines the cplex optimization parameters. **/
	void setCplexParams(const Instance &instance);

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

	void implementFormulation(const std::vector<Variable> &vars, const std::vector<Constraint> &constraints, const ObjectiveFunction &obj);
	
	void exportFormulation(const Instance &instance);
	
	void solve(const std::vector<ObjectiveFunction> &myObjectives);

	IloExpr to_IloExpr(const Expression &e);
	
	/** Returns the total number of CPLEX default cuts applied during optimization. **/
	IloInt getNbCutsFromCplex();
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
	~SolverCplex();

};


#endif