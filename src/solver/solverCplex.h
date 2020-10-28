#ifndef __solverCplex__h
#define __solverCplex__h

#include <ilcplex/ilocplex.h>
#include "abstractSolver.h"

/************************************************************************************
 * This is the class implementing the generic callback interface. It has two main 
 * functions: addUserCuts and addLazyConstraints.												
 ************************************************************************************/
class CplexCallback: public IloCplex::Callback::Function {
private:
	IloNumVarArray var;
	AbstractFormulation *formulation;

public:
	// Constructor with data.
	CplexCallback(const IloNumVarArray _var, AbstractFormulation* &_formulation);

	void addUserCuts (const IloCplex::Callback::Context &context) const; 
    
    void addLazyConstraints(const IloCplex::Callback::Context &context) const;
    
	virtual void invoke (const IloCplex::Callback::Context &context);

	std::vector<double> getIntegerSolution(const IloCplex::Callback::Context &context) const;
	std::vector<double> getFractionalSolution(const IloCplex::Callback::Context &context) const;
	IloExpr to_IloExpr(const IloCplex::Callback::Context &context, const Expression &e) const;
    // Destructor
    virtual ~CplexCallback();
};

/*********************************************************************************************
* This class implements the Online Routing and Spectrum Allocation through a Flow MIP 
* formulation using CPLEX.	
*********************************************************************************************/
class SolverCplex : public AbstractSolver{

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
    SolverCplex(const Instance &instance);

	/****************************************************************************************/
	/*										Getters											*/
	/****************************************************************************************/
	IloCplex getCplex(){ return cplex; }

	

	AbstractSolver::Status getStatus() override;

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
	void setCplexParams(const Input &input);
	/****************************************************************************************/
	/*										Methods											*/
	/****************************************************************************************/

	/** Recovers the obtained MIP solution and builds a path for each demand on its associated graph from RSA. **/
    void updatePath();

	void implementFormulation() override;

	void updateRSA(Instance &instance) override;
	
	void exportFormulation(const Instance &instance);
	
	void solve() override;

	IloExpr to_IloExpr(const Expression &e);
	
	/** Returns the total number of CPLEX default cuts applied during optimization. **/
	IloInt getNbCutsFromCplex();
	
	/* Builds file results.csv containing information about the main obtained results. */
	void outputLogResults(std::string fileName) override;
	/****************************************************************************************/
	/*										Display											*/
	/****************************************************************************************/

	/** Displays the obtained paths. **/
    void displayOnPath();

	/** Displays the value of each variable in the obtained solution. **/
    void displaySolution();

	/****************************************************************************************/
	/*										Destructor										*/
	/****************************************************************************************/

	/** Destructor. Clears the variable matrices, cplex model and environment. **/
	~SolverCplex();

};

#endif