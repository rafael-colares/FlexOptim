#ifndef __cplexForm__h
#define __cplexForm__h

#include "solver.h"

typedef IloArray<IloNumVarArray> IloNumVarMatrix;


/*********************************************************************************************
* This class implements and solve the Online Routing and Spectrum Allocation MIP using CPLEX.	
*********************************************************************************************/
class CplexForm : public Solver{

private:
    IloEnv env;			/**< The CPLEX environment. **/
    IloModel model;		/**< The CPLEX model. **/
    IloCplex cplex;		/**< The CPLEX engine. **/
    IloNumVarMatrix x;	/**< The matrix of variables used in the MIP. x[i][j]=1 if the i-th demand to be routed is assigned to the arc from id j. **/
    static int count;	/**< Counts how many times CPLEX is called. **/

public:

	/****************************************************************************************/
	/*										Constructors									*/
	/****************************************************************************************/
	/** Constructor. Builds the Online RSA mixed-integer program and solves it using CPLEX.  @param instance The instance to be solved. **/
    CplexForm(const Instance &instance);

	/****************************************************************************************/
	/*										Getters											*/
	/****************************************************************************************/
	/** Returns the cplex engine in use. **/
    IloCplex getCplex(){ return cplex; }

	/** Returns the total number of CPLEX default cuts applied during optimization. **/
	IloInt getNbCutsFromCplex();

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

};


#endif