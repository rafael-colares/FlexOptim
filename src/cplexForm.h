#ifndef __cplexForm__h
#define __cplexForm__h

#include "solver.h"

typedef IloArray<IloBoolVarArray> IloBoolVarMatrix;


/*********************************************************************************************
* This class implements and solve the Online Routing and Spectrum Allocation MIP using CPLEX.	
*********************************************************************************************/
class CplexForm : public Solver{

private:
    IloEnv env;						/**< The CPLEX environment. **/
    IloModel model;					/**< The CPLEX model. **/
    IloCplex cplex;					/**< The CPLEX engine. **/
    IloBoolVarMatrix x;				/**< The matrix of assignement variables used in the MIP. x[d][a]=1 if the d-th demand is routed through the arc from index a. **/
    IloIntVarArray maxSlicePerLink;	/**< The array of variables used in the MIP for verifying the max used slice position for each link in the topology network. maxSlicePerLink[i]=p if p is the max used slice position from the link with id i. **/
	IloIntVar maxSliceOverall;		/**< The max used slice position throughout all the network. **/
	static int count;				/**< Counts how many times CPLEX is called. **/

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

	/****************************************************************************************/
	/*										Destructor										*/
	/****************************************************************************************/

	/** Destructor. Clears the variable matrices, cplex model and environment. **/
	~CplexForm();

};


#endif