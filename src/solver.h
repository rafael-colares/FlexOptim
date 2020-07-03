#ifndef __solver__h
#define __solver__h

#include "RSA.h"


//typedef IloArray<IloNumVarArray> IloNumVarMatrix;
typedef IloArray<IloBoolVarArray> IloBoolVarMatrix;


/**********************************************************************************************
 * This class specifies which MIP solver to be used The methods needed for solving 
 * the Routing and Spectrum Allocation problem via CPLEX. 
 * \note It uses the LEMON library to build the arc maps and CPLEX concert library to build 
 * expressions and constraints. 
 * *******************************************************************************************/
class Solver : public RSA{
protected:
    IloEnv env;						/**< The CPLEX environment. **/
    IloModel model;					/**< The CPLEX model. **/
    IloCplex cplex;					/**< The CPLEX engine. **/
	static int count;				/**< Counts how many times the solver is called. **/

public:
	/****************************************************************************************/
	/*										Constructor										*/
	/****************************************************************************************/

	/** Constructor. The RSA constructor is called and the arc map storing the index of the preprocessed graphs associated is built. @param inst The instance to be solved. **/
    Solver(const Instance &inst);

	/****************************************************************************************/
	/*										 Getters										*/
	/****************************************************************************************/

	/** Returns the cplex engine in use. **/
    IloCplex getCplex(){ return cplex; }

	/** Returns the total number of CPLEX default cuts applied during optimization. **/
	IloInt getNbCutsFromCplex();

	/** Returns the cplex status. **/
    Status getStatus();

	/****************************************************************************************/
	/*									Destructor											*/
	/****************************************************************************************/
	/** Destructor. Free solver memory. **/
	~Solver();
};    
#endif