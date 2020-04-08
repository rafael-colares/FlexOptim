#ifndef __cplexForm__h
#define __cplexForm__h

#include "solver.h"


typedef IloArray<IloNumVarArray> IloNumVarMatrix;


class CplexForm : public Solver{

private:
    IloEnv env;
    IloModel model;
    IloCplex cplex;
    IloNumVarMatrix x;

public:
    static int count;
	/************************************************/
	/*					Constructors				*/
	/************************************************/
    CplexForm(const Instance &inst);

	/************************************************/
	/*					   Getters 		    		*/
	/************************************************/
    IloCplex getCplex(){ return cplex; }

	/************************************************/
	/*					   Setters 		    		*/
	/************************************************/
    static void setCount(int i){ count = i; }

	/************************************************/
	/*					   Methods 		    		*/
	/************************************************/
    void updatePath();
	IloInt getNbCutsFromCplex();

	/************************************************/
	/*					   Display 		    		*/
	/************************************************/
    void displayOnPath();
    void displayVariableValues();

};


#endif