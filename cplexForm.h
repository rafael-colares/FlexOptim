#ifndef __cplexForm__h
#define __cplexForm__h

#include "RSA.h"


typedef IloArray<IloNumVarArray> IloNumVarMatrix;


class CplexForm : public RSA{

private:
    IloEnv env;
    IloModel model;
    IloCplex cplex;
    IloNumVarMatrix x;
    std::vector<Demand> toBeRouted; 

public:
    static int count;
	/************************************************/
	/*					Constructors				*/
	/************************************************/
    CplexForm(const Instance &inst);
    void defineVariables();

	/************************************************/
	/*					   Getters 		    		*/
	/************************************************/
    IloExpr getObjFunction();
    IloRange getFlowConservationConstraint_i_d(ListDigraph::Node v, const Demand & demand, int d);
    IloRange getSourceConstraint_d(const Demand & demand, int d, int i);
    IloRange getTargetConstraint_d(const Demand & demand, int d);
    IloRange getLengthConstraint(const Demand &demand, int d);
    IloRange getNonOverlappingConstraint(int linkLabel, int slice, const Demand & demand1, int d1, const Demand & demand2, int d2);
    std::vector<Demand> getToBeRouted() { return toBeRouted; } 
    int getNbDemandsToBeRouted() { return toBeRouted.size();}
    IloCplex getCplex(){ return cplex; }

	/************************************************/
	/*					   Setters 		    		*/
	/************************************************/
    void setToBeRouted(const std::vector<Demand> &d){this->toBeRouted = d;}
    void setObjective();
    static void setCount(int i){ count = i; }
    void addSourceConstraints();
    void addTargetConstraints();
    void addFlowConservationConstraints();
    void addLengthConstraints();
    void addNonOverlappingConstraints();

    void updatePath();
    //void displayPathOnEnv();
    void displayOnPath();
    void displayToBeRouted();
    void displayVariableValues();

};


#endif