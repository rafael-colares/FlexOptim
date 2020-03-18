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
	/************************************************/
	/*					Constructors				*/
	/************************************************/
    CplexForm(const Instance &inst);
    void defineVariables();

	/************************************************/
	/*					   Getters 		    		*/
	/************************************************/
    IloExpr getObjFunction();
    IloRange getShortestPathConstraint_i_d(ListDigraph::Node v, const Demand & demand, int d);
    IloRange getSourceConstraint_d(const Demand & demand, int d, int i);
    IloRange getTargetConstraint_d(const Demand & demand, int d);
    IloRange getLengthConstraint(const Demand &demand, int d);
    IloRange getSubcycleConstraint(const ListDigraph::Arc &a, const Demand & demand, int d);
    std::vector<Demand> getToBeRouted() { return toBeRouted; } 
    int getNbDemandsToBeRouted() { return toBeRouted.size();}
    IloCplex getCplex(){ return cplex; }

	/************************************************/
	/*					   Setters 		    		*/
	/************************************************/
    void setToBeRouted(const std::vector<Demand> &d){this->toBeRouted = d;}

    void updatePath();
    //void displayPathOnEnv();
    void displayOnPath();
    void displayToBeRouted();

};



#endif