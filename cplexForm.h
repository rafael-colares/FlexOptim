#ifndef __Formulation__h
#define __Formulation__h

#include "ShortestPath.h"



class CplexForm : public ShortestPath{

private:
    IloEnv env;
    IloModel model;
    IloCplex cplex;
    IloNumVarArray x;

public:
    CplexForm(Instance &instance, const Demand &demand);
    void defineVariables();
    IloExpr getObjFunction();
    IloRange getShortestPathConstraint_i(ListDigraph::Node v);
    IloRange getLengthConstraint(const Demand &demand);
    void updatePath();
    void displayPathOnEnv();
    void displayOnPath();

};



#endif