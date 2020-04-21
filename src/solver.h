#ifndef __solver__h
#define __solver__h

#include "RSA.h"


typedef IloArray<IloNumVarArray> IloNumVarMatrix;


class Solver : public RSA{
protected:
	std::vector< std::shared_ptr<ArcMap> > vecArcIndex;

public:
	/************************************************/
	/*				    Constructors 		   		*/
	/************************************************/
    Solver(const Instance &inst);

	/************************************************/
	/*					   Getters 		    		*/
	/************************************************/
    int getArcIndex(const ListDigraph::Arc &a, int d) const { return (*vecArcIndex[d])[a]; }

    IloExpr getObjFunction(IloNumVarMatrix &var, IloModel &mod);
    IloRange getSourceConstraint_d(IloNumVarMatrix &var, IloModel &mod, const Demand & demand, int d, int i);
    IloRange getFlowConservationConstraint_i_d(IloNumVarMatrix &var, IloModel &mod, ListDigraph::Node &v, const Demand & demand, int d);
    IloRange getTargetConstraint_d(IloNumVarMatrix &var, IloModel &mod, const Demand & demand, int d);
    IloRange getLengthConstraint(IloNumVarMatrix &var, IloModel &mod, const Demand &demand, int d);
    IloRange getNonOverlappingConstraint(IloNumVarMatrix &var, IloModel &mod, int linkLabel, int slice, const Demand & demand1, int d1, const Demand & demand2, int d2);
    
	/************************************************/
	/*					   Setters 		    		*/
	/************************************************/

	void setArcIndex(const ListDigraph::Arc &a, int d, int val) { (*vecArcIndex[d])[a] = val; }

	/* Define variables x[a][d] for every arc a in the extedend graph and every demand d to be routed. */
    void setVariables(IloNumVarMatrix &var, IloModel &mod);

    void setObjective(IloNumVarMatrix &var, IloModel &mod);

	/* Source constraints. At most 1 leaves each node. Exactly 1 leaves the Source. */
    void setSourceConstraints(IloNumVarMatrix &var, IloModel &mod);

	/* Flow constraints. Everything that enters must go out. */
    void setFlowConservationConstraints(IloNumVarMatrix &var, IloModel &mod);
	
	/* Target constraints. Only 1 enters the Target */
    void setTargetConstraints(IloNumVarMatrix &var, IloModel &mod);
	
	/* Length Constraints. Demands must be routed within a length limit */
    void setLengthConstraints(IloNumVarMatrix &var, IloModel &mod);

	/* Non-Overlapping constraints. Demands must not overlap eachother's slices */
    void setNonOverlappingConstraints(IloNumVarMatrix &var, IloModel &mod);
};    
#endif