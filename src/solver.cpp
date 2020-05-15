#include "solver.h"


/****************************************************************************************/
/*										Constructor										*/
/****************************************************************************************/

/* Constructor. The RSA constructor is called and the arc map storing the index of the preprocessed graphs associated is built. */
Solver::Solver(const Instance &inst) : RSA(inst) {
    std::cout << "--- Solver has been initalized ---" << std::endl;
}

/****************************************************************************************/
/*										Variables    									*/
/****************************************************************************************/

/* Define variables x[d][a] for every arc a in the extedend graph #d. */
void Solver::setVariables(IloBoolVarMatrix &var, IloIntVarArray &maxSliceFromLink, IloIntVar &maxSlice, IloModel &mod){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){ 
        var[d] = IloBoolVarArray(mod.getEnv(), countArcs(*vecGraph[d]));  
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            int arc = getArcIndex(a, d); 
            int label = getArcLabel(a, d); 
            int labelSource = getNodeLabel((*vecGraph[d]).source(a), d);
            int labelTarget = getNodeLabel((*vecGraph[d]).target(a), d);
            int slice = getArcSlice(a, d);
            std::ostringstream varName;
            varName << "x";
            varName << "(" + std::to_string(getToBeRouted_k(d).getId() + 1) + "," ;
            varName <<  std::to_string(labelSource + 1) + "," + std::to_string(labelTarget + 1) + ",";
            varName <<  std::to_string(slice + 1) + ")";
            IloInt upperBound = 1;
            if (instance.hasEnoughSpace(label, slice, getToBeRouted_k(d)) == false){
                upperBound = 0;
                std::cout << "STILL REMOVING VARIABLES IN CPLEX. \n" ;
            }
            var[d][arc] = IloBoolVar(mod.getEnv(), 0, upperBound, varName.str().c_str());
            mod.add(var[d][arc]);
            // std::cout << "Created variable: " << var[d][arc].getName() << std::endl;
        }
    }

    if(instance.getInput().getChosenObj() == Input::OBJECTIVE_METRIC_1p){
        for (int i = 0; i < instance.getNbEdges(); i++){
            std::string varName = "maxSlice(" + std::to_string(instance.getPhysicalLinkFromIndex(i).getId() + 1) + ")";
            IloInt lowerBound = instance.getPhysicalLinkFromIndex(i).getMaxUsedSlicePosition();
            IloInt upperBound = instance.getPhysicalLinkFromIndex(i).getNbSlices();
            maxSliceFromLink[i] = IloIntVar(mod.getEnv(), lowerBound, upperBound, varName.c_str());
            mod.add(maxSliceFromLink[i]);
        }
    }

    if(instance.getInput().getChosenObj() == Input::OBJECTIVE_METRIC_8){
        std::string varName = "maxSlice";
        IloInt lowerBound = instance.getMaxUsedSlicePosition();
        IloInt upperBound = instance.getMaxSlice();
        maxSlice = IloIntVar(mod.getEnv(), lowerBound, upperBound, varName.c_str()); 
        mod.add(maxSlice);
    }
}

/****************************************************************************************/
/*									Objective Function    								*/
/****************************************************************************************/

/* Set the objective Function */
void Solver::setObjective(IloBoolVarMatrix &var, IloIntVarArray &maxSliceFromLink, IloIntVar &maxSlice, IloModel &mod){
    IloExpr objective = getObjFunction(var, maxSliceFromLink, maxSlice, mod);
    mod.add(IloMinimize(mod.getEnv(), objective));
    objective.end();
}

/* Returns the objective function expression. */
IloExpr Solver::getObjFunction(IloBoolVarMatrix &var, IloIntVarArray &maxSliceFromLink, IloIntVar &maxSlice, IloModel &mod){
    IloExpr obj(mod.getEnv());
    
    if(instance.getInput().getChosenObj() == Input::OBJECTIVE_METRIC_1p){
        for (int i = 0; i < instance.getNbEdges(); i++){
            obj += maxSliceFromLink[i];
        }
        return obj;
    }

    if(instance.getInput().getChosenObj() == Input::OBJECTIVE_METRIC_8){
        obj += maxSlice;
        return obj;
    }


    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            int arc = getArcIndex(a, d); 
            double coeff = getCoeff(a, d);
            //coeff += (instance.getInput().getInitialLagrangianMultiplier() * getArcLength(a, 0) );
            obj += coeff*var[d][arc];
        }
    }
    return obj;
}

/****************************************************************************************/
/*										Constraints    									*/
/****************************************************************************************/

/* Defines Source constraints. At most one arc leaves each node and exactly one arc leaves the source. */
void Solver::setSourceConstraints(IloBoolVarMatrix &var, IloModel &mod){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){  
        for (ListDigraph::NodeIt v(*vecGraph[d]); v != INVALID; ++v){
            int label = getNodeLabel(v, d);
            IloRange sourceConstraint = getSourceConstraint_d_n(var, mod, getToBeRouted_k(d), d, label);
            mod.add(sourceConstraint);
        } 
    }
}

/* Returns the source constraint associated with a demand and a node. */
IloRange Solver::getSourceConstraint_d_n(IloBoolVarMatrix &var, IloModel &mod, const Demand & demand, int d, int i){
    IloExpr exp(mod.getEnv());
    IloInt upperBound = 1;
    IloInt lowerBound = 0;
    for (ListDigraph::NodeIt v(*vecGraph[d]); v != INVALID; ++v){
        if (getNodeLabel(v, d) == i){
            for (ListDigraph::OutArcIt a((*vecGraph[d]), v); a != INVALID; ++a){
                int arc = getArcIndex(a, d); 
                exp += var[d][arc];
            }
        }
    }
    std::ostringstream constraintName;
    constraintName << "Source(" << i+1 << "," << demand.getId()+1 << ")";
    if (i == demand.getSource()){
        lowerBound = 1;
    }
    if (i == demand.getTarget()){
        upperBound = 0;
    }
    IloRange constraint(mod.getEnv(), lowerBound, exp, upperBound, constraintName.str().c_str());
    exp.end();
    return constraint;
}

/* Defines Flow Conservation constraints. If an arc enters a node, then an arc must leave it. */
void Solver::setFlowConservationConstraints(IloBoolVarMatrix &var, IloModel &mod){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){   
        for (ListDigraph::NodeIt v(*vecGraph[d]); v != INVALID; ++v){
            int label = getNodeLabel(v, d);
            if( (label != getToBeRouted_k(d).getSource()) && (label != getToBeRouted_k(d).getTarget()) ){
                IloRange st = getFlowConservationConstraint_i_d(var, mod, v, getToBeRouted_k(d), d);
                mod.add(st);
            }
        }
    }
}

/* Returns the flow conservation constraint associated with a demand and a node. */
IloRange Solver::getFlowConservationConstraint_i_d(IloBoolVarMatrix &var, IloModel &mod, ListDigraph::Node &v, const Demand & demand, int d){
    IloExpr exp(mod.getEnv());
    IloInt rhs = 0;
    for (ListDigraph::OutArcIt a((*vecGraph[d]), v); a != INVALID; ++a){
        int arc = getArcIndex(a, d); 
        exp += var[d][arc];
    }
    for (ListDigraph::InArcIt a((*vecGraph[d]), v); a != INVALID; ++a){
        int arc = getArcIndex(a, d); 
        exp += (-1)*var[d][arc];
    }
    std::ostringstream constraintName;
    int label = getNodeLabel(v, d);
    int slice = getNodeSlice(v, d);
    constraintName << "Flow(" << label+1 << "," << slice+1 << "," << demand.getId()+1 << ")";
    IloRange constraint(mod.getEnv(), rhs, exp, rhs, constraintName.str().c_str());
    exp.end();
    return constraint;
}

/* Defines Target constraints. Exactly one arc enters the target. */
void Solver::setTargetConstraints(IloBoolVarMatrix &var, IloModel &mod){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){   
        IloRange targetConstraint = getTargetConstraint_d(var, mod, getToBeRouted_k(d), d);
        mod.add(targetConstraint);
    }
}

/* Returns the target constraint associated with a demand. */
IloRange Solver::getTargetConstraint_d(IloBoolVarMatrix &var, IloModel &mod, const Demand & demand, int d){
    IloExpr exp(mod.getEnv());
    IloInt rhs = 1;
    for (ListDigraph::NodeIt v(*vecGraph[d]); v != INVALID; ++v){
        int label = getNodeLabel(v, d);
        if (label == demand.getTarget()){
            for (ListDigraph::InArcIt a((*vecGraph[d]), v); a != INVALID; ++a){
                int arc = getArcIndex(a, d); 
                exp += var[d][arc];
            }
        }
    }
    std::ostringstream constraintName;
    constraintName << "Target(" << demand.getId()+1 << ")";
    IloRange constraint(mod.getEnv(), rhs, exp, rhs, constraintName.str().c_str());
    exp.end();
    return constraint;
}

/* Defines Length constraints. Demands must be routed within a length limit. */
void Solver::setLengthConstraints(IloBoolVarMatrix &var, IloModel &mod){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){   
        IloRange lengthConstraint = getLengthConstraint(var, mod, getToBeRouted_k(d), d);
        mod.add(lengthConstraint);
    }
}

/* Returns the length constraint associated with a demand. */
IloRange Solver::getLengthConstraint(IloBoolVarMatrix &var, IloModel &mod, const Demand &demand, int d){
    IloExpr exp(mod.getEnv());
    double rhs = demand.getMaxLength();
    for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
        int arc = getArcIndex(a, d); 
        double coeff = getArcLength(a, d);
        exp += coeff*var[d][arc];
    }
    std::ostringstream constraintName;
    constraintName << "Length(" << demand.getId()+1 << ")";
    IloRange constraint(mod.getEnv(), -IloInfinity, exp, rhs, constraintName.str().c_str());
    exp.end();
    return constraint;
}

/* Defines Non-Overlapping constraints. Demands must not overlap eachother's slices. */
void Solver::setNonOverlappingConstraints(IloBoolVarMatrix &var, IloModel &mod){
    for (int d1 = 0; d1 < getNbDemandsToBeRouted(); d1++){
        for (int d2 = 0; d2 < getNbDemandsToBeRouted(); d2++){
            if(d1 != d2){
                for (int i = 0; i < instance.getNbEdges(); i++){
                    for (int s = 0; s < instance.getPhysicalLinkFromIndex(i).getNbSlices(); s++){
                        IloRange nonOverlap = getNonOverlappingConstraint(var, mod, instance.getPhysicalLinkFromIndex(i).getId(), s, getToBeRouted_k(d1), d1, getToBeRouted_k(d2), d2);
                        mod.add(nonOverlap);
                    }
                }   
            }
        }
    }
}

/* Returns the non-overlapping constraint associated with an arc and a pair of demands. */
IloRange Solver::getNonOverlappingConstraint(IloBoolVarMatrix &var, IloModel &mod, int linkLabel, int slice, const Demand & demand1, int d1, const Demand & demand2, int d2){
    IloExpr exp(mod.getEnv());
    IloNum rhs = 1;
    for (ListDigraph::ArcIt a(*vecGraph[d1]); a != INVALID; ++a){
        if( (getArcLabel(a, d1) == linkLabel) && (getArcSlice(a, d1) == slice) ){
            int id = getArcIndex(a, d1);
            exp += var[d1][id];
        }
    }
    for (ListDigraph::ArcIt a(*vecGraph[d2]); a != INVALID; ++a){
        if( (getArcLabel(a, d2) == linkLabel) && (getArcSlice(a, d2) >= slice - demand1.getLoad() + 1) && (getArcSlice(a, d2) <= slice + demand2.getLoad() - 1) ){
            int id = getArcIndex(a, d2);
            exp += var[d2][id];
        }
    }
    std::ostringstream constraintName;
    constraintName << "Subcycle(" << linkLabel+1 << "," << slice+1 << "," << demand1.getId()+1 << "," << demand2.getId()+1 << ")";
    IloRange constraint(mod.getEnv(), -IloInfinity, exp, rhs, constraintName.str().c_str());
    exp.end();
    return constraint;
}

/* Defines the Link's Max Used Slice Position constraints. The max used slice position on each link must be greater than every slice position used in the link. */
void Solver::setMaxUsedSlicePerLinkConstraints(IloBoolVarMatrix &var, IloIntVarArray &maxSlicePerLink, IloModel &mod){
    for (int i = 0; i < instance.getNbEdges(); i++){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            IloRange maxUsedSlicePerLinkConst = getMaxUsedSlicePerLinkConstraints(var, maxSlicePerLink, i, d, mod);
            mod.add(maxUsedSlicePerLinkConst);
        }
    }
}

IloRange Solver::getMaxUsedSlicePerLinkConstraints(IloBoolVarMatrix &var, IloIntVarArray &maxSlicePerLink, int linkIndex, int d, IloModel &mod){
    IloExpr exp(mod.getEnv());
    IloInt rhs = 0;
    int linkLabel = instance.getPhysicalLinkFromIndex(linkIndex).getId();
    for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
        if (getArcLabel(a, d) == linkLabel){
            int index = getArcIndex(a, d);
            int slice = getArcSlice(a, d);
            exp += slice*var[d][index];
        }
    }
    exp += -maxSlicePerLink[linkIndex];
    
    std::ostringstream constraintName;
    constraintName << "MaxUsedSlicePerLink(" << linkLabel+1 << "," << getToBeRouted_k(d).getId()+1 << ")";
    IloRange constraint(mod.getEnv(), -IloInfinity, exp, rhs, constraintName.str().c_str());
    exp.end();
    return constraint;
}

/* Defines the Overall Max Used Slice Position constraints. */
void Solver::setMaxUsedSliceOverallConstraints(IloBoolVarMatrix &var, IloIntVar maxSliceOverall, IloModel &mod){
    for (int i = 0; i < instance.getNbEdges(); i++){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            IloRange maxUsedSliceOverallConst = getMaxUsedSliceOverallConstraints(var, maxSliceOverall, i, d, mod);
            mod.add(maxUsedSliceOverallConst);
        }
    }
}


IloRange Solver::getMaxUsedSliceOverallConstraints(IloBoolVarMatrix &var, IloIntVar &maxSlice, int linkIndex, int d, IloModel &mod){
    IloExpr exp(mod.getEnv());
    IloInt rhs = 0;
    int linkLabel = instance.getPhysicalLinkFromIndex(linkIndex).getId();
    for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
        if (getArcLabel(a, d) == linkLabel){
            int index = getArcIndex(a, d);
            int slice = getArcSlice(a, d);
            exp += slice*var[d][index];
        }
    }
    exp += -maxSlice;
    
    std::ostringstream constraintName;
    constraintName << "MaxUsedSliceOverall(" << linkLabel+1 << "," << getToBeRouted_k(d).getId()+1 << ")";
    IloRange constraint(mod.getEnv(), -IloInfinity, exp, rhs, constraintName.str().c_str());
    exp.end();
    return constraint;
}