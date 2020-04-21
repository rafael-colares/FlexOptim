#include "solver.h"


Solver::Solver(const Instance &inst) : RSA(inst) {
    std::cout << "--- Solver has been initalized ---" << std::endl;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){ 
        vecArcIndex.emplace_back(new ArcMap((*vecGraph[d]), -1));
        int index=0;
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            setArcIndex(a, d, index);
            index++;
        }
    }
}

/* Define variables x[a][d] for every arc a in the extedend graph and every demand d to be routed. */
void Solver::setVariables(IloNumVarMatrix &var, IloModel &mod){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){ 
        var[d] = IloNumVarArray(mod.getEnv(), countArcs(*vecGraph[d]));  
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
            IloNum upperBound = 1.0;
            if (instance.hasEnoughSpace(label, slice, getToBeRouted_k(d)) == false){
                upperBound = 0.0;
                std::cout << "STILL REMOVING VARIABLES IN CPLEX. \n" ;
            }
            var[d][arc] = IloNumVar(mod.getEnv(), 0.0, upperBound, ILOINT, varName.str().c_str());
            mod.add(var[d][arc]);    
            // std::cout << "Created variable: " << var[d][arc].getName() << std::endl;
        }
    }
}

/* Set the objective Function */
void Solver::setObjective(IloNumVarMatrix &var, IloModel &mod){
    IloExpr objective = getObjFunction(var, mod);
    mod.add(IloMinimize(mod.getEnv(), objective));
    objective.end();
}

/* Get an Objective Function */
IloExpr Solver::getObjFunction(IloNumVarMatrix &var, IloModel &mod){
    IloExpr obj(mod.getEnv());
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

/* Source constraints. At most 1 leaves each node. Exactly 1 leaves the Source. */
void Solver::setSourceConstraints(IloNumVarMatrix &var, IloModel &mod){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){  
        for (ListDigraph::NodeIt v(*vecGraph[d]); v != INVALID; ++v){
            int label = getNodeLabel(v, d);
            IloRange sourceConstraint = getSourceConstraint_d(var, mod, getToBeRouted_k(d), d, label);
            mod.add(sourceConstraint);
        } 
    }
}

/* Get an specific Source constraint */
IloRange Solver::getSourceConstraint_d(IloNumVarMatrix &var, IloModel &mod, const Demand & demand, int d, int i){
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

/* Flow constraints. Everything that enters must go out. */
void Solver::setFlowConservationConstraints(IloNumVarMatrix &var, IloModel &mod){
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

/* Get an specific Flow Conservation constraint */
IloRange Solver::getFlowConservationConstraint_i_d(IloNumVarMatrix &var, IloModel &mod, ListDigraph::Node &v, const Demand & demand, int d){
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

/* Target constraints. Only 1 enters the Target */
void Solver::setTargetConstraints(IloNumVarMatrix &var, IloModel &mod){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){   
        IloRange targetConstraint = getTargetConstraint_d(var, mod, getToBeRouted_k(d), d);
        mod.add(targetConstraint);
    }
}

/* Get an specific Target constraint */
IloRange Solver::getTargetConstraint_d(IloNumVarMatrix &var, IloModel &mod, const Demand & demand, int d){
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

/* Length Constraints. Demands must be routed within a length limit */
void Solver::setLengthConstraints(IloNumVarMatrix &var, IloModel &mod){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){   
        IloRange lengthConstraint = getLengthConstraint(var, mod, getToBeRouted_k(d), d);
        mod.add(lengthConstraint);
    }
}

/* Get an specific Length constraint */
IloRange Solver::getLengthConstraint(IloNumVarMatrix &var, IloModel &mod, const Demand &demand, int d){
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

/* Non-Overlapping constraints. Demands must not overlap eachother's slices */
void Solver::setNonOverlappingConstraints(IloNumVarMatrix &var, IloModel &mod){
    for (int d1 = 0; d1 < getNbDemandsToBeRouted(); d1++){
        for (ListDigraph::ArcIt a(*vecGraph[d1]); a != INVALID; ++a){
            for (int d2 = 0; d2 < getNbDemandsToBeRouted(); d2++){
                if(d1 != d2){
                    IloRange nonOverlap = getNonOverlappingConstraint(var, mod, getArcLabel(a, d1), getArcSlice(a, d1), getToBeRouted_k(d1), d1, getToBeRouted_k(d2), d2);
                    mod.add(nonOverlap);
                }   
            }
        }
    }
}

/* Get an specific Non-Overlapping constraint */
IloRange Solver::getNonOverlappingConstraint(IloNumVarMatrix &var, IloModel &mod, int linkLabel, int slice, const Demand & demand1, int d1, const Demand & demand2, int d2){
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


