#include "solver.h"


Solver::Solver(const Instance &inst) : RSA(inst) {
    std::cout << "--- Solver has been initalized ---" << std::endl;
}

/* Define variables x[a][d] for every arc a in the extedend graph and every demand d to be routed. */
void Solver::setVariables(IloNumVarMatrix &var, IloModel &mod){
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a){
        int arc = arcId[a];
        int nbDemandsToBeRouted = getNbDemandsToBeRouted();
        var[arc] = IloNumVarArray(mod.getEnv(), nbDemandsToBeRouted);
        for (int d = 0; d < nbDemandsToBeRouted; d++){    
            std::ostringstream varName;
            varName << "x";
            varName << "(" + std::to_string(nodeLabel[g.source(a)]+1) + "," + std::to_string(nodeLabel[g.target(a)]+1);
            varName << "," + std::to_string(arcSlice[a]+1) + "," + std::to_string(getToBeRouted()[d].getId()+1);
            varName << ")";
            IloNum upperBound = 1.0;
            if (instance.isRoutable(arcLabel[a], arcSlice[a], getToBeRouted()[d]) == false){
                upperBound = 0.0;
            }
            var[arc][d] = IloNumVar(mod.getEnv(), 0.0, upperBound, ILOINT, varName.str().c_str());
            mod.add(var[arc][d]);
        }
        //std::cout << "Creating variable " << id << ": " << varName.str().c_str() << std::endl;
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
        for (ListDigraph::ArcIt a(g); a != INVALID; ++a){
            int arc = arcId[a];
            double coeff = getCoeff(a, getToBeRouted()[d]);
            obj += coeff*var[arc][d];
        }
    }
    return obj;
}


/* Source constraints. At most 1 leaves each node. Exactly 1 leaves the Source. */
void Solver::setSourceConstraints(IloNumVarMatrix &var, IloModel &mod){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){  
        for (ListDigraph::NodeIt v(g); v != INVALID; ++v){
            IloRange sourceConstraint = getSourceConstraint_d(var, mod, getToBeRouted()[d], d, nodeLabel[v]);
            mod.add(sourceConstraint);
        } 
    }
}

/* Get an specific Source constraint */
IloRange Solver::getSourceConstraint_d(IloNumVarMatrix &var, IloModel &mod, const Demand & demand, int d, int i){
    IloExpr exp(mod.getEnv());
    IloInt upperBound = 1;
    IloInt lowerBound = 0;
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v){
        if (nodeLabel[v] == i){
            for (ListDigraph::OutArcIt a(g, v); a != INVALID; ++a){
                int arc = arcId[a];
                exp += var[arc][d];
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
        for (ListDigraph::NodeIt v(g); v != INVALID; ++v){
            if( (nodeLabel[v] != getToBeRouted()[d].getSource()) && (nodeLabel[v] != getToBeRouted()[d].getTarget()) ){
                IloRange st = getFlowConservationConstraint_i_d(var, mod, v, getToBeRouted()[d], d);
                mod.add(st);
            }
        }
    }
}

/* Get an specific Flow Conservation constraint */
IloRange Solver::getFlowConservationConstraint_i_d(IloNumVarMatrix &var, IloModel &mod, ListDigraph::Node &v, const Demand & demand, int d){
    IloExpr exp(mod.getEnv());
    IloInt rhs = 0;
    for (ListDigraph::OutArcIt a(g, v); a != INVALID; ++a){
        int arc = arcId[a];
        exp += var[arc][d];
    }
    for (ListDigraph::InArcIt a(g, v); a != INVALID; ++a){
        int arc = arcId[a];
        exp += (-1)*var[arc][d];
    }
    std::ostringstream constraintName;
    constraintName << "Flow(" << nodeLabel[v]+1 << "," << nodeSlice[v]+1 << "," << demand.getId()+1 << ")";
    IloRange constraint(mod.getEnv(), rhs, exp, rhs, constraintName.str().c_str());
    exp.end();
    return constraint;
}

/* Target constraints. Only 1 enters the Target */
void Solver::setTargetConstraints(IloNumVarMatrix &var, IloModel &mod){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){   
        IloRange targetConstraint = getTargetConstraint_d(var, mod, getToBeRouted()[d], d);
        mod.add(targetConstraint);
    }
}

/* Get an specific Target constraint */
IloRange Solver::getTargetConstraint_d(IloNumVarMatrix &var, IloModel &mod, const Demand & demand, int d){
    IloExpr exp(mod.getEnv());
    IloInt rhs = 1;
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v){
        int node = nodeLabel[v];
        if (node == demand.getTarget()){
            for (ListDigraph::InArcIt a(g, v); a != INVALID; ++a){
                int arc = arcId[a];
                exp += var[arc][d];
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
        IloRange lengthConstraint = getLengthConstraint(var, mod, getToBeRouted()[d], d);
        mod.add(lengthConstraint);
    }
}

/* Get an specific Length constraint */
IloRange Solver::getLengthConstraint(IloNumVarMatrix &var, IloModel &mod, const Demand &demand, int d){
    IloExpr exp(mod.getEnv());
    double rhs = demand.getMaxLength();
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a){
        int arc = arcId[a];
        double coeff = length[a];
        exp += coeff*var[arc][d];
    }
    std::ostringstream constraintName;
    constraintName << "Length(" << demand.getId()+1 << ")";
    IloRange constraint(mod.getEnv(), -IloInfinity, exp, rhs, constraintName.str().c_str());
    exp.end();
    return constraint;
}

/* Non-Overlapping constraints. Demands must not overlap eachother's slices */
void Solver::setNonOverlappingConstraints(IloNumVarMatrix &var, IloModel &mod){
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a){
        for (int d1 = 0; d1 < getNbDemandsToBeRouted(); d1++){
            for (int d2 = 0; d2 < getNbDemandsToBeRouted(); d2++){
                if(d1 != d2){
                    IloRange nonOverlap = getNonOverlappingConstraint(var, mod, arcLabel[a], arcSlice[a], getToBeRouted()[d1], d1, getToBeRouted()[d2], d2);
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
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a){
        if( (arcLabel[a] == linkLabel) && (arcSlice[a] == slice) ){
            int id = arcId[a];
            exp += var[id][d1];
        }
    }
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a){
        if( (arcLabel[a] == linkLabel) && (arcSlice[a] >= slice - demand1.getLoad() + 1) && (arcSlice[a] <= slice + demand2.getLoad() - 1) ){
            int id = arcId[a];
            exp += var[id][d2];
        }
    }
    std::ostringstream constraintName;
    constraintName << "Subcycle(" << linkLabel+1 << "," << slice+1 << "," << demand1.getId()+1 << "," << demand2.getId()+1 << ")";
    IloRange constraint(mod.getEnv(), -IloInfinity, exp, rhs, constraintName.str().c_str());
    exp.end();
    return constraint;
}


