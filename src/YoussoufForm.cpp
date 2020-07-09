#include "YoussoufForm.h"

#include <sstream>


/* Constructor. Builds the Online RSA mixed-integer program and solves it using CPLEX. */
YoussoufForm::YoussoufForm(const Instance &inst) : 
                            Solver(inst), x(env, countEdges(compactGraph)), 
                            z(env, instance.getPhysicalLinkFromIndex(0).getNbSlices()), 
                            t(env, countEdges(compactGraph)), 
                            maxSlicePerLink(env, instance.getNbEdges()), maxSliceOverall(env){
    std::cout << "--- Youssouf formulation has been chosen ---" << std::endl;
    /************************************************/
	/*				    SET VARIABLES				*/
	/************************************************/
    this->setVariables();
    std::cout << "Variables have been defined..." << std::endl;

	/************************************************/
	/*			    SET OBJECTIVE FUNCTION			*/
	/************************************************/
    this->setObjective();
    std::cout << "Objective function has been defined..." << std::endl;

	/************************************************/
	/*			      SET CONSTRAINTS				*/
	/************************************************/
    this->setOriginConstraints();
    std::cout << "Origin constraints have been defined..." << std::endl;
    this->setDestinationConstraints();
    std::cout << "Destination constraints have been defined..." << std::endl;
    this->setDegreeConstraints();
    std::cout << "Degree constraints have been defined..." << std::endl;
    this->setTransmissionReachConstraints();
    std::cout << "Transmission-reach constraints have been defined..." << std::endl;
    this->setChannelSelectionConstraints();
    std::cout << "Channel-selection constraints have been defined..." << std::endl;
    this->setForbiddenSlotConstraints();
    std::cout << "Forbidden-slot constraints have been defined..." << std::endl;
    this->setEdgeSlotConstraints();
    std::cout << "Edge-slot constraints have been defined..." << std::endl;
    this->setDemandEdgeSlotConstraints();
    std::cout << "Demand-edge-slot constraints have been defined..." << std::endl;
    this->setNonOverlappingConstraints();
    std::cout << "Non-overlapping constraints have been defined..." << std::endl;
    
    if(instance.getInput().getChosenObj() == Input::OBJECTIVE_METRIC_1p){
        //this->setMaxUsedSlicePerLinkConstraints(x, maxSlicePerLink, model);    
        std::cout << "Max Used Slice Per Link constraints have been defined..." << std::endl;
    }

    if(instance.getInput().getChosenObj() == Input::OBJECTIVE_METRIC_8){
        //this->setMaxUsedSliceOverallConstraints(x, maxSliceOverall, model);    
        std::cout << "Max Used Slice Overall constraints have been defined..." << std::endl;
    }
    
    
	/************************************************/
	/*		    EXPORT LINEAR PROGRAM TO .LP		*/
	/************************************************/
    std::string file = getInstance().getInput().getOutputPath() + "LP/model" + std::to_string(count) + ".lp";
    //cplex.exportModel(file.c_str());
    std::cout << "LP model has been exported..." << std::endl;
    
	/************************************************/
	/*             DEFINE CPLEX PARAMETERS   		*/
	/************************************************/
    cplex.setParam(IloCplex::Param::MIP::Display, 3);
    cplex.setParam(IloCplex::Param::TimeLimit, getInstance().getInput().getIterationTimeLimit());
    std::cout << "CPLEX parameters have been defined..." << std::endl;


	/************************************************/
	/*           SET UP THE GENERIC CALLBACK   		*/
	/************************************************/
    GenericCallback myGenericCallback(x, instance, compactGraph, toBeRouted, compactNodeLabel, compactNodeId, compactEdgeLabel, compactEdgeId);
    // The generic callback will be used in relaxation and candidate contexts.
    CPXLONG contextMask = 0;
    contextMask |= IloCplex::Callback::Context::Id::Relaxation;
    contextMask |= IloCplex::Callback::Context::Id::Candidate;
    cplex.use(&myGenericCallback, contextMask);

	/************************************************/
	/*		         SOLVE LINEAR PROGRAM   		*/
	/************************************************/
    IloNum timeStart = cplex.getCplexTime();
    std::cout << "Solving..." << std::endl;
    cplex.solve();
    std::cout << "Solved!" << std::endl;
    IloNum timeFinish = cplex.getCplexTime();

	/************************************************/
	/*		    GET OPTIMAL SOLUTION FOUND        	*/
	/************************************************/
    if ((cplex.getStatus() == IloAlgorithm::Optimal) || (cplex.getStatus() == IloAlgorithm::Feasible)){
        std::cout << "Optimization done in " << timeFinish - timeStart << " secs." << std::endl;
        std::cout << "Objective Function Value: " << cplex.getObjValue() << std::endl;
        //displayVariableValues();
        updatePath();
        //displayOnPath();
        
        std::cout << "Number of cplex cuts: " << getNbCutsFromCplex() << std::endl;
    }
    else{
        std::cout << "Could not find a feasible solution..." << std::endl;
    }
}




/****************************************************************************************/
/*										Variables    									*/
/****************************************************************************************/

/* Define variables needed in the MIP. */
void YoussoufForm::setVariables(){
    // Variables x[e][k]
    for (ListGraph::EdgeIt e(compactGraph); e != INVALID; ++e){
        int edge = getCompactEdgeLabel(e);
        x[edge] = IloBoolVarArray(model.getEnv(), getNbDemandsToBeRouted());  
        for (int k = 0; k < getNbDemandsToBeRouted(); k++){
            std::ostringstream varName;
            varName << "x";
            varName << "(" + std::to_string(edge) + "," ;
            varName <<  std::to_string(getToBeRouted_k(k).getId() + 1) + ")";
            IloInt upperBound = 1;
            x[edge][k] = IloBoolVar(model.getEnv(), 0, upperBound, varName.str().c_str());
            model.add(x[edge][k]);
        }
    }
    std::cout << "X variables were created." << std::endl;
    // Variables z[s][k]
    for (int s = 0; s < instance.getPhysicalLinkFromIndex(0).getNbSlices(); s++){
        z[s] = IloBoolVarArray(model.getEnv(), getNbDemandsToBeRouted());  
        for (int k = 0; k < getNbDemandsToBeRouted(); k++){
            std::ostringstream varName;
            varName << "z";
            varName << "(" + std::to_string(s) + "," ;
            varName <<  std::to_string(getToBeRouted_k(k).getId() + 1) + ")";
            IloInt upperBound = 1;
            z[s][k] = IloBoolVar(model.getEnv(), 0, upperBound, varName.str().c_str());
            model.add(z[s][k]);
            // std::cout << "Created variable: " << var[d][arc].getName() << std::endl;
        }
    }
    std::cout << "Z variables were created." << std::endl;
    // Variables t[e][s][k]
    for (ListGraph::EdgeIt e(compactGraph); e != INVALID; ++e){
        int edge = getCompactEdgeLabel(e);
        t[edge] = IloBoolVarMatrix(model.getEnv(), instance.getPhysicalLinkFromIndex(edge).getNbSlices());  
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(edge).getNbSlices(); s++){
            t[edge][s] = IloBoolVarArray(model.getEnv(), getNbDemandsToBeRouted());  
            for (int k = 0; k < getNbDemandsToBeRouted(); k++){
                std::ostringstream varName;
                varName << "t";
                varName << "(" + std::to_string(edge) + "," + std::to_string(s) + "," ;
                varName <<  std::to_string(getToBeRouted_k(k).getId() + 1) + ")";
                IloInt upperBound = 1;
                t[edge][s][k] = IloBoolVar(model.getEnv(), 0, upperBound, varName.str().c_str());
                model.add(t[edge][s][k]);
                // std::cout << "Created variable: " << var[d][arc].getName() << std::endl;
            }
        }
    }
    std::cout << "T variables were created." << std::endl;

    if(instance.getInput().getChosenObj() == Input::OBJECTIVE_METRIC_1p){
        for (ListGraph::EdgeIt e(compactGraph); e != INVALID; ++e){
            int edge = getCompactEdgeLabel(e);
            std::string varName = "maxSlice(" + std::to_string(edge) + ")";
            IloInt lowerBound = instance.getPhysicalLinkFromIndex(edge).getMaxUsedSlicePosition();
            IloInt upperBound = instance.getPhysicalLinkFromIndex(edge).getNbSlices();
            maxSlicePerLink[edge] = IloIntVar(model.getEnv(), lowerBound, upperBound, varName.c_str());
            model.add(maxSlicePerLink[edge]);
        }
    }

    if(instance.getInput().getChosenObj() == Input::OBJECTIVE_METRIC_8){
        std::string varName = "maxSlice";
        IloInt lowerBound = instance.getMaxUsedSlicePosition();
        IloInt upperBound = instance.getMaxSlice();
        maxSliceOverall = IloIntVar(model.getEnv(), lowerBound, upperBound, varName.c_str()); 
        model.add(maxSliceOverall);
    }
}

/****************************************************************************************/
/*									Objective Function    								*/
/****************************************************************************************/

/* Set the objective Function */
void YoussoufForm::setObjective(){
    IloExpr objective = getObjFunction();
    model.add(IloMinimize(model.getEnv(), objective));
    objective.end();
}

/* Returns the objective function expression. */
IloExpr YoussoufForm::getObjFunction(){
    IloExpr obj(model.getEnv());
    
    if(instance.getInput().getChosenObj() == Input::OBJECTIVE_METRIC_0){
        IloInt CONSTANT = 0;
        obj += CONSTANT;
        return obj;
    }
    /*
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
    */

    for (int k = 0; k < getNbDemandsToBeRouted(); k++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(0).getNbSlices(); s++){
            //int arc = getArcIndex(a, d);
            //double coeff = getCoeff(a, d);
            //coeff += (instance.getInput().getInitialLagrangianMultiplier() * getArcLength(a, 0) );
            obj += z[s][k];
        }
    }
    return obj;
}

/****************************************************************************************/
/*									Base Constraints    								*/
/****************************************************************************************/

/* Defines Origin constraints. At most one arc leaves each node and exactly one arc leaves the source. */
void YoussoufForm::setOriginConstraints(){
    for (int k = 0; k < getNbDemandsToBeRouted(); k++){  
        IloRange originConstraint = getOriginConstraint_k(k);
        model.add(originConstraint);
    }
}

/* Returns the source constraint associated with a demand and a node. */
IloRange YoussoufForm::getOriginConstraint_k(int k){
    IloExpr exp(model.getEnv());
    IloInt upperBound = 1;
    IloInt lowerBound = 1;
    int origin = getToBeRouted_k(k).getSource();
    for (ListGraph::EdgeIt e(compactGraph); e != INVALID; ++e){
        int labelU = getCompactNodeLabel(compactGraph.u(e));
        int labelV = getCompactNodeLabel(compactGraph.v(e));
        if (labelU == origin || labelV == origin){
            int edge = getCompactEdgeLabel(e);
            exp += x[edge][k];
        }
    }
    std::ostringstream constraintName;
    constraintName << "Origin(" << getToBeRouted_k(k).getId()+1 << ")";
    IloRange constraint(model.getEnv(), lowerBound, exp, upperBound, constraintName.str().c_str());
    exp.end();
    return constraint;
}

/* Defines Path Continuity constraints. If an edge enters a node, then an arc must leave it. *
void YoussoufForm::setPathContinuityConstraints(){
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

/* Returns the flow conservation constraint associated with a demand and a node. 
IloRange FlowForm::getFlowConservationConstraint_i_d(IloBoolVarMatrix &var, IloModel &mod, ListDigraph::Node &v, const Demand & demand, int d){
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

/* Defines Destination constraints. Exactly one arc enters the target. */
void YoussoufForm::setDestinationConstraints(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){   
        IloRange destinationConstraint = getDestinationConstraint_k(d);
        model.add(destinationConstraint);
    }
}

/* Returns the destination constraint associated with a demand. */
IloRange YoussoufForm::getDestinationConstraint_k(int k){
    IloExpr exp(model.getEnv());
    IloInt upperBound = 1;
    IloInt lowerBound = 1;
    int destination = getToBeRouted_k(k).getTarget();
    for (ListGraph::EdgeIt e(compactGraph); e != INVALID; ++e){
        int labelU = getCompactNodeLabel(compactGraph.u(e));
        int labelV = getCompactNodeLabel(compactGraph.v(e));
        if (labelU == destination || labelV == destination){
            int edge = getCompactEdgeLabel(e);
            exp += x[edge][k];
        }
    }
    std::ostringstream constraintName;
    constraintName << "Destination(" << getToBeRouted_k(k).getId()+1 << ")";
    IloRange constraint(model.getEnv(), lowerBound, exp, upperBound, constraintName.str().c_str());
    exp.end();
    return constraint;
}

/* Defines Degree constraints. At most two edges are adjacent to any node. */
void YoussoufForm::setDegreeConstraints(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){  
        for (ListGraph::NodeIt v(compactGraph); v != INVALID; ++v){
            int node = getCompactNodeLabel(v);
            if (node != getToBeRouted_k(d).getSource() && node != getToBeRouted_k(d).getTarget()){
                IloRange degreeConstraint = getDegreeConstraint_k(d, v);
                model.add(degreeConstraint);
            }
        }
    }
}

/* Returns the degree constraint associated with a demand k and a node v. */
IloRange YoussoufForm::getDegreeConstraint_k(int k, ListGraph::Node &v){
    IloExpr exp(model.getEnv());
    IloInt upperBound = 2;
    IloInt lowerBound = 0;
    for (ListGraph::IncEdgeIt e(compactGraph, v); e != INVALID; ++e){
        int edge = getCompactEdgeLabel(e);
        exp += x[edge][k];
    }
    std::ostringstream constraintName;
    constraintName << "Degree(" << getToBeRouted_k(k).getId()+1 << "," << getCompactNodeLabel(v)+1 << ")";
    IloRange constraint(model.getEnv(), lowerBound, exp, upperBound, constraintName.str().c_str());
    exp.end();
    return constraint;
}

/* Defines Transmission-Reach constraints. Demands must be routed within a length limit. */
void YoussoufForm::setTransmissionReachConstraints(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){   
        IloRange transmissionReach = getTransmissionReachConstraint_k(d);
        model.add(transmissionReach);
    }
}

/* Returns the Transmission-Reach constraint associated with a demand. */
IloRange YoussoufForm::getTransmissionReachConstraint_k(int k){
    IloExpr exp(model.getEnv());
    IloNum upperBound = getToBeRouted_k(k).getMaxLength();
    IloNum lowerBound = 0;
    for (ListGraph::EdgeIt e(compactGraph); e != INVALID; ++e){
        int edge = getCompactEdgeLabel(e);
        double coeff = getCompactLength(e);
        exp += coeff*x[edge][k];
    }
    std::ostringstream constraintName;
    constraintName << "TrasmissionReach(" << getToBeRouted_k(k).getId()+1 << ")";
    IloRange constraint(model.getEnv(), lowerBound, exp, upperBound, constraintName.str().c_str());
    exp.end();
    return constraint;
}


/* Defines Channel-Selection constraints. Only one slot can be the last assigned to a demand. */
void YoussoufForm::setChannelSelectionConstraints(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){   
        IloRange channelSelection = getChannelSelectionConstraint_k(d);
        model.add(channelSelection);
    }
}

/* Returns the channel-selection constraint associated with a demand. */
IloRange YoussoufForm::getChannelSelectionConstraint_k(int k){
    IloExpr exp(model.getEnv());
    IloNum upperBound = 1;
    IloNum lowerBound = 1;
    int load_k = getToBeRouted_k(k).getLoad();
    for (int s = load_k-1; s < instance.getPhysicalLinkFromIndex(0).getNbSlices(); s++){
        exp += z[s][k];
    }
    std::ostringstream constraintName;
    constraintName << "ChannelSelection(" << getToBeRouted_k(k).getId()+1 << ")";
    IloRange constraint(model.getEnv(), lowerBound, exp, upperBound, constraintName.str().c_str());
    exp.end();
    return constraint;
}

/* Defines Forbidden-Slot constraints. */
void YoussoufForm::setForbiddenSlotConstraints(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){   
        IloRange forbiddenSlot = getForbiddenSlotConstraint_k(d);
        model.add(forbiddenSlot);
    }
}

/* Returns the Forbidden-Slot constraint associated with a demand. */
IloRange YoussoufForm::getForbiddenSlotConstraint_k(int k){
    IloExpr exp(model.getEnv());
    IloNum upperBound = 0;
    IloNum lowerBound = 0;
    int load_k = getToBeRouted_k(k).getLoad();
    for (int s = 0; s < load_k-1; s++){
        exp += z[s][k];
    }
    std::ostringstream constraintName;
    constraintName << "ForbiddenSlot(" << getToBeRouted_k(k).getId()+1 << ")";
    IloRange constraint(model.getEnv(), lowerBound, exp, upperBound, constraintName.str().c_str());
    exp.end();
    return constraint;
}

/* Defines Edge-Slot constraints. */
void YoussoufForm::setEdgeSlotConstraints(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListGraph::EdgeIt e(compactGraph); e != INVALID; ++e){
            int edge = getCompactEdgeLabel(e);
            IloRange edgeSlot = getEdgeSlotConstraint_k_e(d, edge);
            model.add(edgeSlot);
        }
    }
}

/* Returns the Edge-Slot constraint associated with a demand and an edge. */
IloRange YoussoufForm::getEdgeSlotConstraint_k_e(int k, int e){
    IloExpr exp(model.getEnv());
    IloNum upperBound = 0;
    IloNum lowerBound = 0;
    int load_k = getToBeRouted_k(k).getLoad();
    for (int s = 0; s < instance.getPhysicalLinkFromIndex(e).getNbSlices(); s++){
        exp += t[e][s][k];
    }
    exp -= load_k*x[e][k];

    std::ostringstream constraintName;
    constraintName << "EdgeSlot(" << getToBeRouted_k(k).getId()+1 << "," << e+1 << ")";
    IloRange constraint(model.getEnv(), lowerBound, exp, upperBound, constraintName.str().c_str());
    exp.end();
    return constraint;
}

/* Defines Demand-Edge-Slot constraints. */
void YoussoufForm::setDemandEdgeSlotConstraints(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListGraph::EdgeIt e(compactGraph); e != INVALID; ++e){
            int edge = getCompactEdgeLabel(e);
            for (int s = 0; s < instance.getPhysicalLinkFromIndex(edge).getNbSlices(); s++){
                IloRange demandEdgeSlot = getDemandEdgeSlotConstraint_k_e_s(d, edge, s);
                model.add(demandEdgeSlot); 
            }
        }
    }
}

/* Returns the Demand-Edge-Slot constraint associated with a demand, an edge and a slice. */
IloRange YoussoufForm::getDemandEdgeSlotConstraint_k_e_s(int k, int e, int s){
    IloExpr exp(model.getEnv());
    IloNum upperBound = 1;
    IloNum lowerBound = 0;
    int load_k = getToBeRouted_k(k).getLoad();
    int maxS = s + load_k - 1;
    if (instance.getPhysicalLinkFromIndex(e).getNbSlices() - 1 < s + load_k - 1){
        maxS = instance.getPhysicalLinkFromIndex(e).getNbSlices() - 1;
    }
    exp += x[e][k];
    for (int slice = s; slice <= maxS; slice++){
        exp += z[slice][k];
    }
    exp -= t[e][s][k];
    
    std::ostringstream constraintName;
    constraintName << "DemandEdgeSlot(" << getToBeRouted_k(k).getId()+1 << "," << e+1 << "," << s+1 << ")";
    std::cout << "Creating " << constraintName.str() << std::endl;
    IloRange constraint(model.getEnv(), lowerBound, exp, upperBound, constraintName.str().c_str());
    exp.end();
    return constraint;
}

/* Defines Non-Overlapping constraints. Demands must not overlap eachother's slices. */ 
void YoussoufForm::setNonOverlappingConstraints(){
    for (ListGraph::EdgeIt e(compactGraph); e != INVALID; ++e){
        int edge = getCompactEdgeLabel(e);
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(edge).getNbSlices(); s++){
            IloRange nonOverlap = getNonOverlappingConstraint_e_s(edge, s);
            model.add(nonOverlap);
        }
    }
}

/* Returns the non-overlapping constraint associated with an edge and a slice. */
IloRange YoussoufForm::getNonOverlappingConstraint_e_s(int e, int s){
    IloExpr exp(model.getEnv());
    IloNum upperBound = 1;
    IloNum lowerBound = 0;
    for (int k = 0; k < getNbDemandsToBeRouted(); k++){
        exp += t[e][s][k];
    }
    std::ostringstream constraintName;
    constraintName << "NonOverlap(" << e+1 << "," << s+1 << ")";
    IloRange constraint(model.getEnv(), lowerBound, exp, upperBound, constraintName.str().c_str());
    exp.end();
    return constraint;
}

/****************************************************************************************/
/*						Objective function related constraints    						*/
/****************************************************************************************/

/* Defines the Link's Max Used Slice Position constraints. The max used slice position on each link must be greater than every slice position used in the link. 
void FlowForm::setMaxUsedSlicePerLinkConstraints(IloBoolVarMatrix &var, IloIntVarArray &maxSlicePerLink, IloModel &mod){
    for (int i = 0; i < instance.getNbEdges(); i++){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            IloRange maxUsedSlicePerLinkConst = getMaxUsedSlicePerLinkConstraints(var, maxSlicePerLink, i, d, mod);
            mod.add(maxUsedSlicePerLinkConst);
        }
    }
}

IloRange FlowForm::getMaxUsedSlicePerLinkConstraints(IloBoolVarMatrix &var, IloIntVarArray &maxSlicePerLink, int linkIndex, int d, IloModel &mod){
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

/* Defines the Overall Max Used Slice Position constraints. 
void FlowForm::setMaxUsedSliceOverallConstraints(IloBoolVarMatrix &var, IloIntVar maxSliceOverall, IloModel &mod){
    for (int i = 0; i < instance.getNbEdges(); i++){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            IloRange maxUsedSliceOverallConst = getMaxUsedSliceOverallConstraints(var, maxSliceOverall, i, d, mod);
            mod.add(maxUsedSliceOverallConst);
        }
    }
}


IloRange FlowForm::getMaxUsedSliceOverallConstraints(IloBoolVarMatrix &var, IloIntVar &maxSlice, int linkIndex, int d, IloModel &mod){
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


/****************************************************************************************/
/*                                    Methods    	                 					*/
/****************************************************************************************/

/* Recovers the obtained MIP solution and builds a path for each demand on its associated graph from RSA. */
void YoussoufForm::updatePath(){
    for(int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            int edge = getArcLabel(a, d);
            int slice = getArcSlice(a, d);
            if ((cplex.getValue(x[edge][d]) >= 0.9) && (cplex.getValue(z[slice][d]) >= 0.9)){
                (*vecOnPath[d])[a] = getToBeRouted_k(d).getId();
            }
            else{
                (*vecOnPath[d])[a] = -1;
            }
        }
    }
}

/* Displays the value of each variable in the obtained solution. 
void FlowForm::displayVariableValues(){
    for(int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            int arc = getArcIndex(a, d);
            std::cout << x[d][arc].getName() << " = " << cplex.getValue(x[d][arc]) << "   ";
        }
        std::cout << std::endl;
    }
}

/** Displays the obtained paths. 
void FlowForm::displayOnPath(){
    for(int d = 0; d < getNbDemandsToBeRouted(); d++){ 
        std::cout << "For demand " << getToBeRouted_k(d).getId() + 1 << " : " << std::endl;
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            if ((*vecOnPath[d])[a] == getToBeRouted_k(d).getId()){
                displayArc(d, a);
            }
        }
    }
}

/****************************************************************************************/
/*										Destructor										*/
/****************************************************************************************/

/* Destructor. Clears the vectors of demands and links. */
YoussoufForm::~YoussoufForm(){
    x.end();
    z.end();
    t.end();
    maxSlicePerLink.end();
    maxSliceOverall.end();
}