#include "cplexForm.h"

#define EPS 1e-6

/* Constructor. Builds the Online RSA mixed-integer program and solves it using CPLEX. */
FlowForm::FlowForm(const Instance &inst) : Solver(inst), x(env, getNbDemandsToBeRouted()), maxSlicePerLink(env, instance.getNbEdges()), maxSliceOverall(env){
    std::cout << "--- Flow formulation has been chosen ---" << std::endl;
    /************************************************/
	/*				    SET VARIABLES				*/
	/************************************************/
    this->setVariables();
    std::cout << "Variables have been defined..." << std::endl;


	/************************************************/
	/*			      SET CONSTRAINTS				*/
	/************************************************/
    this->setSourceConstraints();
    std::cout << "Source constraints have been defined..." << std::endl;

    this->setFlowConservationConstraints();
    std::cout << "Flow conservation constraints have been defined..." << std::endl;

    this->setTargetConstraints();
    std::cout << "Target constraints have been defined..." << std::endl;

    this->setLengthConstraints();
    std::cout << "Length constraints have been defined..." << std::endl;

    //this->setNonOverlappingConstraints(x, model);    
    //std::cout << "Non-Overlapping constraints have been defined..." << std::endl;

    //this->setImprovedNonOverlappingConstraints_1();    
    std::cout << "First set of Improved Non-Overlapping constraints has been defined..." << std::endl;

    this->setImprovedNonOverlappingConstraints_2();    
    std::cout << "Second set of Improved Non-Overlapping constraints has been defined..." << std::endl;


    this->setMaxUsedSlicePerLinkConstraints();    
    std::cout << "Max Used Slice Per Link constraints have been defined..." << std::endl;

    this->setMaxUsedSliceOverallConstraints();    
    std::cout << "Max Used Slice Overall constraints have been defined..." << std::endl;
    this->setMaxUsedSliceOverallConstraints2();    
    std::cout << "Max Used Slice Overall2 constraints have been defined..." << std::endl;
    this->setMaxUsedSliceOverallConstraints3();    
    std::cout << "Max Used Slice Overall3 constraints have been defined..." << std::endl;

	/************************************************/
    /*             DEFINE OBJECTIVE FUNCTION   		*/
	/************************************************/
    
    Input::ObjectiveMetric chosenObjective = instance.getInput().getChosenObj_k(0);
    IloObjective myObjective(model.getEnv());
    myObjective = IloMinimize(model.getEnv(), this->getObjFunction(chosenObjective));
    model.add(myObjective);
    //this->setObjective(instance.getInput().getChosenObj_k(0));
    //std::cout << "Objective function has been defined..." << std::endl;
    
	/************************************************/
	/*             DEFINE CPLEX PARAMETERS   		*/
	/************************************************/
    cplex.setParam(IloCplex::Param::MIP::Display, 4);
    cplex.setParam(IloCplex::Param::TimeLimit, getInstance().getInput().getIterationTimeLimit());
    std::cout << "CPLEX parameters have been defined..." << std::endl;
    
	/************************************************/
	/*		    EXPORT LINEAR PROGRAM TO .LP		*/
	/************************************************/
    std::string file = getInstance().getInput().getOutputPath() + "LP/model" + std::to_string(count) + ".lp";
    //cplex.exportModel(file.c_str());
    std::cout << "LP model has been exported..." << std::endl;
    

	/************************************************/
	/*		         SOLVE LINEAR PROGRAM   		*/
	/************************************************/
    IloNum timeStart = cplex.getCplexTime();
    std::cout << "Solving..." << std::endl;
    for (unsigned int i = 0; i < instance.getInput().getChosenObj().size(); i++){
        if (i >= 1){
            model.remove(myObjective);
            chosenObjective = instance.getInput().getChosenObj_k(i);
            myObjective = IloMinimize(model.getEnv(), this->getObjFunction(chosenObjective));
            model.add(myObjective);
        }
        
        std::cout << "Chosen objective: " << chosenObjective << std::endl;
        /*
        else{
            for (int d = 0; d < getNbDemandsToBeRouted(); d++){
                for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
                    int arc = getArcIndex(a, d);
                    switch (chosenObjective)
                    {
                    case Input::OBJECTIVE_METRIC_0:
                        {
                        myObjective.setLinearCoef(x[d][arc], 0);
                        break;
                        }
                    case Input::OBJECTIVE_METRIC_1:
                        {
                        double coeff = getCoeffObj1(a, d);
                        myObjective.setLinearCoef(x[d][arc], coeff);
                        break;
                        }
                    case Input::OBJECTIVE_METRIC_1p:
                        {
                        myObjective.setLinearCoef(x[d][arc], 0);
                        break;
                        }
                    case Input::OBJECTIVE_METRIC_2:
                        {
                        double coeff = getCoeffObj2(a, d);
                        myObjective.setLinearCoef(x[d][arc], coeff);
                        break;
                        }
                    case Input::OBJECTIVE_METRIC_2p:
                        {
                        double coeff = getCoeffObj2p(a, d);
                        myObjective.setLinearCoef(x[d][arc], coeff);
                        break;
                        }
                    case Input::OBJECTIVE_METRIC_4:
                        {
                        double coeff = getCoeffObj4(a, d);
                        myObjective.setLinearCoef(x[d][arc], coeff);
                        break;
                        }
                    case Input::OBJECTIVE_METRIC_8:
                        {
                        myObjective.setLinearCoef(x[d][arc], 0);
                        break;
                        }
                    default:
                        {
                        std::cout << "ERROR: Objective '" << chosenObjective << "' is not known." << std::endl;
                        exit(0);
                        break;
                        }
                    }
                }
            }
            for (int i = 0; i < instance.getNbEdges(); i++){
                if (chosenObjective == Input::OBJECTIVE_METRIC_1p){
                    myObjective.setLinearCoef(maxSlicePerLink[i], 1);
                }
                else{
                    myObjective.setLinearCoef(maxSlicePerLink[i], 0);
                }
            }
            
            if (chosenObjective == Input::OBJECTIVE_METRIC_8){
                myObjective.setLinearCoef(maxSliceOverall, 1);
            }
            else{
                myObjective.setLinearCoef(maxSliceOverall, 0);
            }
        }
        */
        //cplex.extract(model);
        cplex.solve();
        
        if ((cplex.getStatus() == IloAlgorithm::Optimal) || (cplex.getStatus() == IloAlgorithm::Feasible)){
            
            std::cout << "Objective Function Value: " << cplex.getObjValue() << std::endl;
            if (i < instance.getInput().getChosenObj().size() - 1){
                IloExpr objectiveConst = getObjFunction(chosenObjective);
                IloRange constraint(model.getEnv(), cplex.getObjValue(), objectiveConst, cplex.getObjValue());
                //std::cout << "Add constraint: " << objectiveConst << " = " << cplex.getObjValue() << std::endl;
                model.add(constraint);
                objectiveConst.end();
            }
        }
        else{
            std::cout << "Could not find a feasible solution..." << std::endl;
            i = instance.getInput().getChosenObj().size();
        }
    }
    IloNum timeFinish = cplex.getCplexTime();

	/************************************************/
	/*		    GET OPTIMAL SOLUTION FOUND        	*/
	/************************************************/
    if ((cplex.getStatus() == IloAlgorithm::Optimal) || (cplex.getStatus() == IloAlgorithm::Feasible)){    
        std::cout << "Solved!" << std::endl;
        std::cout << "Optimization done in " << timeFinish - timeStart << " secs." << std::endl;
        std::cout << "Objective Function Value: " << cplex.getObjValue() << std::endl;
        //displayVariableValues();
        updatePath();
        displayOnPath();
        
        std::cout << "Number of cplex cuts: " << getNbCutsFromCplex() << std::endl;
    }
    else{
        std::cout << "Could not find a feasible solution..." << std::endl;
    }
}




/****************************************************************************************/
/*										Variables    									*/
/****************************************************************************************/

/* Define variables x[d][a] for every arc a in the extedend graph #d. */
void FlowForm::setVariables(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){ 
        x[d] = IloBoolVarArray(model.getEnv(), countArcs(*vecGraph[d]));  
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
            x[d][arc] = IloBoolVar(model.getEnv(), 0, upperBound, varName.str().c_str());
            model.add(x[d][arc]);
            // std::cout << "Created variable: " << var[d][arc].getName() << std::endl;
        }
    }

    for (int i = 0; i < instance.getNbEdges(); i++){
        std::string varName = "maxSlice(" + std::to_string(instance.getPhysicalLinkFromIndex(i).getId() + 1) + ")";
        IloInt lowerBound = instance.getPhysicalLinkFromIndex(i).getMaxUsedSlicePosition();
        IloInt upperBound = instance.getPhysicalLinkFromIndex(i).getNbSlices();
        maxSlicePerLink[i] = IloIntVar(model.getEnv(), lowerBound, upperBound, varName.c_str());
        model.add(maxSlicePerLink[i]);
    }

    std::string varName = "maxSlice";
    IloInt lowerBound = instance.getMaxUsedSlicePosition();
    IloInt upperBound = instance.getMaxSlice();
    maxSliceOverall = IloIntVar(model.getEnv(), lowerBound, upperBound, varName.c_str()); 
    model.add(maxSliceOverall);

}

/****************************************************************************************/
/*									Objective Function    								*/
/****************************************************************************************/

/* Set the objective Function */
void FlowForm::setObjective(Input::ObjectiveMetric chosenObjective){
    IloExpr objective = getObjFunction(chosenObjective);
    model.add(IloMinimize(model.getEnv(), objective));
    objective.end();
}

/* Returns the objective function expression. */
IloExpr FlowForm::getObjFunction(Input::ObjectiveMetric chosenObjective){
    IloExpr obj(model.getEnv());
    switch (chosenObjective)
    {
    case Input::OBJECTIVE_METRIC_0:
        {
        IloInt CONSTANT = 0;
        obj += CONSTANT;
        break;
        }
    case Input::OBJECTIVE_METRIC_1:
        {
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
                int arc = getArcIndex(a, d);
                double coeff = getCoeffObj1(a, d);
                obj += coeff*x[d][arc];
            }
        }
        break;
        }
    case Input::OBJECTIVE_METRIC_1p:
        {
        for (int i = 0; i < instance.getNbEdges(); i++){
            obj += maxSlicePerLink[i];
        }
        break;
        }
    case Input::OBJECTIVE_METRIC_2:
        {
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
                int arc = getArcIndex(a, d);
                double coeff = getCoeffObj2(a, d);
                obj += coeff*x[d][arc];
            }
        }
        break;
        }
    case Input::OBJECTIVE_METRIC_2p:
        {
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
                int arc = getArcIndex(a, d);
                double coeff = getCoeffObj2p(a, d);
                obj += coeff*x[d][arc];
            }
        }
        break;
        }
    case Input::OBJECTIVE_METRIC_4:
        {
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
                int arc = getArcIndex(a, d);
                double coeff = getCoeffObj4(a, d);
                obj += coeff*x[d][arc];
            }
        }
        break;
        }
    case Input::OBJECTIVE_METRIC_8:
        {
        obj += maxSliceOverall;
        break;
        }
    default:
        {
        std::cout << "ERROR: Objective '" << chosenObjective << "' is not known." << std::endl;
        exit(0);
        break;
        }
    }
    return obj;
}

/****************************************************************************************/
/*									Base Constraints    								*/
/****************************************************************************************/

/* Defines Source constraints. At most one arc leaves each node and exactly one arc leaves the source. */
void FlowForm::setSourceConstraints(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){  
        for (ListDigraph::NodeIt v(*vecGraph[d]); v != INVALID; ++v){
            int label = getNodeLabel(v, d);
            IloRange sourceConstraint = getSourceConstraint_d_n(getToBeRouted_k(d), d, label);
            model.add(sourceConstraint);
        } 
    }
}

/* Returns the source constraint associated with a demand and a node. */
IloRange FlowForm::getSourceConstraint_d_n(const Demand & demand, int d, int nodeLabel){
    IloExpr exp(model.getEnv());
    IloInt upperBound = 1;
    IloInt lowerBound = 0;
    for (ListDigraph::NodeIt v(*vecGraph[d]); v != INVALID; ++v){
        if (getNodeLabel(v, d) == nodeLabel){
            for (ListDigraph::OutArcIt a((*vecGraph[d]), v); a != INVALID; ++a){
                int arc = getArcIndex(a, d); 
                exp += x[d][arc];
            }
        }
    }
    std::ostringstream constraintName;
    constraintName << "Source(" << nodeLabel+1 << "," << demand.getId()+1 << ")";
    if (nodeLabel == demand.getSource()){
        lowerBound = 1;
    }
    if (nodeLabel == demand.getTarget()){
        upperBound = 0;
    }
    IloRange constraint(model.getEnv(), lowerBound, exp, upperBound, constraintName.str().c_str());
    exp.end();
    return constraint;
}

/* Defines Flow Conservation constraints. If an arc enters a node, then an arc must leave it. */
void FlowForm::setFlowConservationConstraints(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){   
        for (ListDigraph::NodeIt v(*vecGraph[d]); v != INVALID; ++v){
            int label = getNodeLabel(v, d);
            if( (label != getToBeRouted_k(d).getSource()) && (label != getToBeRouted_k(d).getTarget()) ){
                IloRange st = getFlowConservationConstraint_i_d(v, getToBeRouted_k(d), d);
                model.add(st);
            }
        }
    }
}

/* Returns the flow conservation constraint associated with a demand and a node. */
IloRange FlowForm::getFlowConservationConstraint_i_d(ListDigraph::Node &v, const Demand & demand, int d){
    IloExpr exp(model.getEnv());
    IloInt rhs = 0;
    for (ListDigraph::OutArcIt a((*vecGraph[d]), v); a != INVALID; ++a){
        int arc = getArcIndex(a, d); 
        exp += x[d][arc];
    }
    for (ListDigraph::InArcIt a((*vecGraph[d]), v); a != INVALID; ++a){
        int arc = getArcIndex(a, d); 
        exp += (-1)*x[d][arc];
    }
    std::ostringstream constraintName;
    int label = getNodeLabel(v, d);
    int slice = getNodeSlice(v, d);
    constraintName << "Flow(" << label+1 << "," << slice+1 << "," << demand.getId()+1 << ")";
    IloRange constraint(model.getEnv(), rhs, exp, rhs, constraintName.str().c_str());
    exp.end();
    return constraint;
}

/* Defines Target constraints. Exactly one arc enters the target. */
void FlowForm::setTargetConstraints(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){   
        IloRange targetConstraint = getTargetConstraint_d(getToBeRouted_k(d), d);
        model.add(targetConstraint);
    }
}

/* Returns the target constraint associated with a demand. */
IloRange FlowForm::getTargetConstraint_d(const Demand & demand, int d){
    IloExpr exp(model.getEnv());
    IloInt rhs = 1;
    for (ListDigraph::NodeIt v(*vecGraph[d]); v != INVALID; ++v){
        int label = getNodeLabel(v, d);
        if (label == demand.getTarget()){
            for (ListDigraph::InArcIt a((*vecGraph[d]), v); a != INVALID; ++a){
                int arc = getArcIndex(a, d); 
                exp += x[d][arc];
            }
        }
    }
    std::ostringstream constraintName;
    constraintName << "Target(" << demand.getId()+1 << ")";
    IloRange constraint(model.getEnv(), rhs, exp, rhs, constraintName.str().c_str());
    exp.end();
    return constraint;
}

/* Defines Length constraints. Demands must be routed within a length limit. */
void FlowForm::setLengthConstraints(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){   
        IloRange lengthConstraint = getLengthConstraint(getToBeRouted_k(d), d);
        model.add(lengthConstraint);
    }
}

/* Returns the length constraint associated with a demand. */
IloRange FlowForm::getLengthConstraint(const Demand &demand, int d){
    IloExpr exp(model.getEnv());
    double rhs = demand.getMaxLength();
    for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
        int arc = getArcIndex(a, d); 
        double coeff = getArcLength(a, d);
        exp += coeff*x[d][arc];
    }
    std::ostringstream constraintName;
    constraintName << "Length(" << demand.getId()+1 << ")";
    IloRange constraint(model.getEnv(), -IloInfinity, exp, rhs, constraintName.str().c_str());
    exp.end();
    return constraint;
}

/* Defines Non-Overlapping constraints. Demands must not overlap eachother's slices. */
void FlowForm::setNonOverlappingConstraints(){
    for (int d1 = 0; d1 < getNbDemandsToBeRouted(); d1++){
        for (int d2 = 0; d2 < getNbDemandsToBeRouted(); d2++){
            if(d1 != d2){
                for (int i = 0; i < instance.getNbEdges(); i++){
                    for (int s = 0; s < instance.getPhysicalLinkFromIndex(i).getNbSlices(); s++){
                        IloRange nonOverlap = getNonOverlappingConstraint(instance.getPhysicalLinkFromIndex(i).getId(), s, getToBeRouted_k(d1), d1, getToBeRouted_k(d2), d2);
                        model.add(nonOverlap);
                    }
                }
            }
        }
    }
}

/* Returns the non-overlapping constraint associated with an arc and a pair of demands. */
IloRange FlowForm::getNonOverlappingConstraint(int linkLabel, int slice, const Demand & demand1, int d1, const Demand & demand2, int d2){
    IloExpr exp(model.getEnv());
    IloNum rhs = 1;
    for (ListDigraph::ArcIt a(*vecGraph[d1]); a != INVALID; ++a){
        if( (getArcLabel(a, d1) == linkLabel) && (getArcSlice(a, d1) == slice) ){
            int id = getArcIndex(a, d1);
            exp += x[d1][id];
        }
    }
    for (ListDigraph::ArcIt a(*vecGraph[d2]); a != INVALID; ++a){
        if( (getArcLabel(a, d2) == linkLabel) && (getArcSlice(a, d2) >= slice - demand1.getLoad() + 1) && (getArcSlice(a, d2) <= slice + demand2.getLoad() - 1) ){
            int id = getArcIndex(a, d2);
            exp += x[d2][id];
        }
    }
    std::ostringstream constraintName;
    constraintName << "NonOverlap(" << linkLabel+1 << "," << slice+1 << "," << demand1.getId()+1 << "," << demand2.getId()+1 << ")";
    IloRange constraint(model.getEnv(), -IloInfinity, exp, rhs, constraintName.str().c_str());
    exp.end();
    return constraint;
}

/****************************************************************************************/
/*						Objective function related constraints    						*/
/****************************************************************************************/

/* Defines the Link's Max Used Slice Position constraints. The max used slice position on each link must be greater than every slice position used in the link. */
void FlowForm::setMaxUsedSlicePerLinkConstraints(){
    for (int i = 0; i < instance.getNbEdges(); i++){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            IloRange maxUsedSlicePerLinkConst = getMaxUsedSlicePerLinkConstraints(i, d);
            model.add(maxUsedSlicePerLinkConst);
        }
    }
}

IloRange FlowForm::getMaxUsedSlicePerLinkConstraints(int linkIndex, int d){
    IloExpr exp(model.getEnv());
    IloInt rhs = 0;
    int linkLabel = instance.getPhysicalLinkFromIndex(linkIndex).getId();
    for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
        if (getArcLabel(a, d) == linkLabel){
            int index = getArcIndex(a, d);
            int slice = getArcSlice(a, d);
            exp += slice*x[d][index];
        }
    }
    exp += -maxSlicePerLink[linkIndex];
    
    std::ostringstream constraintName;
    constraintName << "MaxUsedSlicePerLink(" << linkLabel+1 << "," << getToBeRouted_k(d).getId()+1 << ")";
    IloRange constraint(model.getEnv(), -IloInfinity, exp, rhs, constraintName.str().c_str());
    exp.end();
    return constraint;
}

/* Defines the Overall Max Used Slice Position constraints. */
void FlowForm::setMaxUsedSliceOverallConstraints(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        IloRange maxUsedSliceOverallConst = getMaxUsedSliceOverallConstraints(d);
        model.add(maxUsedSliceOverallConst);
    }
}

IloRange FlowForm::getMaxUsedSliceOverallConstraints(int d){
    IloExpr exp(model.getEnv());
    IloInt rhs = 0;
    for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
        if (getToBeRouted_k(d).getSource() == getNodeLabel((*vecGraph[d]).source(a), d)){
            int index = getArcIndex(a, d);
            int slice = getArcSlice(a, d);
            exp += slice*x[d][index];
        }
    }
    exp += -maxSliceOverall;
    
    std::ostringstream constraintName;
    constraintName << "MaxUsedSliceOverall(" << getToBeRouted_k(d).getId()+1 << ")";
    IloRange constraint(model.getEnv(), -IloInfinity, exp, rhs, constraintName.str().c_str());
    exp.end();
    return constraint;
}

/* Defines the Overall Max Used Slice Position constraints. */
void FlowForm::setMaxUsedSliceOverallConstraints2(){
    for (int i = 0; i < instance.getNbEdges(); i++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(i).getNbSlices(); s++){
            IloRange maxUsedSliceOverallConst = getMaxUsedSliceOverallConstraints2(instance.getPhysicalLinkFromIndex(i).getId(), s);
            model.add(maxUsedSliceOverallConst);
        }
    }
}



IloRange FlowForm::getMaxUsedSliceOverallConstraints2(int linkLabel, int s){
    IloExpr exp(model.getEnv());
    IloInt rhs = 0;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        int demandLoad = getToBeRouted_k(d).getLoad();
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            if ((getArcLabel(a, d) == linkLabel) && (getArcSlice(a, d) >= s) && (getArcSlice(a, d) <= s + demandLoad - 1)){
                int index = getArcIndex(a, d);
                int slice = getArcSlice(a, d);
                exp += slice*x[d][index];
            }
        }
    }
    exp += -maxSliceOverall;
    
    std::ostringstream constraintName;
    constraintName << "MaxUsedSliceOverall2(" << linkLabel+1 << "," << s+1 << ")";
    IloRange constraint(model.getEnv(), -IloInfinity, exp, rhs, constraintName.str().c_str());
    exp.end();
    return constraint;
}


/* Defines the Overall Max Used Slice Position constraints 3. */
void FlowForm::setMaxUsedSliceOverallConstraints3(){
    for (int i = 0; i < instance.getNbNodes(); i++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(0).getNbSlices(); s++){
            IloRange maxUsedSliceOverallConst = getMaxUsedSliceOverallConstraints3(i, s);
            model.add(maxUsedSliceOverallConst);
        }
    }
}



IloRange FlowForm::getMaxUsedSliceOverallConstraints3(int nodeLabel, int s){
    IloExpr exp(model.getEnv());
    IloInt rhs = 0;
    int degree = 0;
    for (ListGraph::NodeIt v(compactGraph); v != INVALID; ++v){
        if ((getCompactNodeLabel(v) == nodeLabel)){
            for (ListGraph::IncEdgeIt a(compactGraph, v); a != INVALID; ++a){
                degree++;
            }
        }
    }

    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        int demandLoad = getToBeRouted_k(d).getLoad();
        for (ListDigraph::NodeIt v(*vecGraph[d]); v != INVALID; ++v){
            if ((getNodeLabel(v, d) == nodeLabel)){
                for (ListDigraph::OutArcIt a(*vecGraph[d], v); a != INVALID; ++a){
                    if ( (getArcSlice(a, d) >= s) && (getArcSlice(a, d) <= s + demandLoad - 1) ){

                        int index = getArcIndex(a, d);
                        int slice = getArcSlice(a, d);
                        int coeff = ceil(((double)slice)/((double)degree-1));
                        exp += coeff*x[d][index];
                    }
                }
            }
        }
    }
    exp += -maxSliceOverall;
    
    std::ostringstream constraintName;
    constraintName << "MaxUsedSliceOverall2(" << nodeLabel+1 << "," << s+1 << ")";
    IloRange constraint(model.getEnv(), -IloInfinity, exp, rhs, constraintName.str().c_str());
    exp.end();
    return constraint;
}

/****************************************************************************************/
/*						Improved Non-Overlapping constraints    						*/
/****************************************************************************************/

/* Defines the first set of Improved Non-Overlapping constraints. */
void FlowForm::setImprovedNonOverlappingConstraints_1(){
    for (int k = 0; k < getNbLoadsToBeRouted(); k++){
        int load_k = getLoadsToBeRouted_k(k);
        for (int d2 = 0; d2 < getNbDemandsToBeRouted(); d2++){
            for (int i = 0; i < instance.getNbEdges(); i++){
                for (int s = 0; s < instance.getPhysicalLinkFromIndex(i).getNbSlices(); s++){
                    IloRange improvedNonOverlap1 = getImprovedNonOverlappingConstraint_1(instance.getPhysicalLinkFromIndex(i).getId(), s, load_k, getToBeRouted_k(d2), d2);
                    model.add(improvedNonOverlap1);
                }
            }
        }
    }
}

/* Returns the first improved non-overlapping constraint associated with an arc, a demand and a load. */
IloRange FlowForm::getImprovedNonOverlappingConstraint_1(int linkLabel, int slice, int min_load, const Demand & demand2, int d2){
	IloExpr exp(model.getEnv());
    IloNum rhs = 1;
    
    for (int d1 = 0; d1 < getNbDemandsToBeRouted(); d1++){
        if ((d1 != d2) && (getToBeRouted_k(d1).getLoad() >= min_load)){
            for (ListDigraph::ArcIt a(*vecGraph[d1]); a != INVALID; ++a){
                if( (getArcLabel(a, d1) == linkLabel) && (getArcSlice(a, d1) == slice) ){
                    int index = getArcIndex(a, d1);
                    exp += x[d1][index];
                }
            }
        }
    }
    
    for (ListDigraph::ArcIt a(*vecGraph[d2]); a != INVALID; ++a){
        if( (getArcLabel(a, d2) == linkLabel) && (getArcSlice(a, d2) >= slice - min_load + 1) && (getArcSlice(a, d2) <= slice + demand2.getLoad() - 1) ){
            int index = getArcIndex(a, d2);
            exp += x[d2][index];
        }
    }
    
    std::ostringstream constraintName;
    constraintName << "ImprNonOverlap_1(" << linkLabel+1 << "," << slice+1 << "," << min_load << "," << demand2.getId()+1 << ")";
    IloRange constraint(model.getEnv(), -IloInfinity, exp, rhs, constraintName.str().c_str());
    exp.end();
    return constraint;
}

/* Defines the second set of Improved Non-Overlapping constraints. */
void FlowForm::setImprovedNonOverlappingConstraints_2(){
    for (int i = 0; i < instance.getNbEdges(); i++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(i).getNbSlices(); s++){
            IloRange improvedNonOverlap2 = getImprovedNonOverlappingConstraint_2(instance.getPhysicalLinkFromIndex(i).getId(), s);
            model.add(improvedNonOverlap2);
        }
    }
}

/* Returns the second improved non-overlapping constraint associated with an edge-slice. */
IloRange FlowForm::getImprovedNonOverlappingConstraint_2(int linkLabel, int slice){
	IloExpr exp(model.getEnv());
    IloNum rhs = 1;
    
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        int demandLoad = getToBeRouted_k(d).getLoad();
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            if( (getArcLabel(a, d) == linkLabel) && (getArcSlice(a, d) >= slice)  && (getArcSlice(a, d) <= slice + demandLoad - 1) ){
                int index = getArcIndex(a, d);
                exp += x[d][index];
            }
        }
    }
    
    std::ostringstream constraintName;
    constraintName << "ImprNonOverlap_2(" << linkLabel+1 << "," << slice+1 << ")";
    IloRange constraint(model.getEnv(), -IloInfinity, exp, rhs, constraintName.str().c_str());
    exp.end();
    return constraint;
}



/* Recovers the obtained MIP solution and builds a path for each demand on its associated graph from RSA. */
void FlowForm::updatePath(){
    // Reinitialize OnPath
    //std::cout << "Enter update." << std::endl;
    for(int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
                (*vecOnPath[d])[a] = -1;
        }
    }
    // Fill the mapping with the corresponding demands id.
    for(int d = 0; d < getNbDemandsToBeRouted(); d++){
        int origin = getToBeRouted_k(d).getSource();
        int destination = getToBeRouted_k(d).getTarget();
        ListDigraph::Node SOURCE = INVALID;
        ListDigraph::Node TARGET = INVALID;
        for (ListDigraph::NodeIt v(*vecGraph[d]); v != INVALID; ++v){
            if (getNodeLabel(v, d) == origin){
                SOURCE = v;
            }
            if (getNodeLabel(v, d) == destination){
                TARGET = v;
            }
        }
        if (TARGET == INVALID || SOURCE == INVALID){
            throw IloCplex::Exception(-1, "Could not find source or target from demand.");
        }
        
        //std::cout << "Search the path from origin to destination." << std::endl;
        ListDigraph::Node currentNode = TARGET;
        while (currentNode != SOURCE){
            ListDigraph::Arc previousArc = INVALID;
            for (ListDigraph::InArcIt a(*vecGraph[d], currentNode); a != INVALID; ++a){
                int arc = getArcIndex(a, d);
                if (cplex.getValue(x[d][arc]) >= 1 - EPS){
                    (*vecOnPath[d])[a] = getToBeRouted_k(d).getId();
                    previousArc = a;
                }
            }
            if (previousArc == INVALID){
                throw IloCplex::Exception(-1, "Could not find path continuity.");
            }
            currentNode = (*vecGraph[d]).source(previousArc);
        }
    }
    
    //std::cout << "Leave update." << std::endl;
}

/* Displays the value of each variable in the obtained solution. */
void FlowForm::displayVariableValues(){
    for(int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            int arc = getArcIndex(a, d);
            std::cout << x[d][arc].getName() << " = " << cplex.getValue(x[d][arc]) << "   ";
        }
        std::cout << std::endl;
    }
}

/** Displays the obtained paths. */
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
FlowForm::~FlowForm(){
    x.end();
    maxSlicePerLink.end();
    maxSliceOverall.end();
}