#include "flowForm.h"


/* Constructor. Builds the Online RSA mixed-integer program and solves it using CPLEX. */
FlowForm::FlowForm(const Instance &inst) : AbstractFormulation(inst){
    std::cout << "--- Flow formulation has been chosen ---" << std::endl;
    this->setVariables();
    this->setConstraints();
    this->setObjectives();
    std::cout << "--- Flow formulation has been defined ---" << std::endl;
}



/****************************************************************************************/
/*										Variables    									*/
/****************************************************************************************/

/* Define variables x[d][a] for every arc a in the extedend graph #d. */
void FlowForm::setVariables(){
    x.resize(getNbDemandsToBeRouted());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){ 
        x[d].resize(countArcs(*vecGraph[d]));  
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
            int upperBound = 1;
            if (instance.hasEnoughSpace(label, slice, getToBeRouted_k(d)) == false){
                upperBound = 0;
                std::cout << "STILL REMOVING VARIABLES IN FORMULATION. \n" ;
            }
            int varId = getNbVar();
            x[d][arc] = Variable(varId, 0, upperBound, Variable::TYPE_BOOLEAN, 0, varName.str());
            incNbVar();
            // std::cout << "Created variable: " << var[d][arc].getName() << std::endl;
        }
    }
    std::cout << "Flow variables have been defined..." << std::endl;

    maxSlicePerLink.resize(instance.getNbEdges());
    for (int i = 0; i < instance.getNbEdges(); i++){
        std::string varName = "maxSlice(" + std::to_string(instance.getPhysicalLinkFromIndex(i).getId() + 1) + ")";
        int lowerBound = instance.getPhysicalLinkFromIndex(i).getMaxUsedSlicePosition();
        int upperBound = instance.getPhysicalLinkFromIndex(i).getNbSlices();
        int varId = getNbVar();
        maxSlicePerLink[i] = Variable(varId, lowerBound, upperBound, Variable::TYPE_REAL, 0, varName);
        incNbVar();
    }
    std::cout << "Max slice variables have been defined..." << std::endl;

    std::string varName = "maxSliceOverall";
    int lowerBound = instance.getMaxUsedSlicePosition();
    int upperBound = instance.getMaxSlice();
    int varId = getNbVar();
    maxSliceOverall = Variable(varId, lowerBound, upperBound, Variable::TYPE_REAL, 0, varName);
    incNbVar();
    std::cout << "Max slice overall variable has been defined..." << std::endl;
}

VarArray FlowForm::getVariables(){
    VarArray vec;
    vec.resize(getNbVar());
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){ 
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            int arc = getArcIndex(a, d);
            int pos = x[d][arc].getId();
            vec[pos] = x[d][arc];
        }
    }
    for (int i = 0; i < instance.getNbEdges(); i++){
        int pos = maxSlicePerLink[i].getId();
        vec[pos] = maxSlicePerLink[i];
    }
    
    int pos = maxSliceOverall.getId();
    vec[pos] = maxSliceOverall;
    return vec;
}


void FlowForm::setVariableValues(const std::vector<double> &vals){
    
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){ 
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            int arc = getArcIndex(a, d);
            int pos = x[d][arc].getId();
            double newValue = vals[pos];
            x[d][arc].setVal(newValue);
        }
    }
    for (int i = 0; i < instance.getNbEdges(); i++){
        int pos = maxSlicePerLink[i].getId();
        double newValue = vals[pos];
        maxSlicePerLink[i].setVal(newValue);
    }
    
    int pos = maxSliceOverall.getId();
    double newValue = vals[pos];
    maxSliceOverall.setVal(newValue);
}
/****************************************************************************************/
/*									Objective Function    								*/
/****************************************************************************************/

/* Set the objective Function */
void FlowForm::setObjectives(){
    std::vector<Input::ObjectiveMetric> chosenObjectives = instance.getInput().getChosenObj();
    objectiveSet.resize(chosenObjectives.size());
    for (unsigned int i = 0; i < chosenObjectives.size(); i++){
        Expression myObjective = this->getObjFunctionFromMetric(chosenObjectives[i]);
        objectiveSet[i].setExpression(myObjective);
        objectiveSet[i].setDirection(ObjectiveFunction::DIRECTION_MIN);
        objectiveSet[i].setName(instance.getInput().getObjName(chosenObjectives[i]));
        std::cout << "Objective " << chosenObjectives[i] << " has been defined." << std::endl;
    }
}

/* Returns the objective function expression. */
Expression FlowForm::getObjFunctionFromMetric(Input::ObjectiveMetric chosenObjective){
    Expression obj;
    switch (chosenObjective)
    {
    case Input::OBJECTIVE_METRIC_0:
        {
        break;
        }
    case Input::OBJECTIVE_METRIC_1:
        {
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
                int arc = getArcIndex(a, d);
                double coeff = getCoeffObj1(a, d);
                Term term(x[d][arc], coeff);
                obj.addTerm(term);
            }
        }
        break;
        }
    case Input::OBJECTIVE_METRIC_1p:
        {
        for (int i = 0; i < instance.getNbEdges(); i++){
            Term term(maxSlicePerLink[i], 1);
            obj.addTerm(term);
        }
        break;
        }
    case Input::OBJECTIVE_METRIC_2:
        {
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
                int arc = getArcIndex(a, d);
                double coeff = getCoeffObj2(a, d);
                Term term(x[d][arc], coeff);
                obj.addTerm(term);
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
                Term term(x[d][arc], coeff);
                obj.addTerm(term);
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
                Term term(x[d][arc], coeff);
                obj.addTerm(term);
            }
        }
        break;
        }
    case Input::OBJECTIVE_METRIC_8:
        {
        Term term(maxSliceOverall, 1);
        obj.addTerm(term);
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

void FlowForm::setConstraints(){
    this->setSourceConstraints();
    this->setFlowConservationConstraints();
    this->setTargetConstraints();
    this->setLengthConstraints();
    this->setNonOverlappingConstraints();    

    this->setMaxUsedSlicePerLinkConstraints();    
    this->setMaxUsedSliceOverallConstraints();    
    this->setMaxUsedSliceOverallConstraints2();    
    this->setMaxUsedSliceOverallConstraints3();   
} 

/* Defines Source constraints. At most one arc leaves each node and exactly one arc leaves the source. */

void FlowForm::setSourceConstraints(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){  
        for (int label = 0; label < instance.getNbNodes(); label++){
            Constraint sourceConstraint = getSourceConstraint_d_n(getToBeRouted_k(d), d, label);
            constraintSet.push_back(sourceConstraint);
        }
    }
    
    std::cout << "Source constraints have been defined..." << std::endl;
}

/* Returns the source constraint associated with a demand and a node. */

Constraint FlowForm::getSourceConstraint_d_n(const Demand & demand, int d, int nodeLabel){
    Expression exp;
    double upperBound = 1;
    double lowerBound = 0;
    for (ListDigraph::NodeIt v(*vecGraph[d]); v != INVALID; ++v){
        if (getNodeLabel(v, d) == nodeLabel){
            for (ListDigraph::OutArcIt a((*vecGraph[d]), v); a != INVALID; ++a){
                int arc = getArcIndex(a, d); 
                Term term(x[d][arc], 1);
                exp.addTerm(term);
            }
        }
    }
    std::ostringstream constraintName;
    constraintName << "Source(" << nodeLabel+1 << "," << demand.getId()+1 << "):" << constraintSet.size();
    if (nodeLabel == demand.getSource()){
        lowerBound = 1;
    }
    if (nodeLabel == demand.getTarget()){
        upperBound = 0;
    }
    Constraint constraint(lowerBound, exp, upperBound, constraintName.str());
    return constraint;
}

/* Defines Flow Conservation constraints. If an arc enters a node, then an arc must leave it. */
void FlowForm::setFlowConservationConstraints(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){   
        for (ListDigraph::NodeIt v(*vecGraph[d]); v != INVALID; ++v){
            int label = getNodeLabel(v, d);
            if( (label != getToBeRouted_k(d).getSource()) && (label != getToBeRouted_k(d).getTarget()) ){
                Constraint flow = getFlowConservationConstraint_i_d(v, getToBeRouted_k(d), d);
                constraintSet.push_back(flow);
            }
        }
    }
    std::cout << "Flow conservation constraints have been defined..." << std::endl;
}

/* Returns the flow conservation constraint associated with a demand and a node. */
Constraint FlowForm::getFlowConservationConstraint_i_d(ListDigraph::Node &v, const Demand & demand, int d){
    Expression exp;
    int rhs = 0;
    for (ListDigraph::OutArcIt a((*vecGraph[d]), v); a != INVALID; ++a){
        int arc = getArcIndex(a, d);
        Term term(x[d][arc], 1);
        exp.addTerm(term);
    }
    for (ListDigraph::InArcIt a((*vecGraph[d]), v); a != INVALID; ++a){
        int arc = getArcIndex(a, d); 
        Term term(x[d][arc], -1);
        exp.addTerm(term);
    }
    std::ostringstream constraintName;
    int label = getNodeLabel(v, d);
    int slice = getNodeSlice(v, d);
    constraintName << "Flow(" << label+1 << "," << slice+1 << "," << demand.getId()+1 << ")";
    Constraint constraint(rhs, exp, rhs, constraintName.str());
    return constraint;
}

/* Defines Target constraints. Exactly one arc enters the target. */
void FlowForm::setTargetConstraints(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){   
        Constraint targetConstraint = getTargetConstraint_d(getToBeRouted_k(d), d);
        constraintSet.push_back(targetConstraint);
    }
    std::cout << "Target constraints have been defined..." << std::endl;
}

/* Returns the target constraint associated with a demand. */
Constraint FlowForm::getTargetConstraint_d(const Demand & demand, int d){
    Expression exp;
    int rhs = 1;
    for (ListDigraph::NodeIt v(*vecGraph[d]); v != INVALID; ++v){
        int label = getNodeLabel(v, d);
        if (label == demand.getTarget()){
            for (ListDigraph::InArcIt a((*vecGraph[d]), v); a != INVALID; ++a){
                int arc = getArcIndex(a, d); 
                Term term(x[d][arc], 1);
                exp.addTerm(term);
            }
        }
    }
    std::ostringstream constraintName;
    constraintName << "Target(" << demand.getId()+1 << ")";
    Constraint constraint(rhs, exp, rhs, constraintName.str());
    return constraint;
}

/* Defines Length constraints. Demands must be routed within a length limit. */
void FlowForm::setLengthConstraints(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){   
        Constraint lengthConstraint = getLengthConstraint(getToBeRouted_k(d), d);
        constraintSet.push_back(lengthConstraint);
    }
    std::cout << "Length constraints have been defined..." << std::endl;
}

/* Returns the length constraint associated with a demand. */
Constraint FlowForm::getLengthConstraint(const Demand &demand, int d){
    Expression exp;
    double rhs = demand.getMaxLength();
    int source = demand.getSource();
    int hop = instance.getInput().getHopPenalty();
    for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
        int arc = getArcIndex(a, d); 
        double coeff = getArcLength(a, d);
        int tail = getNodeLabel((*vecGraph[d]).source(a), d);
        if (tail != source){
            coeff += hop;
        }
        Term term(x[d][arc], coeff);
        exp.addTerm(term);
    }
    std::ostringstream constraintName;
    constraintName << "Length(" << demand.getId()+1 << ")";
    Constraint constraint(0, exp, rhs, constraintName.str());
    return constraint;
}


/* Defines the second set of Improved Non-Overlapping constraints. */
void FlowForm::setNonOverlappingConstraints(){
    for (int i = 0; i < instance.getNbEdges(); i++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(i).getNbSlices(); s++){
            Constraint nonOverlap = getNonOverlappingConstraint(instance.getPhysicalLinkFromIndex(i).getId(), s);
            constraintSet.push_back(nonOverlap);
        }
    }
    std::cout << "Non-Overlapping constraints has been defined..." << std::endl;
}

/* Returns the non-overlapping constraint associated with an edge and a slice */
Constraint FlowForm::getNonOverlappingConstraint(int linkLabel, int slice){
	Expression exp;
    int rhs = 1;
    
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        int demandLoad = getToBeRouted_k(d).getLoad();
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            if( (getArcLabel(a, d) == linkLabel) && (getArcSlice(a, d) >= slice)  && (getArcSlice(a, d) <= slice + demandLoad - 1) ){
                int index = getArcIndex(a, d);
                Term term(x[d][index], 1);
                exp.addTerm(term);
            }
        }
    }
    
    std::ostringstream constraintName;
    constraintName << "NonOverlap(" << linkLabel+1 << "," << slice+1 << ")";
    Constraint constraint(0, exp, rhs, constraintName.str());
    return constraint;
}



/****************************************************************************************/
/*						Objective function related constraints    						*/
/****************************************************************************************/

/* Defines the Link's Max Used Slice Position constraints. The max used slice position on each link must be greater than every slice position used in the link. */
void FlowForm::setMaxUsedSlicePerLinkConstraints(){
    for (int i = 0; i < instance.getNbEdges(); i++){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){
            Constraint maxUsedSlicePerLinkConst = getMaxUsedSlicePerLinkConstraints(i, d);
            constraintSet.push_back(maxUsedSlicePerLinkConst);
        }
    }
    std::cout << "Max Used Slice Per Link constraints have been defined..." << std::endl;
}

Constraint FlowForm::getMaxUsedSlicePerLinkConstraints(int linkIndex, int d){
    Expression exp;
    int rhs = 0;
    int linkLabel = instance.getPhysicalLinkFromIndex(linkIndex).getId();
    for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
        if (getArcLabel(a, d) == linkLabel){
            int index = getArcIndex(a, d);
            int slice = getArcSlice(a, d);
            Term term(x[d][index], slice);
            exp.addTerm(term);
        }
    }
    Term term(maxSlicePerLink[linkIndex], -1);
    exp.addTerm(term);
    
    std::ostringstream constraintName;
    constraintName << "MaxUsedSlicePerLink(" << linkLabel+1 << "," << getToBeRouted_k(d).getId()+1 << ")";
    Constraint constraint(-INFTY, exp, rhs, constraintName.str());
    return constraint;
}

/* Defines the Overall Max Used Slice Position constraints. */
void FlowForm::setMaxUsedSliceOverallConstraints(){
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        Constraint maxUsedSliceOverallConst = getMaxUsedSliceOverallConstraints(d);
        constraintSet.push_back(maxUsedSliceOverallConst);
    }
    std::cout << "Max Used Slice Overall constraints have been defined..." << std::endl;
}

Constraint FlowForm::getMaxUsedSliceOverallConstraints(int d){
    Expression exp;
    int rhs = 0;
    for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
        if (getToBeRouted_k(d).getSource() == getNodeLabel((*vecGraph[d]).source(a), d)){
            int index = getArcIndex(a, d);
            int slice = getArcSlice(a, d);
            Term term(x[d][index], slice);
            exp.addTerm(term);
        }
    }
    Term term(maxSliceOverall, -1);
    exp.addTerm(term);

    std::ostringstream constraintName;
    constraintName << "MaxUsedSliceOverall(" << getToBeRouted_k(d).getId()+1 << ")";
    Constraint constraint(-INFTY, exp, rhs, constraintName.str());
    return constraint;
}

/* Defines the Overall Max Used Slice Position constraints. */
void FlowForm::setMaxUsedSliceOverallConstraints2(){
    for (int i = 0; i < instance.getNbEdges(); i++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(i).getNbSlices(); s++){
            Constraint maxUsedSliceOverallConst = getMaxUsedSliceOverallConstraints2(instance.getPhysicalLinkFromIndex(i).getId(), s);
            constraintSet.push_back(maxUsedSliceOverallConst);
        }
    }
    std::cout << "Max Used Slice Overall2 constraints have been defined..." << std::endl;
}



Constraint FlowForm::getMaxUsedSliceOverallConstraints2(int linkLabel, int s){
    Expression exp;
    int rhs = 0;
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        int demandLoad = getToBeRouted_k(d).getLoad();
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            if ((getArcLabel(a, d) == linkLabel) && (getArcSlice(a, d) >= s) && (getArcSlice(a, d) <= s + demandLoad - 1)){
                int index = getArcIndex(a, d);
                int slice = getArcSlice(a, d);
                Term term(x[d][index], slice);
                exp.addTerm(term);
            }
        }
    }
    Term term(maxSliceOverall, -1);
    exp.addTerm(term);

    std::ostringstream constraintName;
    constraintName << "MaxUsedSliceOverall2(" << linkLabel+1 << "," << s+1 << ")";
    Constraint constraint(-INFTY, exp, rhs, constraintName.str());
    return constraint;
}


/* Defines the Overall Max Used Slice Position constraints 3. */
void FlowForm::setMaxUsedSliceOverallConstraints3(){
    for (int i = 0; i < instance.getNbNodes(); i++){
        for (int s = 0; s < instance.getPhysicalLinkFromIndex(0).getNbSlices(); s++){
            Constraint maxUsedSliceOverallConst = getMaxUsedSliceOverallConstraints3(i, s);
            constraintSet.push_back(maxUsedSliceOverallConst);
        }
    }
    std::cout << "Max Used Slice Overall3 constraints have been defined..." << std::endl;
}



Constraint FlowForm::getMaxUsedSliceOverallConstraints3(int nodeLabel, int s){
    Expression exp;
    int rhs = 0;
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
                        Term term(x[d][index], coeff);
                        exp.addTerm(term);
                    }
                }
            }
        }
    }
    Term term(maxSliceOverall, -1);
    exp.addTerm(term);
    
    std::ostringstream constraintName;
    constraintName << "MaxUsedSliceOverall2(" << nodeLabel+1 << "," << s+1 << ")";
    Constraint constraint(-INFTY, exp, rhs, constraintName.str());
    return constraint;
}

/** Returns a vector of node id's corresponding to the sequence of nodes that the d-th demand passes through. **/
std::vector<int> FlowForm::getPathNodeSequence(int d){
    std::vector<int> path;
    int origin = getToBeRouted_k(d).getSource();
    int destination = getToBeRouted_k(d).getTarget();
    ListDigraph::Node SOURCE = getFirstNodeFromLabel(d, origin);
    ListDigraph::Node TARGET = getFirstNodeFromLabel(d, destination);
    if (TARGET == INVALID || SOURCE == INVALID){
        std::cout << "ERROR: Could not find source or target from demand " << getToBeRouted_k(d).getId() + 1 << "." << std::endl;
        exit(0);
    }
    path.push_back(getNodeLabel(SOURCE, d));
    //std::cout << "Search the path from origin to destination." << std::endl;
    ListDigraph::Node currentNode = SOURCE;
    while (currentNode != TARGET){
        ListDigraph::Arc nextArc = INVALID;
        for (ListDigraph::OutArcIt a(*vecGraph[d], currentNode); a != INVALID; ++a){
            int arc = getArcIndex(a, d);
            if (x[d][arc].getVal() >= 1 - EPS){
                nextArc = a;
            }
        }
        if (nextArc == INVALID){
            std::cout << "ERROR: Could not find path continuity.." << std::endl;
            exit(0);
        }
        currentNode = (*vecGraph[d]).target(nextArc);
        path.push_back(getNodeLabel(currentNode, d));
    }
    return path;
}

/* Recovers the obtained MIP solution and builds a path for each demand on its associated graph from RSA. */
void FlowForm::updatePath(const std::vector<double> &vals){
    
    setVariableValues(vals);
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
            std::cout << "ERROR: Could not find source or target from demand." << std::endl;
            exit(0);
        }
        
        //std::cout << "Search the path from origin to destination." << std::endl;
        ListDigraph::Node currentNode = TARGET;
        while (currentNode != SOURCE){
            ListDigraph::Arc previousArc = INVALID;
            for (ListDigraph::InArcIt a(*vecGraph[d], currentNode); a != INVALID; ++a){
                int arc = getArcIndex(a, d);
                if (x[d][arc].getVal() >= 1 - EPS){
                    (*vecOnPath[d])[a] = getToBeRouted_k(d).getId();
                    previousArc = a;
                }
            }
            if (previousArc == INVALID){
                std::cout << "ERROR: Could not find path continuity.." << std::endl;
                exit(0);
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
            std::cout << x[d][arc].getName() << " = " << x[d][arc].getVal() << "   ";
        }
        std::cout << std::endl;
    }
}

std::vector<Constraint> FlowForm::solveSeparationProblemInt(const std::vector<double> &solution){
    //std::cout << "Entering separation problem of an integer point for Flow Form." << std::endl;
    //Gnpy is disabled.
    std::vector<Constraint> cuts;
    //cuts = separationGNPY(solution);
    return cuts;
}

std::vector<Constraint> FlowForm::solveSeparationProblemFract(const std::vector<double> &solution){
    //std::cout << "Entering separation problem of a fractional point for Flow Form." << std::endl;
    return std::vector<Constraint>();
}

std::vector<Constraint> FlowForm::separationGNPY(const std::vector<double> &solution){
    std::vector<Constraint> cuts;
    setVariableValues(solution);
    writeServiceFile();

    system("gnpy-path-request ../oopt-gnpy/gnpy/example-data/topo_spain.json ../oopt-gnpy/gnpy/example-data/service_test.json -e ../oopt-gnpy/gnpy/example-data/eqpt_config_spain.json -o result.json");
    std::ifstream ifs("result.json");
    std::string fileContent((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        
    }
    
    return cuts;
}

void FlowForm::writeServiceFile(){
    std::ofstream serviceFile;
    serviceFile.open ("../oopt-gnpy/gnpy/example-data/service_test.json");
    serviceFile << "{\n";
    
    serviceFile << "\t\"path-request\": [\n";
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){
        writePathRequest(serviceFile, d);
    }
    
    serviceFile << "\t],\n";
    serviceFile << "\t\"synchronization\": [\n";
    serviceFile << "\t]\n";
    serviceFile << "}\n";
}


void FlowForm::writePathRequest(std::ofstream &serviceFile, int d){
    std::string source = std::to_string(getToBeRouted_k(d).getSource()+1) + ".1";
    std::string destination = std::to_string(getToBeRouted_k(d).getTarget()+1) + ".1";
    int load = getToBeRouted_k(d).getLoad();
    std::string mode = "mode_" + std::to_string(load);
    std::vector<int> path = getPathNodeSequence(d);
    serviceFile << "\t{\n";
    serviceFile << "\t\t" << "\"request-id\": \"" << std::to_string(getToBeRouted_k(d).getId()+1) << "\",\n";
    serviceFile << "\t\t" << "\"source\": \"" << source << "\",\n";
    serviceFile << "\t\t" << "\"destination\": \"" << destination << "\",\n";
    serviceFile << "\t\t" << "\"src-tp-id\": \"" << source << "\",\n";
    serviceFile << "\t\t" << "\"dst-tp-id\": \"" << destination << "\",\n";
    serviceFile << "\t\t" << "\"bidirectional\": false,\n";
    serviceFile << "\t\t" << "\"path-constraints\": {\n";
    serviceFile << "\t\t\t" << "\"te-bandwidth\": {\n";
    serviceFile << "\t\t\t\t" << "\"technology\": \"flexi-grid\",\n";
    serviceFile << "\t\t\t\t" << "\"trx_type\": \"Voyager\",\n";
    serviceFile << "\t\t\t\t" << "\"trx_mode\": \"" << mode << "\",\n";
    serviceFile << "\t\t\t\t" << "\"effective-freq-slot\": [\n";
    serviceFile << "\t\t\t\t\t" << "{\n";
    serviceFile << "\t\t\t\t\t\t" << "\"N\": \"null\",\n";
    serviceFile << "\t\t\t\t\t\t" << "\"M\": \"null\"\n";
    serviceFile << "\t\t\t\t\t" << "}\n";
    serviceFile << "\t\t\t\t" << "],\n";
    serviceFile << "\t\t\t\t" << "\"spacing\":" << std::to_string(12.5*load) << "e9,\n";
    serviceFile << "\t\t\t\t" << "\"max-nb-of-channel\": null,\n";
    serviceFile << "\t\t\t\t" << "\"output-power\": null,\n";
    serviceFile << "\t\t\t\t" << "\"path_bandwidth\": 400e9\n";
    serviceFile << "\t\t\t" << "}\n";
    serviceFile << "\t\t" << "},\n";
    serviceFile << "\t\t" << "\"explicit-route-objects\": {\n";
    serviceFile << "\t\t\t" << "\"route-object-include-exclude\": [\n";
    for (unsigned int i = 0; i < path.size(); i++){
        serviceFile << "\t\t\t\t" << "{\n";
        serviceFile << "\t\t\t\t\t" << "\"explicit-route-usage\": \"route-include-ero\",\n";
        serviceFile << "\t\t\t\t\t" << "\"index\": " << std::to_string(i) << ",\n";
        serviceFile << "\t\t\t\t\t" << "\"num-unnum-hop\": {\n";
        serviceFile << "\t\t\t\t\t\t" << "\"node-id\": \"" << std::to_string(path[i]+1) << "\",\n";
        serviceFile << "\t\t\t\t\t\t" << "\"link-tp-id\": \"link-tp-id is not used\",\n";
        serviceFile << "\t\t\t\t\t\t" << "\"hop-type\": \"STRICT\"\n";
        serviceFile << "\t\t\t\t\t" << "}\n";
        if (i < path.size()-1){
            serviceFile << "\t\t\t\t" << "},\n";
        }
        else{
            serviceFile << "\t\t\t\t" << "}\n";
        }
    }
    serviceFile << "\t\t\t" << "]\n";
    serviceFile << "\t\t" << "}\n";
    if (d < getNbDemandsToBeRouted()-1){
        serviceFile << "\t" << "},\n";
    }
    else{
        serviceFile << "\t" << "}\n";
    }
}
/****************************************************************************************/
/*										Destructor										*/
/****************************************************************************************/

/* Destructor. Clears the vectors of demands and links. */
FlowForm::~FlowForm(){
    x.clear();
    maxSlicePerLink.clear();
}