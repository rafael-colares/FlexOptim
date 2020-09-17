#include "solverCplex.h"


int SolverCplex::count = 0;

/****************************************************************************************/
/*										Constructor										*/
/****************************************************************************************/

/* Constructor. The RSA constructor is called and the arc map storing the index of the preprocessed graphs associated is built. */
SolverCplex::SolverCplex(Instance &inst) : Solver(STATUS_UNKNOWN), model(env), cplex(model), obj(env){
    std::cout << "--- CPLEX has been initalized ---" << std::endl;
    setCplexParams(inst);
    switch (inst.getInput().getChosenFormulation()){
        case Input::FORMULATION_FLOW:
        {
            FlowForm myFormulation(inst);
            implementFormulation(myFormulation.getVariables(), myFormulation.getConstraints(), myFormulation.getObjFunction(0));
            solve(myFormulation.getObjectiveSet());
            if ((cplex.getStatus() == IloAlgorithm::Optimal) || (cplex.getStatus() == IloAlgorithm::Feasible)){    
                std::cout << "Status: " << cplex.getStatus() << std::endl;
				std::cout << "Objective Function Value: " << cplex.getObjValue() << std::endl;
                myFormulation.updatePath(getSolution());
                myFormulation.displayOnPath();
                myFormulation.updateInstance(inst);
                // update instance.
            }
            else{
                std::cout << "Could not find a feasible solution..." << std::endl;
                std::cout << "Decrease the number of demands to be treated." << std::endl;
                inst.decreaseNbDemandsAtOnce();
                inst.setWasBlocked(true);
            }
            break;
        }
        case Input::FORMULATION_EDGE_NODE:
        {
            std::cout << "WARNING: Formulation edge_node has not been implemented yet." << std::endl; /** @todo implement edge_node **/
            break;
        }
        default:
        {
            std::cout<< "ERROR: Formulation " << inst.getInput().getChosenFormulation() << " is invalid." << std::endl;
            exit(0);
        }
    }
    //exportFormulation(inst);
    count++;
}




std::vector<double> SolverCplex::getSolution(){
    std::vector<double> solution;
    solution.resize(var.getSize());
    for (unsigned int i = 0; i < solution.size(); i++){
        solution[i] = cplex.getValue(var[i]);
    }
    return solution;
}

void SolverCplex::solve(const std::vector<ObjectiveFunction> &myObjectives){
    IloNum timeStart = cplex.getCplexTime();
    std::cout << "Solving..." << std::endl;
    for (unsigned int i = 0; i < myObjectives.size(); i++){
        if (i >= 1){
            model.remove(obj);
            setObjective(myObjectives[i]);
        }
        
        std::cout << "Chosen objective: " << myObjectives[i].getName() << std::endl;
        cplex.solve();
        
        if ((cplex.getStatus() == IloAlgorithm::Optimal) || (cplex.getStatus() == IloAlgorithm::Feasible)){
            IloNum objValue = cplex.getObjValue();
            std::cout << "Objective Function Value: " << objValue << std::endl;
            if (i < myObjectives.size() - 1){
                IloExpr objectiveExpression = to_IloExpr(myObjectives[i].getExpression());
                IloRange constraint(model.getEnv(), objValue, objectiveExpression, objValue);
                std::cout << "Add constraint: " << objectiveExpression << " = " << objValue << std::endl;
                model.add(constraint);
                objectiveExpression.end();
            }
        }
        else{
            // Stop optimizing.
            std::cout << "Could not find a feasible solution..." << std::endl;
            i = myObjectives.size()+1;
        }
    }
    IloNum timeFinish = cplex.getCplexTime();
    std::cout << "Optimization done in " << timeFinish - timeStart << " secs." << std::endl;
}


Solver::Status SolverCplex::getStatus(){
    setStatus(STATUS_ERROR);

    if (cplex.getStatus() == IloAlgorithm::Unknown) {
        setStatus(STATUS_UNKNOWN);
    }
    if (cplex.getStatus() == IloAlgorithm::Feasible) {
        setStatus(STATUS_FEASIBLE);
    }
    if (cplex.getStatus() == IloAlgorithm::Optimal) {
        setStatus(STATUS_OPTIMAL);
    }
    if (cplex.getStatus() == IloAlgorithm::Infeasible) {
        setStatus(STATUS_INFEASIBLE);
    }
    if (cplex.getStatus() == IloAlgorithm::Unbounded) {
        setStatus(STATUS_UNBOUNDED);
    }
    if (cplex.getStatus() == IloAlgorithm::InfeasibleOrUnbounded) {
        setStatus(STATUS_INFEASIBLE_OR_UNBOUNDED);
    }
    
    if (currentStatus == STATUS_ERROR){
		std::cout << "Got an status error." << std::endl;
		exit(0);
	}
    return currentStatus;
}

/* Returns the total number of CPLEX default cuts applied during optimization. */
IloInt SolverCplex::getNbCutsFromCplex(){
    IloInt cutsFromCplex = 0;
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutCover);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutGubCover);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutFlowCover);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutClique);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutFrac);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutMir);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutFlowPath);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutDisj); 
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutImplBd);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutZeroHalf);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutMCF);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutLocalCover);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutTighten);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutObjDisj);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutLiftProj);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutLocalImplBd);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutBQP);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutRLT);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutBenders);
    return cutsFromCplex;
}


void SolverCplex::exportFormulation(const Instance &instance){
    std::string file = instance.getInput().getOutputPath() + "LP/model" + std::to_string(count) + ".lp";
    cplex.exportModel(file.c_str());
    std::cout << "LP model has been exported." << std::endl;
}

void SolverCplex::setCplexParams(const Instance &instance){
    cplex.setParam(IloCplex::Param::MIP::Display, 4);
    cplex.setParam(IloCplex::Param::TimeLimit, instance.getInput().getIterationTimeLimit());
    std::cout << "CPLEX parameters have been defined..." << std::endl;
}

void SolverCplex::implementFormulation(const std::vector<Variable> &vars, const std::vector<Constraint> &constraints, const ObjectiveFunction &myObjective){
    setVariables(vars);
    setConstraints(constraints);
    setObjective(myObjective);
}

IloExpr SolverCplex::to_IloExpr(const Expression &e){
    IloExpr exp(model.getEnv());
    for (unsigned int i = 0; i < e.getTerms().size(); i++){
        int index = e.getTerm_i(i).getVar().getId();
        double coefficient = e.getTerm_i(i).getCoeff();
        exp += coefficient*var[index];
    }
    return exp;
}

/** Defines the objective function. **/
void SolverCplex::setObjective(const ObjectiveFunction &myObjective){
    
    IloExpr exp(model.getEnv());
    exp = to_IloExpr(myObjective.getExpression());
    switch (myObjective.getDirection())
    {
    case ObjectiveFunction::DIRECTION_MIN:
        obj = IloMinimize(model.getEnv(), exp);
        break;
    case ObjectiveFunction::DIRECTION_MAX:
        obj = IloMaximize(model.getEnv(), exp);
        break;
    default:
        std::cout << "ERROR: Invalid direction of objective function." << std::endl;
        exit(0);
        break;
    }
    
    model.add(obj);
    exp.end();
}

/* Defines the constraints needed in the MIP formulation. */
void SolverCplex::setConstraints(const std::vector<Constraint> &myConstraints){
    for (unsigned int i = 0; i < myConstraints.size(); i++){ 
        IloExpr exp(model.getEnv());
        for (unsigned int j = 0; j < myConstraints[i].getExpression().getTerms().size(); j++){
            int index = myConstraints[i].getExpression().getTerm_i(j).getVar().getId();
            double coefficient = myConstraints[i].getExpression().getTerm_i(j).getCoeff();
            exp += coefficient*var[index];
        }
        IloRange constraint(model.getEnv(), myConstraints[i].getLb(), exp, myConstraints[i].getUb(), myConstraints[i].getName().c_str());
        model.add(constraint);
        exp.end();
    }
}

void SolverCplex::setVariables(const std::vector<Variable> &myVars){
    var = IloNumVarArray(model.getEnv(), myVars.size());
    for (unsigned int i = 0; i < myVars.size(); i++){ 
        int pos = myVars[i].getId();
        switch (myVars[i].getType())
        {
        case Variable::TYPE_BOOLEAN:
            var[pos] = IloBoolVar(model.getEnv(), myVars[i].getLb(), myVars[i].getUb(), myVars[i].getName().c_str());
            break;
        
        case Variable::TYPE_INTEGER:
            var[pos] = IloIntVar(model.getEnv(), myVars[i].getLb(), myVars[i].getUb(), myVars[i].getName().c_str());
            break;
        
        case Variable::TYPE_REAL:
            var[pos] = IloFloatVar(model.getEnv(), myVars[i].getLb(), myVars[i].getUb(), myVars[i].getName().c_str());
            break;
        
        default:
            std::cout << "ERROR: Variable type has not been recognized." << std::endl;
            exit(0);
            break;
        }
        model.add(var[pos]);
            // std::cout << "Created variable: " << var[d][arc].getName() << std::endl;
    }
    std::cout << "CPLEX variables have been defined." << std::endl;
}



/****************************************************************************************/
/*										Destructor										*/
/****************************************************************************************/
/* Destructor. Free solver memory. */
SolverCplex::~SolverCplex(){
    obj.end();
    var.end();
    cplex.end();
    model.end();
    env.end();
}