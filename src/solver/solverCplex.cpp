#include "solverCplex.h"


int SolverCplex::count = 0;

/****************************************************************************************/
/*										Constructor										*/
/****************************************************************************************/

/* Constructor. The RSA constructor is called and the arc map storing the index of the preprocessed graphs associated is built. */
SolverCplex::SolverCplex(const Instance &inst) : AbstractSolver(inst, STATUS_UNKNOWN), model(env), cplex(model), obj(env){
    std::cout << "--- CPLEX has been initalized ---" << std::endl;
    setCplexParams(inst.getInput());
    implementFormulation();
    exportFormulation(inst);
    count++;
}

void SolverCplex::updateRSA(Instance &instance){
    if ((cplex.getStatus() == IloAlgorithm::Optimal) || (cplex.getStatus() == IloAlgorithm::Feasible)){   
        formulation->updatePath(getSolution());
        formulation->updateInstance(instance);
        instance.setWasBlocked(false);
        formulation->displayPaths();

    }
    else{
        std::cout << "Decrease the number of demands to be treated." << std::endl;
        instance.decreaseNbDemandsAtOnce();
        instance.setWasBlocked(true);
    }
}


std::vector<double> SolverCplex::getSolution(){
    std::vector<double> solution;
    solution.resize(var.getSize());
    for (unsigned int i = 0; i < solution.size(); i++){
        solution[i] = cplex.getValue(var[i]);
    }
    return solution;
}

void SolverCplex::solve(){
    CplexCallback myGenericCallback(var, formulation);
    CPXLONG contextMask = 0;
    contextMask |= IloCplex::Callback::Context::Id::Candidate;
    contextMask |= IloCplex::Callback::Context::Id::Relaxation;
    cplex.use(&myGenericCallback, contextMask);
    
    IloNum timeStart = cplex.getCplexTime();
    std::cout << "Solving..." << std::endl;
    std::vector<ObjectiveFunction> myObjectives = formulation->getObjectiveSet();
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
    if ((cplex.getStatus() == IloAlgorithm::Optimal) || (cplex.getStatus() == IloAlgorithm::Feasible)){    
        std::cout << "Status: " << cplex.getStatus() << std::endl;
        std::cout << "Objective Function Value: " << cplex.getObjValue() << std::endl;
        displaySolution();
    }
    else{
        std::cout << "Could not find a feasible solution..." << std::endl;
    }
}


AbstractSolver::Status SolverCplex::getStatus(){
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

void SolverCplex::setCplexParams(const Input &input){
    cplex.setParam(IloCplex::Param::MIP::Display, 4);
    cplex.setParam(IloCplex::Param::TimeLimit, input.getIterationTimeLimit());
    
    std::cout << "CPLEX parameters have been defined..." << std::endl;
}

void SolverCplex::implementFormulation(){
    setVariables(formulation->getVariables());
    setConstraints(formulation->getConstraints());
    setObjective(formulation->getObjFunction(0));
    formulation->clearConstraints();
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
        obj = IloMinimize(model.getEnv(), exp, myObjective.getName().c_str());
        break;
    case ObjectiveFunction::DIRECTION_MAX:
        obj = IloMaximize(model.getEnv(), exp, myObjective.getName().c_str());
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



/** Displays the value of each variable in the obtained solution. **/
void SolverCplex::displaySolution(){
    for (IloInt i = 0; i < var.getSize(); i++){
        if (cplex.getValue(var[i]) >= 1 - EPS){
            std::cout << var[i].getName() << " = " << cplex.getValue(var[i]) << std::endl;
        }
    }
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

CplexCallback::CplexCallback(const IloNumVarArray _var, AbstractFormulation* &_formulation): var(_var){ 
    formulation = _formulation;
}

void CplexCallback::invoke (const IloCplex::Callback::Context &context){
    if ( context.inRelaxation() ) {
        addUserCuts(context);
    }
    if ( context.inCandidate() ){
        addLazyConstraints(context);
    }
}

void CplexCallback::addUserCuts(const IloCplex::Callback::Context &context) const{
    
    try {
        Constraint constraint = formulation->solveSeparationProblemFract(getFractionalSolution(context));
        if (constraint.getSize() > 0){
            //std::cout << "A violated cut was found: ";
            if (constraint.getLb() != -INFTY){
                //std::cout << "A >= violated cut was found: ";
                context.addUserCut(to_IloExpr(context, constraint.getExpression()) >= constraint.getLb(), IloCplex::UseCutPurge, IloFalse);
            }
            if (constraint.getUb() != INFTY){
                //std::cout << "A <= violated cut was found: ";
                context.addUserCut(to_IloExpr(context, constraint.getExpression()) <= constraint.getUb(), IloCplex::UseCutPurge, IloFalse);
            }
        }
    }
    catch (...) {
        throw;
    }
}

void CplexCallback::addLazyConstraints(const IloCplex::Callback::Context &context) const{
    
    if ( !context.isCandidatePoint() ){
        throw IloCplex::Exception(-1, "Unbounded solution");
    }
    try {
        Constraint constraint = formulation->solveSeparationProblemInt(getIntegerSolution(context));
        if (constraint.getSize() > 0){
            //std::cout << "A lazy constraint was found:";
            //constraint.display();
            IloRange cut(context.getEnv(), constraint.getLb(), to_IloExpr(context, constraint.getExpression()), constraint.getUb());
            //std::cout << cut << std::endl;
            context.rejectCandidate(cut);
        }
        else{
            //std::cout << "The solution is valid." << std::endl;
        }

    }
    catch (...) {
        throw;
    }
    
}

std::vector<double> CplexCallback::getIntegerSolution(const IloCplex::Callback::Context &context) const{
    const int NB_VAR = var.getSize();
    std::vector<double> solution(NB_VAR);
    for (int i = 0; i < NB_VAR; i++){
        solution[i] = context.getCandidatePoint(var[i]);
    }
    return solution;
}

std::vector<double> CplexCallback::getFractionalSolution(const IloCplex::Callback::Context &context) const{
    const int NB_VAR = var.getSize();
    std::vector<double> solution(NB_VAR);
    for (int i = 0; i < NB_VAR; i++){
        solution[i] = context.getRelaxationPoint(var[i]);
    }
    return solution;
}


IloExpr CplexCallback::to_IloExpr(const IloCplex::Callback::Context &context, const Expression &e) const{
    IloExpr exp(context.getEnv());
    for (int i = 0; i < e.getNbTerms(); i++){
        int index = e.getTerm_i(i).getVar().getId();
        double coefficient = e.getTerm_i(i).getCoeff();
        exp += coefficient*var[index];
    }
    return exp;
}

// Destructor
CplexCallback::~CplexCallback(){}