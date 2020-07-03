#include "solver.h"


int Solver::count = 0;

/****************************************************************************************/
/*										Constructor										*/
/****************************************************************************************/

/* Constructor. The RSA constructor is called and the arc map storing the index of the preprocessed graphs associated is built. */
Solver::Solver(const Instance &inst) : RSA(inst), model(env), cplex(model) {
    std::cout << "--- CPLEX has been initalized ---" << std::endl;
    count++;
}


RSA::Status Solver::getStatus(){
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
    
    return currentStatus;
}

/* Returns the total number of CPLEX default cuts applied during optimization. */
IloInt Solver::getNbCutsFromCplex(){
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



/****************************************************************************************/
/*										Destructor										*/
/****************************************************************************************/
/* Destructor. Free solver memory. */
Solver::~Solver(){
    cplex.end();
    model.end();
    env.end();
}