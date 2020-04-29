#include "cplexForm.h"


int CplexForm::count = 0;

/* Constructor. Builds the Online RSA mixed-integer program and solves it using CPLEX. */
CplexForm::CplexForm(const Instance &inst) : Solver(inst), model(env), cplex(model), x(env, getNbDemandsToBeRouted()){
    std::cout << "--- CPLEX has been chosen ---" << std::endl;
    count++;
    /************************************************/
	/*				    SET VARIABLES				*/
	/************************************************/
    this->setVariables(x, model);
    std::cout << "Variables have been defined..." << std::endl;

	/************************************************/
	/*			    SET OBJECTIVE FUNCTION			*/
	/************************************************/
    this->setObjective(x, model);
    std::cout << "Objective function has been defined..." << std::endl;

	/************************************************/
	/*			      SET CONSTRAINTS				*/
	/************************************************/
    this->setSourceConstraints(x, model);
    std::cout << "Source constraints have been defined..." << std::endl;

    this->setFlowConservationConstraints(x, model);
    std::cout << "Flow conservation constraints have been defined..." << std::endl;

    this->setTargetConstraints(x, model);
    std::cout << "Target constraints have been defined..." << std::endl;

    this->setLengthConstraints(x, model);
    std::cout << "Length constraints have been defined..." << std::endl;

    this->setNonOverlappingConstraints(x, model);    
    std::cout << "Non-Overlapping constraints have been defined..." << std::endl;
    
	/************************************************/
	/*		    EXPORT LINEAR PROGRAM TO .LP		*/
	/************************************************/
    std::string file = getInstance().getInput().getOutputPath() + "LP/model" + std::to_string(count) + ".lp";
    //cplex.exportModel(file.c_str());
    std::cout << "LP model has been exported..." << std::endl;
    
	/************************************************/
	/*             DEFINE CPLEX PARAMETERS   		*/
	/************************************************/
    cplex.setParam(IloCplex::Param::MIP::Display, 2);
    std::cout << "CPLEX parameters have been defined..." << std::endl;

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
    if (cplex.getStatus() == IloAlgorithm::Optimal){
        std::cout << "Optimization done in " << timeFinish - timeStart << " secs." << std::endl;
        std::cout << "Objective Function Value: " << cplex.getObjValue() << std::endl;
        //displayVariableValues();
        updatePath();
        displayOnPath();
        
        std::cout << "Number of cplex cuts: " << getNbCutsFromCplex() << std::endl;
    }
    else{
        std::cout << "Could not find a path!" << std::endl;
        //exit(0);
    }
}

/* Returns the total number of CPLEX default cuts applied during optimization. */
IloInt CplexForm::getNbCutsFromCplex(){
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

/* Recovers the obtained MIP solution and builds a path for each demand on its associated graph from RSA. */
void CplexForm::updatePath(){
    for(int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            int arc = getArcIndex(a, d);
            if (cplex.getValue(x[d][arc]) >= 0.9){
                (*vecOnPath[d])[a] = getToBeRouted_k(d).getId();
            }
            else{
                (*vecOnPath[d])[a] = -1;
            }
        }
    }
}

/* Displays the value of each variable in the obtained solution. */
void CplexForm::displayVariableValues(){
    for(int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            int arc = getArcIndex(a, d);
            std::cout << x[d][arc].getName() << " = " << cplex.getValue(x[d][arc]) << "   ";
        }
        std::cout << std::endl;
    }
}

/** Displays the obtained paths. */
void CplexForm::displayOnPath(){
    for(int d = 0; d < getNbDemandsToBeRouted(); d++){ 
        std::cout << "For demand " << getToBeRouted_k(d).getId() + 1 << " : " << std::endl;
        for (ListDigraph::ArcIt a(*vecGraph[d]); a != INVALID; ++a){
            if ((*vecOnPath[d])[a] == getToBeRouted_k(d).getId()){
                displayArc(d, a);
            }
        }
    }
}
