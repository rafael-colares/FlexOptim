#include "cplexForm.h"


int CplexForm::count = 0;

CplexForm::CplexForm(const Instance &inst) : Solver(inst), model(env), cplex(model), x(env, countArcs(g)){
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
    std::string file = inst.getInput().getOutputPath() + "LP/model" + std::to_string(count) + ".lp";
    cplex.exportModel(file.c_str());
    
	/************************************************/
	/*             DEFINE CPLEX PARAMETERS   		*/
	/************************************************/
    cplex.setParam(IloCplex::Param::MIP::Display, 2);

	/************************************************/
	/*		         SOLVE LINEAR PROGRAM   		*/
	/************************************************/
    IloNum timeStart = cplex.getCplexTime();
    cplex.solve();
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
        exit(0);
    }
}

IloInt CplexForm::getNbCutsFromCplex(){
    IloInt cutsFromCplex = 0;
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutClique);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutCover);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutFlowCover);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutGubCover);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutFrac);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutMir);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutFlowPath);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutImplBd);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutDisj); 
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutLocalImplBd);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutZeroHalf);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutMCF);
    cutsFromCplex += getCplex().getNcuts(IloCplex::CutLiftProj);
    return cutsFromCplex;
}
void CplexForm::updatePath(){
    for(int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(g); a != INVALID; ++a){
            int arc = arcId[a];
            if (cplex.getValue(x[arc][d]) >= 0.9){
                onPath[a] = getToBeRouted()[d].getId();
            }
        }
    }
}
void CplexForm::displayVariableValues(){
    for(int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(g); a != INVALID; ++a){
            int arc = arcId[a];
            std::cout << x[arc][d].getName() << " = " << cplex.getValue(x[arc][d]) << "   ";
        }
        std::cout << std::endl;
    }
}

void CplexForm::displayOnPath(){
    for(int d = 0; d < getNbDemandsToBeRouted(); d++){
        std::cout << "For demand " << getToBeRouted()[d].getId() + 1 << " : " << std::endl;
        for (ListDigraph::ArcIt a(g); a != INVALID; ++a){
            if (onPath[a] == getToBeRouted()[d].getId()){
                displayArc(a);
            }
        }
    }
}

/*
void CplexForm::displayPathOnEnv(){
    ListDigraph::Node source = getNode(SOURCE,-1);
    ListDigraph::Arc current = getNextOnPath(source);
    env.out() << "(" << nodeLabel[g.source(current)] + 1 << ", " << nodeSlice[g.source(current)] << ")";
    while(current != INVALID){
        env.out() << " --> (" << nodeLabel[g.target(current)] + 1 << ", " << nodeSlice[g.target(current)] << ")";
        current = getNextOnPath(g.target(current));
    }
    env.out() << "\n";
}
*/