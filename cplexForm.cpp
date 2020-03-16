#include "cplexForm.h"


CplexForm::CplexForm(Instance &instance, const Demand &demand) : ShortestPath(instance, demand), model(env), cplex(model), x(env, countArcs(g)){
    
	/************************************************/
	/*				    SET VARIABLES				*/
	/************************************************/
    std::cout << "--- CPLEX has been initalized ---" << std::endl;
    defineVariables();
    std::cout << "Variables have been defined..." << std::endl;

	/************************************************/
	/*			    SET OBJECTIVE FUNCTION			*/
	/************************************************/
    IloExpr objective = getObjFunction();
    model.add(IloMinimize(env, objective));
    objective.end();
    std::cout << "Objective function has been defined..." << std::endl;

	/************************************************/
	/*			      SET CONSTRAINTS				*/
	/************************************************/
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v){
        IloRange st = getShortestPathConstraint_i(v);
        model.add(st);
    }
    std::cout << "Shortest path constraints have been defined..." << std::endl;

    IloRange lengthConstraint = getLengthConstraint(demand);
    model.add(lengthConstraint);
    //lengthConstraint.end();
    std::cout << "Length constraint has been defined..." << std::endl;

    
	/************************************************/
	/*		    EXPORT LINEAR PROGRAM TO .LP		*/
	/************************************************/
    std::string file = "model.lp";
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
        updatePath();
        displayOnPath();
    }

	/************************************************/
	/*		            UPDATE MAPPING        		*/
	/************************************************/
    if (cplex.getStatus() == IloAlgorithm::Optimal){
        updateInstance(instance, demand);
    }
    else{
        std::cout << "Last topology:" << std::endl;
		instance.displayDetailedTopology();
        std::cout << "Could not find a path!" << std::endl;
        exit(0);
    }
}

void CplexForm::defineVariables(){
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a){
        int id = arcId[a];
        std::ostringstream varName;
        varName << "x(";
        varName << "(" + std::to_string(nodeLabel[g.source(a)]) + "," + std::to_string(nodeSlice[g.source(a)]) + ")";
        varName << "_";
        varName << "(" + std::to_string(nodeLabel[g.target(a)]) + "," + std::to_string(nodeSlice[g.target(a)]) + ")";
        varName << ")";
        x[id] = IloNumVar(env, 0.0, 1.0, ILOINT, varName.str().c_str());
        model.add(x[id]);
        //std::cout << "Creating variable " << id << ": " << varName.str().c_str() << std::endl;
    }
}

IloExpr CplexForm::getObjFunction(){
    IloExpr obj(env);
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a){
        ListDigraph::Node s = g.source(a);
        if(nodeLabel[s] == SOURCE){
            int id = arcId[a];
            int coeff = arcSlice[a]; 
            obj += coeff*x[id];
        }
        else{
            int id = arcId[a];
            int coeff = 1; 
            obj += coeff*x[id];
        }
    }
    return obj;
}

// Flow constraints. Everything that enters must go out. Only 1 leaves the Source and only 1 enters the Target
IloRange CplexForm::getShortestPathConstraint_i(ListDigraph::Node v){
    IloExpr exp(env);
    IloInt rhs = 0;
    for (ListDigraph::OutArcIt a(g, v); a != INVALID; ++a){
        int id = arcId[a];
        exp += x[id];
    }
    for (ListDigraph::InArcIt a(g, v); a != INVALID; ++a){
        int id = arcId[a];
        exp += (-1)*x[id];
    }
    if (nodeLabel[v] == SOURCE){
        rhs = 1;
    }
    if (nodeLabel[v] == TARGET){
        rhs = -1;
    }
    std::ostringstream constraintName;
    constraintName << "Flow(" << nodeLabel[v] << "," << nodeSlice[v] << ")";
    IloRange constraint(env, rhs, exp, rhs, constraintName.str().c_str());
    exp.end();
    return constraint;
}


IloRange CplexForm::getLengthConstraint(const Demand &demand){
    IloExpr exp(env);
    double rhs = demand.getMaxLength();
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a){
        int id = arcId[a];
        double coeff = length[a];
        exp += coeff*x[id];
    }
    std::ostringstream constraintName;
    constraintName << "Length";
    IloRange constraint(env, -IloInfinity, exp, rhs, constraintName.str().c_str());
    exp.end();
    return constraint;
}

void CplexForm::updatePath(){
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a){
        int id = arcId[a];
        if (cplex.getValue(x[id]) >= 0.9){
            onPath[a] = true;
        }
        else{
            onPath[a] = false;
        }
    }
}

void CplexForm::displayOnPath(){
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a){
        if (onPath[a]){
            displayArc(a);
        }
    }
}
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