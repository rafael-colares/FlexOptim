#include "cplexForm.h"


CplexForm::CplexForm(const Instance &inst) : RSA(inst), model(env), cplex(model), x(env, countArcs(g)){
    
    this->setToBeRouted(inst.getNextDemands());
    displayToBeRouted();
    
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
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){  
        for (int i = 0; i < instance.getNbNodes(); i++){
            // std::cout << "Creating source constraint " << d << " ..." << std::endl;
            IloRange sourceConstraint = getSourceConstraint_d(getToBeRouted()[d], d, i);
            //std::cout << "Created source constraint " << d << " ..." << std::endl;
            model.add(sourceConstraint);
        } 
    }
    std::cout << "Source constraints have been defined..." << std::endl;

    for (int d = 0; d < getNbDemandsToBeRouted(); d++){   
        IloRange targetConstraint = getTargetConstraint_d(getToBeRouted()[d], d);
        model.add(targetConstraint);
    }
    std::cout << "Target constraints have been defined..." << std::endl;

    for (int d = 0; d < getNbDemandsToBeRouted(); d++){   
        for (ListDigraph::NodeIt v(g); v != INVALID; ++v){
            if( (nodeLabel[v] != getToBeRouted()[d].getSource()) && (nodeLabel[v] != getToBeRouted()[d].getTarget()) ){
                IloRange st = getShortestPathConstraint_i_d(v, getToBeRouted()[d], d);
                model.add(st);
            }
        }
    }
    std::cout << "Shortest path constraints have been defined..." << std::endl;
    
    for (int d = 0; d < getNbDemandsToBeRouted(); d++){   
        IloRange lengthConstraint = getLengthConstraint(getToBeRouted()[d], d);
        model.add(lengthConstraint);
    }
    std::cout << "Length constraint has been defined..." << std::endl;

    
    //must be done: x_a + x_{bar(a)} ,= 1
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a){
        for (int d = 0; d < getNbDemandsToBeRouted(); d++){   
            IloRange subcycle = getSubcycleConstraint(a, getToBeRouted()[d], d);
            model.add(subcycle);
       }
    }
    
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
        std::cout << "Objective Function Value: " << cplex.getObjValue() << std::endl;
        updatePath();
        displayOnPath();
    }
    else{
        std::cout << "Could not find a path!" << std::endl;
        exit(0);
    }
}

void CplexForm::defineVariables(){
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a){
        int arc = arcId[a];
        int nbDemandsToBeRouted = getNbDemandsToBeRouted();
        x[arc] = IloNumVarArray(env, nbDemandsToBeRouted);
        for (int d = 0; d < nbDemandsToBeRouted; d++){    
            std::ostringstream varName;
            varName << "x";
            varName << "(" + std::to_string(nodeLabel[g.source(a)]+1) + "," + std::to_string(nodeLabel[g.target(a)]+1);
            varName << "," + std::to_string(arcSlice[a]+1) + "," + std::to_string(getToBeRouted()[d].getId()+1);
            varName << ")";
            int linkLabel = arcLabel[a];
            int linkSlice = arcSlice[a];
            IloNum upperBound = 1.0;
            if (instance.isRoutable(linkLabel, linkSlice, getToBeRouted()[d]) == false){
                upperBound = 0.0;
            }
            x[arc][d] = IloNumVar(env, 0.0, upperBound, ILOINT, varName.str().c_str());
            model.add(x[arc][d]);
        }
        //std::cout << "Creating variable " << id << ": " << varName.str().c_str() << std::endl;
    }
}

IloExpr CplexForm::getObjFunction(){
    IloExpr obj(env);
    int nbDemandsToBeRouted = getNbDemandsToBeRouted();
    for (int d = 0; d < nbDemandsToBeRouted; d++){
        for (ListDigraph::ArcIt a(g); a != INVALID; ++a){
            ListDigraph::Node s = g.source(a);
            if(nodeLabel[s] == getToBeRouted()[d].getSource()){
                int arc = arcId[a];
                int coeff = arcSlice[a]+1; 
                obj += coeff*x[arc][d];
            }
            else{
                int arc = arcId[a];
                int coeff = 0; 
                obj += coeff*x[arc][d];
            }
        }
    }
    return obj;
}

// Flow constraints. At most 1 leaves each node. Exactly 1 leaves the Source.
IloRange CplexForm::getSourceConstraint_d(const Demand & demand, int d, int i){
    IloExpr exp(env);
    IloInt upperBound = 1;
    IloInt lowerBound = 0;
    std::cout << "Source of demand " << demand.getId()+1 << ": " << demand.getSource()+1 << std::endl;
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v){
        if (nodeLabel[v] == i){
            for (ListDigraph::OutArcIt a(g, v); a != INVALID; ++a){
                int arc = arcId[a];
                // std::cout << "Add arc " << arc << "." << std::endl;
                exp += x[arc][d];
                // std::cout << "Added arc " << arc << "." << std::endl;
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
    IloRange constraint(env, lowerBound, exp, upperBound, constraintName.str().c_str());
    exp.end();
    return constraint;
}
// Flow constraints. Only 1 enters the Target
IloRange CplexForm::getTargetConstraint_d(const Demand & demand, int d){
    IloExpr exp(env);
    IloInt rhs = 1;
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v){
        int node = nodeLabel[v];
        if (node == demand.getTarget()){
            for (ListDigraph::InArcIt a(g, v); a != INVALID; ++a){
                int arc = arcId[a];
                exp += x[arc][d];
            }
        }
    }
    std::ostringstream constraintName;
    constraintName << "Target(" << demand.getId()+1 << ")";
    IloRange constraint(env, rhs, exp, rhs, constraintName.str().c_str());
    exp.end();
    return constraint;
}

// Flow constraints. Everything that enters must go out. 
IloRange CplexForm::getShortestPathConstraint_i_d(ListDigraph::Node v, const Demand & demand, int d){
    IloExpr exp(env);
    IloInt rhs = 0;
    for (ListDigraph::OutArcIt a(g, v); a != INVALID; ++a){
        int arc = arcId[a];
        exp += x[arc][d];
    }
    for (ListDigraph::InArcIt a(g, v); a != INVALID; ++a){
        int arc = arcId[a];
        exp += (-1)*x[arc][d];
    }
    std::ostringstream constraintName;
    constraintName << "Flow(" << nodeLabel[v]+1 << "," << nodeSlice[v]+1 << "," << demand.getId()+1 << ")";
    IloRange constraint(env, rhs, exp, rhs, constraintName.str().c_str());
    exp.end();
    return constraint;
}


IloRange CplexForm::getSubcycleConstraint(const ListDigraph::Arc &arc, const Demand & demand, int d){
    IloExpr exp(env);
    IloNum rhs = 1;
    int label = arcLabel[arc];
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a){
        if(arcLabel[a] == label){
            int id = arcId[a];
            exp += x[id][d];
        }
    }
    std::ostringstream constraintName;
    constraintName << "Subcycle(" << arcId[arc] << "," << demand.getId()+1 << ")";
    IloRange constraint(env, -IloInfinity, exp, rhs, constraintName.str().c_str());
    exp.end();
    return constraint;
}

IloRange CplexForm::getLengthConstraint(const Demand &demand, int d){
    IloExpr exp(env);
    double rhs = demand.getMaxLength();
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a){
        int arc = arcId[a];
        double coeff = length[a];
        exp += coeff*x[arc][d];
    }
    std::ostringstream constraintName;
    constraintName << "Length(" << demand.getId()+1 << ")";
    IloRange constraint(env, -IloInfinity, exp, rhs, constraintName.str().c_str());
    exp.end();
    return constraint;
}

void CplexForm::updatePath(){
    for(int d = 0; d < getNbDemandsToBeRouted(); d++){
        for (ListDigraph::ArcIt a(g); a != INVALID; ++a){
            int arc = arcId[a];
            if (cplex.getValue(x[arc][d]) >= 0.9){
                onPath[a] = getToBeRouted()[d].getId();
            }
            else{
                onPath[a] = -1;
            }
        }
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

void CplexForm::displayToBeRouted(){
    std::cout << "--- ROUTING DEMANDS ";
    for (int i = 0; i < getNbDemandsToBeRouted(); i++){
        std::cout << "#" << toBeRouted[i].getId()+1 << " (" << toBeRouted[i].getSource()+1 << ", " << toBeRouted[i].getTarget()+1 << ")";
    }
    std::cout << " ... --- " << std::endl;
	
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