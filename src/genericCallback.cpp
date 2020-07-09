#include "genericCallback.h"

GenericCallback::GenericCallback(const IloBoolVarMatrix &_x, Instance &i, ListGraph &g, std::vector<Demand> &toBeRouted,
         ListGraph::NodeMap<int> &nodeLabelMap, 
         ListGraph::NodeMap<int> &nodeIdMap, 
         ListGraph::EdgeMap<int> &edgeLabelMap, 
         ListGraph::EdgeMap<int> &edgeIdMap): x(_x){ 

    instance = &i;
    graph = & g;
    nodeLabel = & nodeLabelMap;
    nodeId = & nodeIdMap;
    edgeLabel = & edgeLabelMap;
    edgeId = & edgeIdMap;
    demands = & toBeRouted;
}

void GenericCallback::invoke (const IloCplex::Callback::Context &context){
    if ( context.inRelaxation() ) {
        addUserCuts(context);
    }

    if ( context.inCandidate() ){
        addLazyConstraints(context);
    }
}

void GenericCallback::addLazyConstraints(const IloCplex::Callback::Context &context) const{
    std::cout << "Entering callback." << std::endl;
    IloInt const nbEdges = x.getSize();
    IloInt const nbDemands = x[0].getSize();
    if ( !context.isCandidatePoint() ){
        throw IloCplex::Exception(-1, "Unbounded solution");
    }
    for(int d = 0; d < demands->size(); d++){
        
        int origin = (*demands)[d].getSource();
        int destination = (*demands)[d].getTarget();
        ListGraph::Node SOURCE = INVALID;
        ListGraph::Node TARGET = INVALID;
        for (ListGraph::NodeIt v(*graph); v != INVALID; ++v){
            if ((*nodeLabel)[v] == origin){
                SOURCE = v;
            }
            if ((*nodeLabel)[v] == destination){
                TARGET = v;
            }
        }
        if (TARGET == INVALID || SOURCE == INVALID){
            throw IloCplex::Exception(-1, "Could not find source or target from demand inside generic callback.");
        }
        
        ListGraph::Node currentNode = SOURCE;
        std::vector<int> setOfNodes;
        setOfNodes.push_back(origin);
        bool cutFound = false;
        while (currentNode != TARGET && !cutFound){
            ListGraph::Edge nextEdge = INVALID;
            
            for (ListGraph::IncEdgeIt e(*graph, currentNode); e != INVALID; ++e){
                int edge = (*edgeLabel)[e];
                if (edge != (*edgeLabel)[nextEdge]){
                    if (context.getCandidatePoint(x[edge][d]) >= 1 - EPS){
                        setOfNodes.push_back((*nodeLabel)[currentNode]);
                        nextEdge = e;
                    }
                }
            }
            
            if (nextEdge == INVALID){
                // create cut
                std::cout << "Could not find path continuity." << std::endl;
            }
            if((*nodeLabel)[currentNode] == (*nodeLabel)[(*graph).u(nextEdge)]){
                currentNode = (*graph).v(nextEdge);
            }
            else{
                currentNode = (*graph).u(nextEdge);
            }
        }
        
    }
    /*
    IloNum isUsed = context.getCandidatePoint(opened[j]);
    IloNum served = 0.0; // Number of clients currently served from j
    for (IloInt c = 0; c < nbClients; ++c)
        served += context.getCandidatePoint(supply[c][j]);
    if ( served > (nbClients - 1.0) * isUsed + EPS ) {
        IloNumExpr sum = IloExpr(context.getEnv());
        for (IloInt c = 0; c < nbClients; ++c)
            sum += supply[c][j];
        sum -= (nbClients - 1) * opened[j];
        cout << "Adding lazy capacity constraint " << sum << " <= 0" << endl;
        context.rejectCandidate(sum <= 0.0);
        sum.end();
    }*/
}

// Destructor
GenericCallback::~GenericCallback(){}