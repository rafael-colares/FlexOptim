#include "ExtendedGraph.h"
/** @todo Erase this class **/
ExtendedGraph::ExtendedGraph(const Instance &instance) : edgeMap(g), edgeLabel(g), edgeSlice(g), nodeId(g), nodeLabel(g), nodeSlice(g){
    // For each PhysicalLink (source, target) and each non-used Slice s
    
	std::cout << "--- constructing INITIAL MAPPING GRAPH... --- " << std::endl;
    for (int i = 0; i < instance.getNbEdges(); i++){
        int source = instance.getPhysicalLinkFromId(i).getSource();
        int target = instance.getPhysicalLinkFromId(i).getTarget();
        for (int s = 0; s < instance.getPhysicalLinkFromId(i).getNbSlices(); s++){
            if (instance.getPhysicalLinkFromId(i).getSlice_i(s).isUsed() ==  false){
                // Create nodes (source, s) and (target, s) if they still do not exists.
                
                ListGraph::Node source_s = getNode(source, s);
                ListGraph::Node target_s = getNode(target, s);

                if (source_s == INVALID){
                    source_s = g.addNode();
                    nodeId[source_s] = g.id(source_s);
                    nodeLabel[source_s] = source;
                    nodeSlice[source_s] = s;
                    displayNode(source_s);
                }
                
                if (target_s == INVALID){
                    target_s = g.addNode();
                    nodeId[target_s] = g.id(target_s);
                    nodeLabel[target_s] = target;
                    nodeSlice[target_s] = s;
                    displayNode(target_s);
                }

                
                ListGraph::Edge e = g.addEdge(source_s, target_s);
                edgeMap[e] = g.id(e);
                edgeLabel[e] = i;
                edgeSlice[e] = s;
                displayEdge(e);
            }
        }
    }
   
}

std::string ExtendedGraph::getEdgeName(ListGraph::Edge e){
    std::string u_label = "(" + std::to_string(nodeLabel[g.u(e)]) + ", " + std::to_string(nodeSlice[g.u(e)]) + ")";
    std::string v_label = "(" + std::to_string(nodeLabel[g.v(e)]) + ", " + std::to_string(nodeSlice[g.v(e)]) + ")";
    return u_label + "-" +  v_label;
}
std::string ExtendedGraph::getNodeName(ListGraph::Node v){
    return "(" + std::to_string(nodeLabel[v]) + ", " + std::to_string(nodeSlice[v]) + ")";
}
bool ExtendedGraph::checkNodeExistence(int label, int slice){
    for (ListGraph::NodeIt n(g); n != INVALID; ++n){
        if(nodeLabel[n] == label && nodeSlice[n] == slice){
            return true;
        }
    }
    return false;
}

ListGraph::Node ExtendedGraph::getNode(int label, int slice){
    for (ListGraph::NodeIt n(g); n != INVALID; ++n){
        if(nodeLabel[n] == label && nodeSlice[n] == slice){
            return n;
        }
    }
    return INVALID;
}

void ExtendedGraph::displayNode (ListGraph::Node n){
    std::cout << "Node id #" << nodeId[n] << std::endl;
}
void ExtendedGraph::displayEdge (ListGraph::Edge e){
    std::cout << "Edge (" << nodeLabel[g.u(e)] << ", " << nodeSlice[g.u(e)] 
            << ")--(" << nodeLabel[g.v(e)] << ", " << nodeSlice[g.v(e)] << ") id #" << edgeMap[e] << ". " << std::endl;
}