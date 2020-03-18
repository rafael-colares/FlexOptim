#ifndef __Instance__h
#define __Instance__h

#include <vector>
#include <string>
#include <iostream>
#include "PhysicalLink.h"
#include "Demand.h"
#include "CSVReader.h"
#include "input.h"


/************************************************************************************/
/*  An Instance corresponds to the initial mapping for 								*/
/*	the online procedure for RSA.													*/
/*	This consists of a topology graph where 										*/
/*	some slices of some edges are already occupied.									*/												
/************************************************************************************/
class Instance {
private:
	Input input;
	int nbNodes;						// number of nodes
	std::vector<PhysicalLink> tabEdge;	// vector of edges
	std::vector<Demand> tabDemand;		// vector of demands

public:

	/************************************************/
	/*				Constructor						*/
	/************************************************/
	Instance(const Input &i);
	Instance(const Instance & i);

	/************************************************/
	/*					Getters						*/
	/************************************************/
	int getNbDemands() const { return (int)this->tabDemand.size(); }				// return nb of demands
	int getNbRoutedDemands() const;													// return nb of demands already routed
	int getNbNonRoutedDemands() const {return getNbDemands() - getNbRoutedDemands(); }// return nb of demands not yet routed
	int getNbEdges() const { return (int)this->tabEdge.size(); }					// return nb of edges
	int getNbNodes() const { return this->nbNodes; }								// return nb of nodes
	Input getInput() const { return this->input; }									// return input

	// return the PhysicalLink with given id
	PhysicalLink getPhysicalLinkFromId(int id) const { return this->tabEdge[id]; }	
	// return the vector of PhysicalLinks
	std::vector<PhysicalLink> getTabEdge() const { return this->tabEdge; }

	// return the demand with given id
	Demand getDemandFromId(int id) const { return this->tabDemand[id]; }
	// return the vector of demands
	std::vector<Demand> getTabDemand() const { return this->tabDemand; }
	std::vector<Demand> getNextDemands() const;
	
	/************************************************/
	/*					Setters						*/
	/************************************************/
	void setNbDemands(int nb) { this->tabDemand.resize(nb); }	// set nb of demands
	void setNbNodes(int nb) { this->nbNodes = nb; }				// set nb of nodes
	void setTabEdge(std::vector<PhysicalLink> tab) { this->tabEdge = tab; }
	void setTabDemand(std::vector<Demand> tab) { this->tabDemand = tab; }
	void setEdgeFromId(int id, PhysicalLink &edge);				// set the edge with given id to a given edge
	void setDemandFromId(int id, Demand &demand);				// set the demand with given id to a given demand


	/************************************************/
	/*					Display						*/
	/************************************************/

	void displayInstance();
	void displayTopology();
	void displaySlices();
	void displayRoutedDemands();
	void displayDetailedTopology();
	void displayNonRoutedDemands();

	/************************************************/
	/*					Methods						*/
	/************************************************/
	void createInitialMapping();
	void readTopology(std::string file);
	void readDemands(std::string file);
	void readDemandAssignment(std::string file);


	void generateRandomDemandsFromFile();
	void generateRandomDemands(const int N);
	bool isRoutable(const int i, const int s, const Demand &demand);
	void assignSlicesOfLink(int linkLabel, int slice, const Demand &demand);
};

#endif