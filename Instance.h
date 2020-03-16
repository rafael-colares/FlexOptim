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
	std::vector<Demand> tabOnlineDemand;// vector of new demands

	std::string linkFile = "Parameters/Link.csv";
	std::string demandFile = "Parameters/Demand.csv";
	std::string assignmentFile = "Parameters/Demand_edges_slices.csv";


public:

	/************************************************/
	/*				Constructor						*/
	/************************************************/
	Instance(const Input &i);

	/************************************************/
	/*					Getters						*/
	/************************************************/
	int getNbDemands() const { return (int)this->tabDemand.size(); }				// return nb of demands
	int getNbOnlineDemands() const { return (int)this->tabOnlineDemand.size(); }	// return nb of demands
	int getNbEdges() const { return (int)this->tabEdge.size(); }					// return nb of edges
	int getNbNodes() const { return this->nbNodes; }								// return nb of nodes
	Input getInput() const { return this->input; }									// return input

	// return the PhysicalLink with given id
	PhysicalLink getPhysicalLinkFromId(int id) const { return this->tabEdge[id]; }	
	// return the vector of PhysicalLinks
	std::vector<PhysicalLink> getPhysicalLinks() const { return this->tabEdge; }

	// return the demand with given id
	Demand getDemandFromId(int id) const { return this->tabDemand[id]; }
	// return the vector of demands
	std::vector<Demand> getDemands() const { return this->tabDemand; }

	// return the online demand with given id
	Demand getOnlineDemandFromId(int id) const { return this->tabOnlineDemand[id]; }
	// return the vector of demands
	std::vector<Demand> getOnlineDemands() const { return this->tabOnlineDemand; }	

	/************************************************/
	/*					Setters						*/
	/************************************************/
	void setNbDemands(int nb) { this->tabDemand.resize(nb); }	// set nb of demands
	void setNbNodes(int nb) { this->nbNodes = nb; }				// set nb of nodes
	void setEdgeFromId(int id, PhysicalLink &edge);				// set the edge with given id to a given edge
	void setDemandFromId(int id, Demand &demand);				// set the demand with given id to a given demand


	/************************************************/
	/*					Methods						*/
	/************************************************/
	void createInitialMapping();
	void readTopology(std::string file);
	void readDemands(std::string file);
	void readDemandAssignment(std::string file);
	void displayInstance();
	void displayTopology();
	void displaySlices();
	void displayRoutedDemands();
	void generateRandomDemandsFromFile();
	void generateRandomDemands(const int N);
	bool isRoutable(const int i, const int s, const Demand &demand);
	void assignSlicesOfLink(int linkLabel, int slice, const Demand &demand);
	void displayDetailedTopology();
};

#endif