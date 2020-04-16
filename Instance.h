#ifndef __Instance__h
#define __Instance__h

#include <vector>
#include <string>
#include <iostream>
#include "PhysicalLink.h"
#include "Demand.h"
#include "CSVReader.h"
#include "input.h"


/********************************************************************************************
 * This class stores the initial mapping that serves as input for the Online Routing and 
 * Spectrum Allocation problem. This consists of a topology graph where some slices of some 
 * edges are already occupied by some given demands.												
********************************************************************************************/
class Instance {
private:
	Input input;						/**< An instance needs an input. **/
	int nbNodes;						/**< Number of nodes in the physical network. **/
	std::vector<PhysicalLink> tabEdge;	/**< A set of PhysicalLink. **/
	std::vector<Demand> tabDemand;		/**< A set of Demand (already routed or not). **/

public:

	/****************************************************************************************/
	/*										Constructor										*/
	/****************************************************************************************/
	/** Constructor initializes the object with the information of an Input. **/
	Instance(const Input &i);
	/** Copy constructor. **/
	Instance(const Instance & i);

	/****************************************************************************************/
	/*										Getters											*/
	/****************************************************************************************/
	/** Returns the total number of demands. **/
	int getNbDemands() const { return (int)this->tabDemand.size(); }

	/** Returns the number of demands already routed. **/
	int getNbRoutedDemands() const;	

	/** Returns the number of non-routed demands. **/
	int getNbNonRoutedDemands() const {return getNbDemands() - getNbRoutedDemands(); }

	/** Returns the number of links in the physical network. **/
	int getNbEdges() const { return (int)this->tabEdge.size(); }

	/** Returns the number of nodes in the physical network. **/
	int getNbNodes() const { return this->nbNodes; }

	/** Returns the instance's input. **/
	Input getInput() const { return this->input; }

	/** Returns the PhysicalLink with given index. **/
	PhysicalLink getPhysicalLinkFromId(int index) const { return this->tabEdge[index]; }	

	/** Returns the first PhysicalLink with source s and target t. Should only be called if method hasLink returns true. @warning If there is no such link, the program is aborted! **/
	PhysicalLink getPhysicalLinkBetween(int s, int t);	

	/** Returns the vector of PhysicalLink. **/
	std::vector<PhysicalLink> getTabEdge() const { return this->tabEdge; }

	/** Returns the demand with given index. **/
	Demand getDemandFromIndex(int id) const { return this->tabDemand[id]; }

	/** Returns the vector of Demand. **/
	std::vector<Demand> getTabDemand() const { return this->tabDemand; }
	
	/** Returns the vector of demands to be routed in the next optimization. **/
	std::vector<Demand> getNextDemands() const;
	
	/****************************************************************************************/
	/*										Setters											*/
	/****************************************************************************************/
	/** Change the total number of demands. **/
	void setNbDemands(int nb) { this->tabDemand.resize(nb); }

	/** Change the number of nodes in the physical network. **/
	void setNbNodes(int nb) { this->nbNodes = nb; }	

	/** Change the set of links in the physical network. **/
	void setTabEdge(std::vector<PhysicalLink> tab) { this->tabEdge = tab; }

	/** Change the set of demands. **/
	void setTabDemand(std::vector<Demand> tab) { this->tabDemand = tab; }

	/** Changes the attributes of the PhysicalLink from the given index according to the attributes of the given link. **/
	void setEdgeFromId(int index, PhysicalLink &edge);

	/** Changes the attributes of the Demand from the given index according to the attributes of the given demand. **/
	void setDemandFromId(int id, Demand &demand);

	/****************************************************************************************/
	/*										Methods											*/
	/****************************************************************************************/
	/** Builds the initial mapping based on the information retrived from the Input. **/
	void createInitialMapping();
	
	/** Reads the topology information from file. @warning File should be structured as in Link.csv. **/
	void readTopology(std::string file);
	
	/** Reads the routed demand information from file. @warning File should be structured as in Demand.csv. **/
	void readDemands(std::string file);

	/** Reads the assignment information from file. Sets the demands to routed and update the slices of the edges. @warning File should be structured as in Demand_edges_slices.csv. **/
	void readDemandAssignment(std::string file);

	/** Adds non-routed demands to the pool by reading the information from onlineDemands Input file. **/
	void generateRandomDemandsFromFile();

	/** Adds non-routed demands to the pool by generating N random demands. **/
	void generateRandomDemands(const int N);

	/** Assigns the given demand to the j-th slice of the i-th link. **/
	void assignSlicesOfLink(int i, int j, const Demand &demand);

	/** Verifies if there is enough place for a given demand to be routed through link i on last slice position s. **/
	bool hasEnoughSpace(const int i, const int s, const Demand &demand);

	/** Verifies if there exists a link between nodes of id u and v. **/
	bool hasLink(int u, int v);

	/** Call the methods allowing the build of output files. **/
	void output(std::string i = "");
	
	/** Builds file Edge_Slice_Holes_i.csv containing information about the mapping after n optimizations. **/
	void outputEdgeSliceHols(std::string n);

	/** Builds file Demand.csv containing information about the routed demands. **/
	void outputDemand();
	
	/** Builds file Demand_edges_slices.csv containing information about the assignment of routed demands. **/
	void outputDemandEdgeSlices();

	/****************************************************************************************/
	/*										Display											*/
	/****************************************************************************************/
	/** Displays overall information about the current instance. **/
	void displayInstance();

	/** Displays information about the physical topology. **/	
	void displayTopology();
	
	/** Displays detailed information about state of the physical topology. **/
	void displayDetailedTopology();
	
	/** Displays summarized information about slice occupation. **/
	void displaySlices();
	
	/** Displays information about the routed demands. **/
	void displayRoutedDemands();

	/** Displays information about the non-routed demands. **/
	void displayNonRoutedDemands();

};

#endif