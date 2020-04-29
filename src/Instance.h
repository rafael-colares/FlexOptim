#ifndef __Instance__h
#define __Instance__h

#include <vector>
#include <string>
#include <iostream>
#include <float.h>

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
	int nbInitialDemands;				/**< The number of demands routed in the first initial mapping. **/

public:

	/****************************************************************************************/
	/*										Constructor										*/
	/****************************************************************************************/
	/** Constructor initializes the object with the information of an Input. @param i The input used for creating the instance.**/
	Instance(const Input &i);

	/** Copy constructor. @param i The instance to be copied. **/
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

	/** Returns the number of demands routed in the first initial mapping. **/
	int getNbInitialDemands() const { return nbInitialDemands; }

	/** Returns the number of links in the physical network. **/
	int getNbEdges() const { return (int)this->tabEdge.size(); }

	/** Returns the number of nodes in the physical network. **/
	int getNbNodes() const { return this->nbNodes; }

	/** Returns the instance's input. **/
	Input getInput() const { return this->input; }

	/** Returns the PhysicalLink with given index. @param index The index of PhysicalLink required in tabEdge. **/
	PhysicalLink getPhysicalLinkFromId(int index) const { return this->tabEdge[index]; }	

	/** Returns the first PhysicalLink with the given source and target.  @warning Should only be called if method hasLink returns true. If there is no such link, the program is aborted! @param s Source node id. @param t Target node id. **/
	PhysicalLink getPhysicalLinkBetween(int s, int t);	

	/** Returns the vector of PhysicalLink. **/
	std::vector<PhysicalLink> getTabEdge() const { return this->tabEdge; }

	/** Returns the demand with given index. @param index The index of Demand required in tabDemand.**/
	Demand getDemandFromIndex(int index) const { return this->tabDemand[index]; }

	/** Returns the vector of Demand. **/
	std::vector<Demand> getTabDemand() const { return this->tabDemand; }
	
	/** Returns the vector of demands to be routed in the next optimization. **/
	std::vector<Demand> getNextDemands() const;
	
	/****************************************************************************************/
	/*										Setters											*/
	/****************************************************************************************/

	/** Change the total number of demands. @param nb New total number of demands. @warning This function resizes vector tabDemand, creating default demands if nb is greater than the previous size.**/
	void setNbDemands(int nb) { this->tabDemand.resize(nb); }

	/** Returns the number of demands routed in the first initial mapping. @param nb New number of demands.**/
	void setNbInitialDemands(int nb) { this->nbInitialDemands = nb; }

	/** Change the number of nodes in the physical network. @param nb New number of nodes. **/
	void setNbNodes(int nb) { this->nbNodes = nb; }	

	/** Change the set of links in the physical network. @param tab New vector of PhysicalLinks. **/
	void setTabEdge(std::vector<PhysicalLink> tab) { this->tabEdge = tab; }

	/** Change the set of demands. @param tab New vector of Demands. **/
	void setTabDemand(std::vector<Demand> tab) { this->tabDemand = tab; }

	/** Changes the attributes of the PhysicalLink from the given index according to the attributes of the given link. @param i The index of the PhysicalLink to be changed. @param link the PhysicalLink to be copied. **/
	void setEdgeFromId(int i, PhysicalLink &link);

	/** Changes the attributes of the Demand from the given index according to the attributes of the given demand. @param i The index of the Demand to be changed. @param demand the Demand to be copied. **/
	void setDemandFromId(int i, Demand &demand);

	/****************************************************************************************/
	/*										Methods											*/
	/****************************************************************************************/
	
	/** Builds the initial mapping based on the information retrived from the Input. **/
	void createInitialMapping();
	
	/** Reads the topology information from input's linkFile. Builds the set of links. @warning File should be structured as in Link.csv. **/
	void readTopology();
	
	/** Reads the routed demand information from input's demandFile. Builds the set of demands. @warning File should be structured as in Demand.csv. **/
	void readDemands();

	/** Reads the assignment information from input's assignmentFile. Sets the demands to routed and update the slices of the edges. @warning File should be structured as in Demand_edges_slices.csv. **/
	void readDemandAssignment();

	/** Adds non-routed demands to the pool by reading the information from file. @param filePath The path of the file to be read. **/
	void generateRandomDemandsFromFile(std::string filePath);

	/** Adds non-routed demands to the pool by generating random demands. @param N The number of random demands to be generated. **/
	void generateRandomDemands(const int N);

	/** Assigns a demand to a slice of a link. @param index The index of the PhysicalLink to be modified. @param pos The last slice position. @param demand The demand to be assigned. **/
	void assignSlicesOfLink(int index, int pos, const Demand &demand);

	/** Verifies if there is enough place for a given demand to be routed through a link on a last slice position. @param index The index of the PhysicalLink to be inspected. @param pos The last slice positon. @param demand The candidate demand to be assigned. **/
	bool hasEnoughSpace(const int index, const int pos, const Demand &demand);

	/** Verifies if there exists a link between two nodes. @param u Source node id. @param v Target node id. **/
	bool hasLink(int u, int v);

	/** Call the methods allowing the build of output files. @param i Indicates the current iteration. **/
	void output(std::string i = "0");
	
	/** Builds file Edge_Slice_Holes_i.csv containing information about the mapping after n optimizations. @param i The i-th output file to be generated. **/
	void outputEdgeSliceHols(std::string i);

	/** Builds file Demand.csv containing information about the routed demands. **/
	void outputDemand();
	
	/** Builds file Demand_edges_slices.csv containing information about the assignment of routed demands. @param i The i-th output file to be generated. **/
	void outputDemandEdgeSlices(std::string i);

	/** Builds file results.csv containing information about the main obtained results. @param fileName The name of demand file being optimized. **/
	void outputLogResults(std::string fileName);
	
	/****************************************************************************************/
	/*										Display											*/
	/****************************************************************************************/
	/** Displays overall information about the current instance. **/
	void displayInstance();

	/** Displays information about the physical topology. **/	
	void displayTopology();
	
	/** Displays detailed information about state of the physical topology. **/
	void displayDetailedTopology();
	
	/** Displays summarized information about slice occupation of each PhysicalLink. **/
	void displaySlices();
	
	/** Displays information about the routed demands. **/
	void displayRoutedDemands();

	/** Displays information about the non-routed demands. **/
	void displayNonRoutedDemands();

};

#endif