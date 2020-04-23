#ifndef __PhysicalLink__h
#define __PhysicalLink__h
#include <vector>
#include "Slice.h"
#include "Demand.h"

/************************************************************************************
 * This class identifies a link on the physical network. It corresponds to an 
 * edge in the topology graph and is defined by its id, a source and a target node, 
 * its length in the physical network, and a cost. A edge also has a frequency 
 * spectrum that is splitted into a given number of slices. \note The id of each 
 * PhysicalLink is considered to be in the range [0, ..., n-1].												
 ************************************************************************************/
class PhysicalLink{
private:
	int id;							/**< Link identifier. **/
	int idSource;					/**< Source node identifier. **/ 
	int idTarget;					/**< Target node identifier. **/ 
	int nbSlices;					/**< Number of slices the frequency spectrum is divided into. **/
	double length;					/**< Length of the link in the physical network. **/
	double cost;					/**< Cost of routing a demand through the link. **/
	std::vector<Slice> spectrum;	/**< Link's spectrum **/

public:
	/****************************************************************************************/
	/*										Constructor										*/
	/****************************************************************************************/

	/** Constructor. @param i Id. @param s Source node's id. @param t Target node's id. @param len Link's length n the physical network. @param nb Number of slices its frequency spectrum is divided into. @param c Cost of routing a demand through the link. **/
	PhysicalLink(int i, int s, int t, double len = 0.0, int nb = 1, double c = 0.0);

	/****************************************************************************************/
	/*										Getters											*/
	/****************************************************************************************/

	/** Returns the link's id. **/
	int getId() const { return id; }
	
	/** Returns the link's source node id. **/
	int getSource() const { return idSource; }
	
	/** Returns the link's target node id. **/
	int getTarget() const { return idTarget; }
	
	/** Returns the number of slices its frequency spectrum is divided into. **/
	int getNbSlices() const { return nbSlices; }
	
	/** Returns the link's physical length. **/
	double getLength() const { return length; }

	/** Returns the cost of routing a demand through the link. **/
	double getCost() const { return cost; }
	
	/** Returns the frequency spectrum as a vector of Slices. **/
	std::vector<Slice> getSlices() const { return spectrum; }
	
	/** Returns a slice of the frequency spectrum. @param i The slice in the i-th position. **/
	Slice getSlice_i(int i) const {return spectrum[i];}
	
	/** Returns a string summarizing the link's from/to information. **/
	std::string getString() const { return "[" + std::to_string(getSource()+1) + "," + std::to_string(getTarget()+1) + "]"; }

	/****************************************************************************************/
	/*										Setters											*/
	/****************************************************************************************/

	/** Changes the link's id. @param i New id. **/
	void setId(int i) { this->id = i; }
	
	/** Changes the link's source node id. @param s New source node's id. **/
	void setSource(int s) { this->idSource = s; }

	/** Changes the link's target node id. @param t New target node's id. **/
	void setTarget(int t) { this->idTarget = t; }

	/** Changes the number of slices its frequency spectrum is divided into. @param nb New number of slices.**/
	void setNbSlices(int nb) { this->nbSlices = nb; }

	/** Changes the link's length. @param l New length. **/
	void setLength(double l) { this->length = l; }
	
	/** Changes the cost of routing a demand through the link. @param c New cost. **/
	void setCost(double c) { this->cost = c; }

	/****************************************************************************************/
	/*										Methods											*/
	/****************************************************************************************/

	/** Copies all information from a given link. @param link The PhysicalLink to be copied. **/
	void copyPhysicalLink(PhysicalLink &link);

	/** Verifies if the current PhysicalLink routes the a demand. @param dem The demand to be checked. **/
	bool contains(const Demand &dem) const;

	/** Assigns a demand to a given position in the spectrum. @param dem The demand to be assigned. @param p The position of the last slice to be assigned. **/
	void assignSlices(const Demand &dem, int p);

	/** Returns the maximal slice position used in the frequency spectrum. **/
	int getMaxUsedSlicePosition();

	/****************************************************************************************/
	/*										Display											*/
	/****************************************************************************************/

	/** Displays summarized information about the link. **/
	void displayPhysicalLink();

	/** Displays detailed information about state of the link. **/
	void displayDetailedPhysicalLink();

	/** Displays summarized information about slice occupation. **/
	void displaySlices();
};

#endif