#ifndef __PhysicalLink__h
#define __PhysicalLink__h
#include <vector>
#include "Slice.h"
#include "Demand.h"

/************************************************************************************/
/*  A PhysicalLink corresponds to an edge in the topology graph.					*/
/*	We consider that the id of each PhysicalLink is in the range [0, ..., m-1]		*/
/*	Each PhysicalLink has a certain number of slices that defines how its 			*/
/* 	spectrum is divided.															*/												
/************************************************************************************/
class PhysicalLink{
private:
	int id;							// Edge identifier. Starts at 0 and goes until n-1.
	int idSource;					// Refers to the source node identifier. 
	int idTarget;					// Refers to the target node identifier. 
	int nbSlices;					// Number of slices in the link
	double length;					// Physical length of the link 
	double cost;					// Cost of routing through this link
	std::vector<Slice> spectrum;	// Link's spectrum


public:
	/************************************************/
	/*				Constructor						*/
	/************************************************/
	PhysicalLink(int i, int s, int t, double l = 0.0, int nb = 1, double c = 0.0);

	/************************************************/
	/*					Getters						*/
	/************************************************/
	int getId() const { return id; }
	int getSource() const { return idSource; }
	int getTarget() const { return idTarget; }
	int getNbSlices() const { return nbSlices; }
	double getLength() const { return length; }
	double getCost() const { return cost; }
	std::vector<Slice> getSlices() const { return spectrum; }
	Slice getSlice_i(int i) const {return spectrum[i];}

	/************************************************/
	/*					Setters						*/
	/************************************************/
	void setId(int i) { this->id = i; }
	void setSource(int s) { this->idSource = s; }
	void setTarget(int t) { this->idTarget = t; }
	void setNbSlices(int nb) { this->nbSlices = nb; }
	void setLength(double l) { this->length = l; }
	void setCost(double c) { this->cost = c; }


	/************************************************/
	/*					Methods						*/
	/************************************************/
	void copyPhysicalLink(PhysicalLink &edge);

	// Assigns a demand d to position p in the spectrum
	void assignSlices(const Demand &d, int p);

	// Displays the link properties
	void displayPhysicalLink();

	// Display the spectrum properties
	void displaySlices();
};

#endif