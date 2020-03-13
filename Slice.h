#ifndef __Slice__h
#define __Slice__h

/********************************************************************************************/
/*  A Slice may have a demand assigned to it. 												*/
/*	If this is the case, int assignedDemand refers to the demandId of the assigned demand.	*/
/*	If not int assignedDemand equals -1. 													*/
/********************************************************************************************/
class Slice
{
private:
	// each slice has an assigned demand; -1 if no demand is asigned
	int assignedDemand;		
public:
	// Constructor
	Slice(int d = -1);

	// Getters
	int getAssignment() const { return assignedDemand; }	// get the assigned demand
	// Setters
	void setAssignment(int d) { this->assignedDemand = d; }	// set the assigned demand

	// Methods
	bool isUsed();			// returns true if slice is occupied; false otherwise
	
};

#endif