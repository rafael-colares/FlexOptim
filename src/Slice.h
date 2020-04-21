#ifndef __Slice__h
#define __Slice__h

/*********************************************************************************************
* This class identifies a Slice of a PhysicalLink.								
* A Demand may be routed through a Slice. This class allows to identify, specify and verify 
* which Demand is routed (or not) through the slice.
*********************************************************************************************/
class Slice
{
private:

	int assignedDemand;	/**< Refers to the id of the routed demand; If no demand is routed, it equals -1. **/

public:
	/** Default constructor. **/
	Slice(int id = -1);

	/** Returns the id of the Demand assigned to the Slice. **/
	int getAssignment() const { return assignedDemand; }

	/** Specifies the id of the Demand assigned to the Slice. **/
	void setAssignment(int d) { this->assignedDemand = d; }

	/** Verifies if the Slice is occupied by some demand. **/
	bool isUsed();
	
};

#endif