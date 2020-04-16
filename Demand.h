#ifndef __Demand__h
#define __Demand__h

#include <string>
/******************************************************************
 * This class identifies a demand. A demand is defined by its id, 
 * its source and target node, its load (i.e., how many slices it 
 * requests) and its maximal length (i.e., how long can be its routing
 * path). It also can be already routed or not.
 *****************************************************************/
class Demand
{
private:
	int id;				/**< The demand's id. **/
	int source;			/**< The demand's source node id. **/
	int target;			/**< The demand's target node id. **/
	int load;			/**< Refers to how many slices the demand requires. **/
	double maxLength;	/**< Refers to the maximum length of the path on which the demand can be routed. **/
	bool routed;		/**< True if the demand is already routed. **/
	int sliceAllocation;/**< The id of the last slice assigned to this demand. @warning Equals -1 if not routed yet. **/

public:
	/************************************************/
	/*					Constructors				*/
	/************************************************/
	/** Constructor. **/
	Demand(int id = -1, int s = -1, int t = -1, int l = -1, double max = 0, bool a=false, int slice=-1);

	/************************************************/
	/*					Getters						*/
	/************************************************/
	/** Returns the demand's id. **/
	int getId() const { return id; }
	
	/** Returns the demand's source node id. **/
	int getSource() const { return source; }

	/** Returns the demand's target node id. **/
	int getTarget() const { return target; }
	
	/** Returns the number of slices requested by the demand. **/
	int getLoad() const { return load; }
	
	/** Returns the id of the last slice assigned to the demand. **/
	int getSliceAllocation() const { return sliceAllocation; }
	
	/** Returns the maximum length of the path on which the demand can be routed. **/
	double getMaxLength() const { return maxLength; }

	/** Returns true if the demand has already been routed. **/
	bool isRouted() const { return routed; }

	/** Returns a compact description of the demand in the form (source, target, load). **/
	std::string getString() const;

	/************************************************/
	/*					Setters						*/
	/************************************************/
	/** Changes the demand's id. **/
	void setId(int i) { this->id = i; }
	
	/** Changes the demand's source node id. **/
	void setSource(int s) { this->source = s; }
	
	/** Changes the demand's target node id. **/
	void setTarget(int t) { this->target = t; }
	
	/** Changes the number of slices requested by the demand. **/
	void setLoad(int l) { this->load = l; }

	/** Changes the id of the last slice assigned to the demand. **/
	void setSliceAllocation(int i) { this->sliceAllocation = i; }
	
	/** Changes the maximum length of the path on which the demand can be routed. **/
	void setMaxLength(double max) { this->maxLength = max; }

	/** Changes the status of a demand. **/
	void setRouted(bool a) { this->routed = a; }

	/************************************************/
	/*					Methods						*/
	/************************************************/
	/** Copies all information from a given demand. **/
	void copyDemand(Demand &demand);

	/** Displays demand information. **/
	void displayDemand();

	/** Verifies if the demand has exactly the given informations. **/
	void checkDemand(int id, int source, int target, int load);
};

#endif