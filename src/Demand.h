#ifndef __Demand__h
#define __Demand__h

#include <string>
/********************************************************************************************
 * This class identifies a demand. A demand is defined by its id, its source and target node, 
 * its load (i.e., how many slices it requests) and its maximal length (i.e., how long can be
 * its routing path). It also can be already routed or not. \note The id of each Demand is 
 * considered to be in the range [0, ..., n-1].
 ********************************************************************************************/
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
	/****************************************************************************************/
	/*									Constructors										*/
	/****************************************************************************************/
	/** Constructor. @param i Id. @param s Source node id. @param t Target node id. @param l Load. @param maxL Maximum length. @param route Whether the demand is already routed. @param pos The last slice position assigned to the demand. **/
	Demand(int i = -1, int s = -1, int t = -1, int l = 0, double maxL = 0, bool route=false, int pos=-1);

	/****************************************************************************************/
	/*										Getters											*/
	/****************************************************************************************/
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

	/****************************************************************************************/
	/*										Setters											*/
	/****************************************************************************************/
	/** Changes the demand's id. @param i New id value. **/
	void setId(int i) { this->id = i; }
	
	/** Changes the demand's source node id. @param s New source node id. **/
	void setSource(int s) { this->source = s; }
	
	/** Changes the demand's target node id. @param s New target node id. **/
	void setTarget(int t) { this->target = t; }
	
	/** Changes the number of slices requested by the demand. @param l New load value.**/
	void setLoad(int l) { this->load = l; }

	/** Changes the id of the last slice assigned to the demand. @param i The new last slice position.**/
	void setSliceAllocation(int i) { this->sliceAllocation = i; }
	
	/** Changes the maximum length of the path on which the demand can be routed. @param max New demand's maximum length value.**/
	void setMaxLength(double max) { this->maxLength = max; }

	/** Changes the status of a demand. @param b If true, demand becomes routed. Otherwise, it becomes non-routed.**/
	void setRouted(bool b) { this->routed = b; }

	/****************************************************************************************/
	/*										Methods											*/
	/****************************************************************************************/
	/** Copies all information from a given demand. @param demand The demand to be copied. **/
	void copyDemand(Demand &demand);

	/** Verifies if the demand has exactly the given informations. @param id The supposed id value. @param source The supposed source node id value. @param target The supposed target node id value. @param load The supposed load value.**/
	void checkDemand(int id, int source, int target, int load);

	/****************************************************************************************/
	/*										Display											*/
	/****************************************************************************************/
	/** Displays demand information. **/
	void displayDemand();
};

#endif