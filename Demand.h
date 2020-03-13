#ifndef __Demand__h
#define __Demand__h
class Demand
{
private:
	int id;
	int source;		// refers to the node source of a demand
	int target;		// refers to the target node of a demand
	int load;		// refers to how many slices the demand requires
	double maxLength;	// refers to the maximum length of the path on which it can be routed
	bool routed;

public:
	// Constructors
	Demand(int id = 0, int s = -1, int t = -1, int l = -1, double max = 0, bool a=false);

	/************************************************/
	/*					Getters						*/
	/************************************************/
	int getId() const { return id; }
	int getSource() const { return source; }
	int getTarget() const { return target; }
	int getLoad() const { return load; }
	double getMaxLength() const { return maxLength; }
	bool isRouted() const { return routed; }

	/************************************************/
	/*					Setters						*/
	/************************************************/
	void setId(int i) { this->id = i; }
	void setSource(int s) { this->source = s; }
	void setTarget(int t) { this->target = t; }
	void setLoad(int l) { this->load = l; }
	void setMaxLength(double max) { this->maxLength = max; }
	void setRouted(bool a) { this->routed = a; }

	/************************************************/
	/*					Methods						*/
	/************************************************/
	void copyDemand(Demand &demand);
	void displayDemand();
	void checkDemand(int i, int s, int t, int l);
};

#endif