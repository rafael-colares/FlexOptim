#include "Demand.h"
#include <iostream>

/* Constructor. */
Demand::Demand(int i, int s, int t, int l, double max, bool a, int slice){
	this->setId(i);
	this->setSource(s);
	this->setTarget(t);
	this->setLoad(l);
	this->setMaxLength(max);
	this->setRouted(a);
	this->setSliceAllocation(slice);
}


/* Copies all information from a given demand. */
void Demand::copyDemand(Demand & demand){
	this->setId(demand.getId());
	this->setSource(demand.getSource());
	this->setTarget(demand.getTarget());
	this->setLoad(demand.getLoad());
	this->setMaxLength(demand.getMaxLength());
	this->setRouted(demand.isRouted());
	this->setSliceAllocation(demand.getSliceAllocation());
}

/* Displays demand information. */
void Demand::displayDemand(){
	std::cout << "#" << this->getId()+1 << ". " << this->getSource()+1 << " -- " << this->getTarget()+1;
	std::cout << ". nbSlices: " << this->getLoad() << ", maxLength: " << this->getMaxLength() << std::endl;
}

/* Verifies if the demand has exactly the given informations. */
void Demand::checkDemand(int i, int s, int t, int l){
	try {
		if (this->getId() != i || this->getSource() != s || this->getTarget() != t || this->getLoad() != l) {
			throw "Demands do not match. Verify files Demand.csv and Demand_edges_slices.csv\n";
		}
	}
	catch (const char* msg) {
		std::cerr << msg << std::endl;
	}
}

/* Returns a compact description of the demand in the form (source, target, load). */
std::string Demand::getString() const{
	std::string str = "(" + std::to_string(this->getSource() + 1) + ",";
	str += std::to_string(this->getTarget() + 1) + ",";
	str += std::to_string(this->getLoad()) + ")";
	return str;
}
