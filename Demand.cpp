#include "Demand.h"
#include <iostream>

Demand::Demand(int i, int s, int t, int l, double max, bool a){
	this->setId(i);
	this->setSource(s);
	this->setTarget(t);
	this->setLoad(l);
	this->setMaxLength(max);
	this->setRouted(a);
}

void Demand::copyDemand(Demand & demand){
	this->setId(demand.getId());
	this->setSource(demand.getSource());
	this->setTarget(demand.getTarget());
	this->setLoad(demand.getLoad());
	this->setMaxLength(demand.getMaxLength());
	this->setRouted(demand.isRouted());
}

void Demand::displayDemand(){
	std::cout << "#" << this->getId()+1 << ". " << this->getSource()+1 << " -- " << this->getTarget()+1;
	std::cout << ". nbSlices: " << this->getLoad() << ", maxLength: " << this->getMaxLength() << std::endl;
}

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
