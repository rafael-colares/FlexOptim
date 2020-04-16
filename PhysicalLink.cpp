#include "PhysicalLink.h"
#include <iostream>

PhysicalLink::PhysicalLink(int i, int s, int t, double l, int nb, double c) {
	this->setId(i);
	this->setSource(s);
	this->setTarget(t);
	this->setLength(l);
	this->setNbSlices(nb);
	for (int i = 0; i < nb; i++){
		this->spectrum.push_back(Slice());
	}
	this->setCost(c);
}

void PhysicalLink::copyPhysicalLink(PhysicalLink & edge){
	this->setSource(edge.getSource());
	this->setTarget(edge.getTarget());
	this->setNbSlices(edge.getNbSlices());
	this->spectrum.resize(this->getNbSlices());
	for (int i = 0; i < this->getNbSlices(); i++){
		this->spectrum[i].setAssignment(edge.spectrum[i].getAssignment());
	}
	this->setLength(edge.getLength());
	this->setCost(edge.getCost());
}

/* Verifies if the current Physical Link routes the given demand d. */
bool PhysicalLink::contains(const Demand &d) const{
	for (int i = 0; i < this->getNbSlices(); i++){
		if (getSlice_i(i).getAssignment() == d.getId()){
			return true;
		}
	}
	return false;
}

/* Assigns a demand d to position p in the spectrum. */
void PhysicalLink::assignSlices(const Demand &d, int p){
	int demandLoad = d.getLoad();
	// assign demand d to this edge from position p - demandLoad + 1 to position p
	int first = p - demandLoad + 1;
	for (int i = first; i <= p; i++) {
		this->spectrum[i].setAssignment(d.getId());
	}
}

void PhysicalLink::displayPhysicalLink(){
	std::cout << "#" << this->getId()+1 << ". " << this->getSource()+1 << " -- " << this->getTarget()+1;
	std::cout << ". lenght: " << this->getLength() << ", cost: " << this->getCost() << std::endl;
}
void PhysicalLink::displayDetailedPhysicalLink(){
	std::cout << "#" << this->getId()+1 << ". " << this->getSource()+1 << " -- " << this->getTarget()+1;
	std::cout << ". lenght: " << this->getLength() << ", cost: " << this->getCost() << std::endl;
	for (int i = 0; i < this->getNbSlices(); i++){
		std::cout << "\tSlice #" << i+1 << ". ";
		if (this->spectrum[i].isUsed()) {
			std::cout << this->spectrum[i].getAssignment()+1 << std::endl;
		}
		else{
			std::cout << "--" << std::endl;
		}
	}
}

void PhysicalLink::displaySlices(){
	for (int i = 0; i < this->getNbSlices(); i++){
		if (this->spectrum[i].isUsed()) {
			std::cout << "*";
		}
		else {
			std::cout << " ";
		}
	}
	std::cout << std::endl;
}
