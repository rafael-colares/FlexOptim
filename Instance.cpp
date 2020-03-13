#include "Instance.h"

Instance::Instance(const Input &i) : input(i){
	this->setNbNodes(0);
}

// Sets the edge with the given id with the attributes of a given edge. Notice that index in tabEdge = id-1!
void Instance::setEdgeFromId(int id, PhysicalLink & edge){
	this->tabEdge[id].copyPhysicalLink(edge);
}

void Instance::setDemandFromId(int id, Demand & demand){
	this->tabDemand[id].copyDemand(demand);
}

void Instance::createInitialMapping(){
	readTopology(input.getLinkFile());
	readDemands(input.getDemandFile());
	readDemandAssignment(input.getAssignmentFile());
}

// Reads the topology from file Link.cvs
void Instance::readTopology(std::string file){
	std::cout << "Reading " << input.getLinkFile() << " ..."  << std::endl;
	CSVReader reader(file);
	/* dataList is a vector of vectors of strings. */
	/* dataList[0] corresponds to the first line of the document and dataList[0][i] to the i-th word.*/
	std::vector<std::vector<std::string> > dataList = reader.getData();
	int numberOfLines = (int)dataList.size();
	// The number of nodes is given by the max index of sources and targets
	int maxNode = 0;
	//skip the first line (headers)
	// edges and nodes id starts on 1 in the input files. In this program ids will be in the range [0,n-1]!
	for (int i = 1; i < numberOfLines; i++)	{
		int idEdge = std::stoi(dataList[i][0]) - 1;
		int edgeSource = std::stoi(dataList[i][1]) - 1;
		int edgeTarget = std::stoi(dataList[i][2]) - 1;
		double edgeLength = std::stod(dataList[i][3]);
		int edgeNbSlices = std::stoi(dataList[i][4]);
		double edgeCost = std::stod(dataList[i][5]);
		PhysicalLink edge(idEdge, edgeSource, edgeTarget, edgeLength, edgeNbSlices, edgeCost);
		this->tabEdge.push_back(edge);
		if (edgeSource > maxNode) {
			maxNode = edgeSource;
		}
		if (edgeTarget > maxNode) {
			maxNode = edgeTarget;
		}
		std::cout << "Creating edge ";
		edge.displayPhysicalLink();
	}
	this->setNbNodes(maxNode+1);
}

// Reads the demands from file Demand.cvs
void Instance::readDemands(std::string file){
	std::cout << "Reading " << input.getDemandFile() << " ..." << std::endl;
	CSVReader reader(file);
	/* dataList is a vector of vectors of strings. */
	/* dataList[0] corresponds to the first line of the document and dataList[0][i] to the i-th word.*/
	std::vector<std::vector<std::string> > dataList = reader.getData();
	int numberOfLines = (int)dataList.size();
	//skip the first line (headers)
	for (int i = 1; i < numberOfLines; i++) {
		int idDemand = std::stoi(dataList[i][0]) - 1;
		int demandSource = std::stoi(dataList[i][1]) - 1;
		int demandTarget = std::stoi(dataList[i][2]) - 1;
		int demandLoad = std::stoi(dataList[i][3]);
		double DemandMaxLength = std::stod(dataList[i][4]);
		Demand demand(idDemand, demandSource, demandTarget, demandLoad, DemandMaxLength, false);
		this->tabDemand.push_back(demand);
	}
}

//Sets the demands to routed and update the slices of the edges
void Instance::readDemandAssignment(std::string file){
	CSVReader reader(file);
	std::cout << "Reading " << input.getAssignmentFile() << " ..." << std::endl;

	/* dataList is a vector of vectors of strings. */
	/* dataList[0] corresponds to the first line of the document and dataList[0][0] to the first word.*/
	std::vector<std::vector<std::string> > dataList = reader.getData();
	int numberOfColumns = (int)dataList[0].size();
	int numberOfLines = (int)dataList.size();

	//check if the demands in this file are the same as the ones read in Demand.csv
	//skip the first word (headers) and the last one (empty)
	for (int i = 1; i < numberOfColumns-1; i++) {
		int demandId = stoi(getInBetweenString(dataList[0][i], "_", "=")) - 1;
		std::string demandStr = getInBetweenString(dataList[0][i], "(", ")");
		std::vector<std::string> demand = splitBy(demandStr, ",");
		int demandSource = std::stoi(demand[0]) - 1;
		int demandTarget = std::stoi(demand[1]) - 1;
		int demandLoad = std::stoi(demand[2]);
		this->tabDemand[demandId].checkDemand(demandId, demandSource, demandTarget, demandLoad);
	}

	//search for slice allocation line
	for (int alloc = 0; alloc < numberOfLines; alloc++)	{
		if (dataList[alloc][0].find("slice allocation") != std::string::npos) {
			// for each demand
			for (int d = 0; d < this->getNbDemands(); d++) {
				int demandMaxSlice = std::stoi(dataList[alloc][d+1]);
				this->tabDemand[d].setRouted(true);
				// look for which edges the demand is routed
				for (int i = 0; i < this->getNbEdges(); i++) {
					if (dataList[i+1][d+1] == "1") {
						this->tabEdge[i].assignSlices(this->tabDemand[d], demandMaxSlice);
					}
				}
			}
		}
	}
}

void Instance::displayInstance() {
	std::cout << "**********************************" << std::endl;
	std::cout << "*      Constructed Instance      *" << std::endl;
	std::cout << "**********************************" << std::endl;
	std::cout << "Number of nodes : " << this->getNbNodes() << std::endl;
	std::cout << "Number of edges : " << this->getNbEdges() << std::endl;

	displayTopology();
	displaySlices();
	displayRoutedDemands();

}

void Instance::displayTopology(){
	std::cout << std::endl << "--- The Physical Topology ---" << std::endl;
	for (int i = 0; i < this->getNbEdges(); i++) {
		tabEdge[i].displayPhysicalLink();
	}
	std::cout << std::endl;
}

void Instance::displaySlices() {
	std::cout << std::endl << "--- Slice occupation ---" << std::endl;
	for (int i = 0; i < this->getNbEdges(); i++) {
		std::cout << "#" << i+1 << ". ";
		tabEdge[i].displaySlices();
	}
	std::cout << std::endl;
}

void Instance::displayRoutedDemands(){
	std::cout << std::endl << "--- The Routed Demands ---" << std::endl;
	for (int i = 0; i < this->getNbDemands(); i++) {
		if (tabDemand[i].isRouted()) {
			tabDemand[i].displayDemand();
		}
	}
	std::cout << std::endl;
}

void Instance::generateRandomDemands(){
	
	std::cout << "Reading " << input.getOnlineDemandFile() << " ..." << std::endl;
	CSVReader reader(input.getOnlineDemandFile());
	/* dataList is a vector of vectors of strings. */
	/* dataList[0] corresponds to the first line of the document and dataList[0][i] to the i-th word.*/
	std::vector<std::vector<std::string> > dataList = reader.getData();
	int numberOfLines = (int)dataList.size();
	//skip the first line (headers)
	for (int i = 1; i < numberOfLines; i++) {
		int idDemand = std::stoi(dataList[i][0]) - 1 + getNbDemands();
		int demandSource = std::stoi(dataList[i][1]) - 1;
		int demandTarget = std::stoi(dataList[i][2]) - 1;
		int demandLoad = std::stoi(dataList[i][3]);
		double DemandMaxLength = std::stod(dataList[i][4]);
		Demand demand(idDemand, demandSource, demandTarget, demandLoad, DemandMaxLength, false);
		this->tabOnlineDemand.push_back(demand);
	}
}

bool Instance::isRoutable(const int i, const int s, const Demand &demand){
	const int LOAD = demand.getLoad();
	if (s < LOAD - 1){
		return false;
	}
	for (int slice = s - LOAD + 1; slice <= s; slice++){
		if (getPhysicalLinkFromId(i).getSlice_i(slice).isUsed() ==  true){
			return false;
		}
	}
	return true;
}

void Instance::assignSlicesOfLink(int link, const Demand &demand, int slice){
	this->tabEdge[link].assignSlices(demand, slice);
}