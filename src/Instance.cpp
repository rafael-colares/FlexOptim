#include "Instance.h"

/************************************************/
/*				Constructors					*/
/************************************************/
/* Constructor initializes the object with the information of an Input. */
Instance::Instance(const Input &i) : input(i){
	this->setNbNodes(0);
}

/* Copy constructor. */
Instance::Instance(const Instance &i) : input(i.getInput()){
	this->setNbNodes(i.getNbNodes());
	this->setTabEdge(i.getTabEdge());
	this->setTabDemand(i.getTabDemand());
}

/************************************************/
/*					Methods						*/
/************************************************/
/* Returns the number of demands already routed. */
int Instance::getNbRoutedDemands() const{
	int counter = 0;
	for(int i = 0; i < getNbDemands(); i++){
    	if (tabDemand[i].isRouted() == true){
			counter++;
		}
	}
	return counter;
}

/* Returns the vector of demands to be routed in the next optimization. */
std::vector<Demand> Instance::getNextDemands() const { 
	std::vector<Demand> toBeRouted;
	for(int i = 0; i < getNbDemands(); i++){
		if( (tabDemand[i].isRouted() == false) && ((int)toBeRouted.size() < getInput().getNbDemandsAtOnce()) ){
			toBeRouted.push_back(tabDemand[i]);
		}
	}
	return toBeRouted;
}
	

/* Changes the attributes of the PhysicalLink from the given index according to the attributes of the given link. */
void Instance::setEdgeFromId(int id, PhysicalLink & edge){
	this->tabEdge[id].copyPhysicalLink(edge);
}

/* Changes the attributes of the Demand from the given index according to the attributes of the given demand. */
void Instance::setDemandFromId(int id, Demand & demand){
	this->tabDemand[id].copyDemand(demand);
}

/* Builds the initial mapping based on the information retrived from the Input. */
void Instance::createInitialMapping(){
	readTopology();
	readDemands();
	readDemandAssignment();
	setNbInitialDemands(getNbRoutedDemands());
}

/* Reads the topology information from file. */
void Instance::readTopology(){
	std::cout << "Reading " << input.getLinkFile() << " ..."  << std::endl;
	CSVReader reader(input.getLinkFile());
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

/* Reads the routed demand information from file. */
void Instance::readDemands(){
	std::cout << "Reading " << input.getDemandFile() << " ..." << std::endl;
	CSVReader reader(input.getDemandFile());
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

/* Reads the assignment information from file. */
void Instance::readDemandAssignment(){
	CSVReader reader(input.getAssignmentFile());
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
	std::cout << "Checking done." << std::endl;

	//search for slice allocation line
	for (int alloc = 0; alloc < numberOfLines; alloc++)	{
		if (dataList[alloc][0].find("slice allocation") != std::string::npos) {
			// for each demand
			for (int d = 0; d < this->getNbDemands(); d++) {
				int demandMaxSlice = std::stoi(dataList[alloc][d+1]) - 1;
				this->tabDemand[d].setRouted(true);
				this->tabDemand[d].setSliceAllocation(demandMaxSlice);
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

/* Displays overall information about the current instance. */
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

/* Displays information about the physical topology. */
void Instance::displayTopology(){
	std::cout << std::endl << "--- The Physical Topology ---" << std::endl;
	for (int i = 0; i < this->getNbEdges(); i++) {
		tabEdge[i].displayPhysicalLink();
	}
	std::cout << std::endl;
}


/* Displays detailed information about state of the physical topology. */
void Instance::displayDetailedTopology(){
	std::cout << std::endl << "--- The Detailed Physical Topology ---" << std::endl;
	for (int i = 0; i < this->getNbEdges(); i++) {
		tabEdge[i].displayDetailedPhysicalLink();
	}
	std::cout << std::endl;

}

/* Displays summarized information about slice occupation of each PhysicalLink. */
void Instance::displaySlices() {
	std::cout << std::endl << "--- Slice occupation ---" << std::endl;
	for (int i = 0; i < this->getNbEdges(); i++) {
		std::cout << "#" << i+1 << ". ";
		tabEdge[i].displaySlices();
	}
	std::cout << std::endl;
}

/* Displays information about the routed demands. */
void Instance::displayRoutedDemands(){
	std::cout << std::endl << "--- The Routed Demands ---" << std::endl;
	for (int i = 0; i < this->getNbDemands(); i++) {
		if (tabDemand[i].isRouted()) {
			tabDemand[i].displayDemand();
		}
	}
	std::cout << std::endl;
}

/* Adds non-routed demands to the pool by reading the information from onlineDemands Input file. */
void Instance::generateRandomDemandsFromFile(std::string filePath){
	
	std::cout << "Reading " << filePath << " ..." << std::endl;
	CSVReader reader(filePath);
	/* dataList is a vector of vectors of strings. */
	/* dataList[0] corresponds to the first line of the document and dataList[0][i] to the i-th word.*/
	std::vector<std::vector<std::string> > dataList = reader.getData();
	int numberOfLines = (int)dataList.size();
	//skip the first line (headers)
	for (int i = 1; i < numberOfLines; i++) {
		int idDemand = std::stoi(dataList[i][0]) - 1 + getNbRoutedDemands();
		int demandSource = std::stoi(dataList[i][1]) - 1;
		int demandTarget = std::stoi(dataList[i][2]) - 1;
		int demandLoad = std::stoi(dataList[i][3]);
		double DemandMaxLength = std::stod(dataList[i][4]);
		Demand demand(idDemand, demandSource, demandTarget, demandLoad, DemandMaxLength, false);
		this->tabDemand.push_back(demand);
	}
}

/* Adds non-routed demands to the pool by generating N random demands. */
void Instance::generateRandomDemands(const int N){
	srand (1234567890);
	for (int i = 0; i < N; i++){
		int idDemand =  i + getNbRoutedDemands();
		int demandSource = rand() % getNbNodes();
		int demandTarget = rand() % getNbNodes();
		while (demandTarget == demandSource){
			demandTarget = rand() % getNbNodes();
		} 
		int demandLoad = 3;
		double DemandMaxLength = 3000;
		Demand demand(idDemand, demandSource, demandTarget, demandLoad, DemandMaxLength, false);
		this->tabDemand.push_back(demand);
	}
}

/* Verifies if there is enough place for a given demand to be routed through link i on last slice position s. */
bool Instance::hasEnoughSpace(const int i, const int s, const Demand &demand){
	// std::cout << "Calling hasEnoughSpace..." << std::endl;
	const int LOAD = demand.getLoad();
	int firstPosition = s - LOAD + 1;
	if (firstPosition < 0){
		return false;
	}
	for (int pos = firstPosition; pos <= s; pos++){
		if (getPhysicalLinkFromId(i).getSlice_i(pos).isUsed() ==  true){
			return false;
		}
	}
	// std::cout << "Called hasEnoughSpace." << std::endl;
	return true;
}

/* Assigns the given demand to the j-th slice of the i-th link. */
void Instance::assignSlicesOfLink(int linkLabel, int slice, const Demand &demand){
	this->tabEdge[linkLabel].assignSlices(demand, slice);
	this->tabDemand[demand.getId()].setRouted(true);
	this->tabDemand[demand.getId()].setSliceAllocation(slice);
}


/* Displays information about the non-routed demands. */
void Instance::displayNonRoutedDemands(){
	std::cout << std::endl << "--- The Non Routed Demands ---" << std::endl;
	for (int i = 0; i < this->getNbDemands(); i++) {
		if (tabDemand[i].isRouted() == false) {
			tabDemand[i].displayDemand();
		}
	}
	std::cout << std::endl;

}

/* Call the methods allowing the build of output files. */
void Instance::output(std::string i){
	std::cout << "Output " << i << std::endl;
	outputEdgeSliceHols(i);
	//outputDemand();
	outputDemandEdgeSlices(i);
}

/* Builds file Demand_edges_slices.csv containing information about the assignment of routed demands. */
void Instance::outputDemandEdgeSlices(std::string counter){
	std::string delimiter = ";";
	std::string filePath = this->input.getOutputPath() + "Demand_edges_slices_" + counter + ".csv";
	std::ofstream myfile(filePath.c_str(), std::ios::out | std::ios::trunc);
	if (myfile.is_open()){
		myfile << "edge_slice_demand" << delimiter;
		for (int i = 0; i < getNbDemands(); i++){
			if (getDemandFromIndex(i).isRouted()){
				myfile << "k_" << getDemandFromIndex(i).getId()+1 << "= " << getDemandFromIndex(i).getString() << delimiter;
			}
		}
		myfile << "\n";
		for (int e = 0; e < getNbEdges(); e++){
			myfile << getPhysicalLinkFromId(e).getString() << delimiter;
			for (int i = 0; i < getNbDemands(); i++){
				if (getDemandFromIndex(i).isRouted()){
					// if demand is routed through edge: 1
					if (getPhysicalLinkFromId(e).contains(getDemandFromIndex(i)) == true){
						myfile << "1" << delimiter;
					}
					else{
						myfile << " " << delimiter;
					}
				}
			}
			myfile << "\n";
		}
		myfile << " slice allocation " << delimiter;
		for (int i = 0; i < getNbDemands(); i++){
			if (getDemandFromIndex(i).isRouted()){
				myfile << getDemandFromIndex(i).getSliceAllocation()+1 << delimiter;
			}
		}
		myfile << "\n";
	}
}

/* Builds file Demand.csv containing information about the routed demands. */
void Instance::outputDemand(){
	std::string delimiter = ";";
	std::string filePath = this->input.getOutputPath() + "Demand" + ".csv";
	std::ofstream myfile(filePath.c_str(), std::ios::out | std::ios::trunc);
	if (myfile.is_open()){
		myfile << "index" << delimiter;
		myfile << "origin" << delimiter;
		myfile << "destination" << delimiter;
		myfile << "slots" << delimiter;
		myfile << "max_length" << "\n";
		for (int i = 0; i < getNbDemands(); i++){
			if (getDemandFromIndex(i).isRouted()){
				myfile << std::to_string(getDemandFromIndex(i).getId()+1) << delimiter;
				myfile << std::to_string(getDemandFromIndex(i).getSource()+1) << delimiter;
				myfile << std::to_string(getDemandFromIndex(i).getTarget()+1) << delimiter;
				myfile << std::to_string(getDemandFromIndex(i).getLoad()) << delimiter;
				myfile << std::to_string(getDemandFromIndex(i).getMaxLength()) << "\n";
			}
		}
	}
}

/* Builds file Edge_Slice_Holes_i.csv containing information about the mapping after n optimizations. */
void Instance::outputEdgeSliceHols(std::string counter){
	std::cout << "Output EdgeSliceHols: " << counter << std::endl;
	std::string delimiter = ";";
	std::string filePath = this->input.getOutputPath() + "Edge_Slice_Holes_" + counter + ".csv";
	std::ofstream myfile(filePath.c_str(), std::ios::out | std::ios::trunc);
	if (myfile.is_open()){
		myfile << " Slice-Edge " << delimiter;
		for (int i = 0; i < getNbEdges(); i++){
			std::string edge = "e_" + std::to_string(i+1);
			myfile << edge << delimiter;
		}
		myfile << "\n";
		for (int s = 0; s < input.getnbSlicesInOutputFile(); s++){
			std::string slice = "s_" + std::to_string(s+1);
			myfile << slice << delimiter;
			for (int i = 0; i < getNbEdges(); i++){
				if (s < getPhysicalLinkFromId(i).getNbSlices() && getPhysicalLinkFromId(i).getSlice_i(s).isUsed() == true){
					myfile << "1" << delimiter;
				}
				else{
					myfile << "0" << delimiter;
				}
			}
			myfile << "\n";
		}
		myfile << "Nb_New_Demands:" << delimiter << getNbRoutedDemands() - getNbInitialDemands() << "\n";
	}
	else{
		std::cerr << "Unable to open file.\n";
	}
  	myfile.close();
}


/* Builds file results.csv containing information about the main obtained results. */
void Instance::outputLogResults(std::string fileName){
	std::string delimiter = ";";
	std::string filePath = this->input.getOutputPath() + "results.csv";
	std::ofstream myfile(filePath.c_str(), std::ios_base::app);
	if (myfile.is_open()){
		myfile << fileName << delimiter;
		int nbRouted = getNbRoutedDemands();
		myfile << nbRouted - getNbInitialDemands() << delimiter;
		myfile << nbRouted << "\n";
	}
}

/* Verifies if there exists a link between nodes of id u and v. */
bool Instance::hasLink(int u, int v){
	for (unsigned int e = 0; e < tabEdge.size(); e++){
		if ((tabEdge[e].getSource() == u) && (tabEdge[e].getTarget() == v)){
			return true;
		}
		if ((tabEdge[e].getSource() == v) && (tabEdge[e].getTarget() == u)){
			return true;
		}
	}
	return false;
}

/* Returns the first PhysicalLink with source s and target t. */
PhysicalLink Instance::getPhysicalLinkBetween(int u, int v){
	for (unsigned int e = 0; e < tabEdge.size(); e++){
		if ((tabEdge[e].getSource() == u) && (tabEdge[e].getTarget() == v)){
			return tabEdge[e];
		}
		if ((tabEdge[e].getSource() == v) && (tabEdge[e].getTarget() == u)){
			return tabEdge[e];
		}
	}
	std::cerr << "Did not found a link between " << u << " and " << v << "!!\n";
	exit(0);
	PhysicalLink link(-1,-1,-1);
	return link;
}