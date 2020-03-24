#include <ilcplex/ilocplex.h>
ILOSTLBEGIN
// cplex10Test.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>

#include <lemon/list_graph.h>
#include <lemon/concepts/graph.h>

#include "PhysicalLink.h"
#include "Instance.h"
#include "ExtendedGraph.h"
#include "input.h"
#include "cplexForm.h"

using namespace lemon;

int main(int argc, char *argv[]) {
	try{
		std::string parameterFile;
		if (argc != 2){
			throw std::invalid_argument( "did not receive an argument" );
		}
		else{
			parameterFile = argv[1];
		}
		std::cout << "PARAMETER FILE: " << parameterFile << std::endl;
		Input input(parameterFile);

		std::cout << "--- READING INSTANCE... --- " << std::endl;
		Instance instance(input);

		std::cout << "--- CREATING INITIAL MAPPING... --- " << std::endl;
		instance.createInitialMapping();
		std::cout << instance.getNbRoutedDemands() << " were routed." << std::endl;

		//instance.displayDetailedTopology();
		std::cout << "--- READING NEW ONLINE DEMANDS... --- " << std::endl;
		//instance.generateRandomDemandsFromFile();
		instance.generateRandomDemands(30);
		instance.displayNonRoutedDemands();
		std::cout << instance.getNbNonRoutedDemands() << " were generated." << std::endl;
		//CplexForm::setCount(0);
		int optimizationCounter = 0;
		instance.output(std::to_string(optimizationCounter));
		while(instance.getNbRoutedDemands() < instance.getNbDemands()){
			optimizationCounter++;
			CplexForm solver(instance);
			/************************************************/
			/*		            UPDATE MAPPING        		*/
			/************************************************/
			if (solver.getCplex().getStatus() == IloAlgorithm::Optimal){
				std::cout << "Update instance" << std::endl;
				solver.updateInstance(instance);
				instance.output(std::to_string(optimizationCounter));
				//instance.displayDetailedTopology();
			}
		}
		//instance.displayInstance(); 
		
	}
	catch(const std::invalid_argument& e){
		std::cerr << std::endl << "ERROR: Caught exception: A parameter file is required in the arguments. PLease run the program as \n./exec parameterFile.par\n" << std::endl;
	}/*
	catch(...){
		std::cerr << std::endl << "BIG FUCKING ERROR !!" << std::endl;
	}*/

	return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
