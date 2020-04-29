#include <ilcplex/ilocplex.h>
ILOSTLBEGIN
// cplex10Test.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <bits/stdc++.h> 
#include <chrono> 

#include <lemon/list_graph.h>
#include <lemon/concepts/graph.h>

#include "PhysicalLink.h"
#include "Instance.h"
#include "input.h"

#include "cplexForm.h"
#include "subgradient.h"

using namespace lemon;

int main(int argc, char *argv[]) {
	try{
		std::string parameterFile;
		if (argc != 2){
			std::cerr << "A parameter file is required in the arguments. PLease run the program as \n./exec parameterFile.par\n";
			throw std::invalid_argument( "did not receive an argument" );
		}
		else{
			parameterFile = argv[1];
		}
		std::cout << "PARAMETER FILE: " << parameterFile << std::endl;
		Input input(parameterFile);
		
		std::cout << "> Number of online demand files: " << input.getNbOnlineDemandFiles() << std::endl;
		for (int i = 0; i < input.getNbOnlineDemandFiles(); i++) {
			
			std::cout << "--- READING INSTANCE... --- " << std::endl;
			Instance instance(input);

			std::cout << "--- CREATING INITIAL MAPPING... --- " << std::endl;
			instance.createInitialMapping();
			std::cout << instance.getNbRoutedDemands() << " demands were routed." << std::endl;
			
			//instance.displayDetailedTopology();
			std::cout << "--- READING NEW ONLINE DEMANDS... --- " << std::endl;
			std::string nextFile = instance.getInput().getOnlineDemandFilesFromIndex(i);
			instance.generateRandomDemandsFromFile(nextFile);
			//instance.generateRandomDemands(1);
			instance.displayNonRoutedDemands();
			std::cout << instance.getNbNonRoutedDemands() << " demands were generated." << std::endl;
			//CplexForm::setCount(0);
			int optimizationCounter = 0;
			std::string outputCode = getInBetweenString(nextFile, "/", ".") + "_" + std::to_string(optimizationCounter);
			instance.output(outputCode);
			bool feasibility = true;
			while(instance.getNbRoutedDemands() < instance.getNbDemands() && feasibility == true){
				optimizationCounter++;
				outputCode = getInBetweenString(nextFile, "/", ".") + "_" + std::to_string(optimizationCounter);
				std::chrono::_V2::system_clock::time_point start = std::chrono::high_resolution_clock::now();

				switch (instance.getInput().getChosenMethod()){
				case Input::METHOD_CPLEX:
					{
						CplexForm solver(instance);			
						if (solver.getCplex().getStatus() == IloAlgorithm::Optimal){
							solver.updateInstance(instance);
							instance.output(outputCode);
							//instance.displayDetailedTopology();
						}
						else{
							feasibility = false;
							instance.outputLogResults(getInBetweenString(nextFile, "/", "."));
						}
						break;
					}
				case Input::METHOD_SUBGRADIENT:
					{
						Subgradient sub(instance);
						sub.updateInstance(instance);
						instance.output(outputCode);
						break;
					}
				default:
					{
						std::cerr << "The parameter \'chosenMethod\' is invalid. " << std::endl;
						throw std::invalid_argument( "did not receive an argument" );
						break;
					}
					
				}
			
				std::chrono::_V2::system_clock::time_point end = std::chrono::high_resolution_clock::now();
				double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); 
				time_taken *= 1e-9; 
		
				std::cout << "Time taken by program is : " << std::fixed  << time_taken << std::setprecision(9); 
				std::cout << " sec" << std::endl; 
			}
		
		}
		

		//instance.displayInstance();
	}
	catch(const std::invalid_argument& e){
		std::cerr << std::endl << "ERROR: Caught exception." << std::endl;
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
