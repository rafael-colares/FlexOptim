#include <ilcplex/ilocplex.h>
ILOSTLBEGIN
// cplex10Test.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

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
		/********************************************************************/
		/* 						Get Parameter file 							*/
		/********************************************************************/
		
		std::chrono::_V2::system_clock::time_point GLOBAL_START = std::chrono::high_resolution_clock::now();
		std::string parameterFile;
		if (argc < 2){
			std::cerr << "A parameter file is required in the arguments. PLease run the program as \n./exec parameterFile.par\n";
			throw std::invalid_argument( "did not receive an argument" );
		}
		else{
			parameterFile = argv[1];
		}
		std::cout << "PARAMETER FILE: " << parameterFile << std::endl;
		Input input(parameterFile);
		
		// For each set of online demands, optimize it!
		std::cout << "> Number of online demand files: " << input.getNbOnlineDemandFiles() << std::endl;
		for (int i = 0; i < input.getNbOnlineDemandFiles(); i++) {
				
			/********************************************************************/
			/* 						Create initial mapping 						*/
			/********************************************************************/
			std::cout << "--- READING INSTANCE... --- " << std::endl;
			Instance instance(input);
			std::cout << instance.getNbRoutedDemands() << " demands were routed." << std::endl;
			
			
			/********************************************************************/
			/* 					Define set of demands to be routed 				*/
			/********************************************************************/
			//instance.displayDetailedTopology();
			std::cout << "--- READING NEW ONLINE DEMANDS... --- " << std::endl;
			std::string nextFile = instance.getInput().getOnlineDemandFilesFromIndex(i);
			instance.generateDemandsFromFile(nextFile);
			//instance.generateRandomDemands(1);
			instance.displayNonRoutedDemands();
			std::cout << instance.getNbNonRoutedDemands() << " demands were generated." << std::endl;
			
			
			/********************************************************************/
			/* 							Start optimizing 						*/
			/********************************************************************/
			int optimizationCounter = 0;
			std::string outputCode = getInBetweenString(nextFile, "/", ".") + "_" + std::to_string(optimizationCounter);
			//instance.output(outputCode);
			bool feasibility = true;
			while(instance.getNbRoutedDemands() < instance.getNbDemands() && feasibility == true){
				optimizationCounter++;
				outputCode = getInBetweenString(nextFile, "/", ".") + "_" + std::to_string(optimizationCounter);
				std::chrono::_V2::system_clock::time_point start = std::chrono::high_resolution_clock::now();

				switch (instance.getInput().getChosenMethod()){
				case Input::METHOD_CPLEX:
					{
						CplexForm solver(instance);	
						
						std::cout << "Status: " << solver.getCplex().getStatus() << std::endl;
						RSA::Status STATUS = solver.getStatus();
						if (STATUS == RSA::STATUS_ERROR){
							std::cout << "Got error on CPLEX status." << std::endl;
							exit(0);
						}
						if (STATUS == RSA::STATUS_FEASIBLE || STATUS == RSA::STATUS_OPTIMAL){
							solver.updateInstance(instance);
						}
						else{
							std::cout << "Decrease the number of demands to be treated." << std::endl;
							instance.decreaseNbDemandsAtOnce();
							//instance.displayDetailedTopology();
						}

						if (instance.getInput().getNbDemandsAtOnce() <= 0){
							std::cout << "There is no room for an additional demand." << std::endl;
							feasibility = false;
						}
						break;
					}
				case Input::METHOD_SUBGRADIENT:
					{
						Subgradient sub(instance);
						if (sub.getIsFeasible() ==  true){
							sub.updateInstance(instance);
						}
						else{
							feasibility = false;
						}
						//instance.output(outputCode);
						break;
					}
				default:
					{
						std::cerr << "The parameter \'chosenMethod\' is invalid. " << std::endl;
						throw std::invalid_argument( "did not receive an argument" );
						break;
					}
					
				}
				if (instance.getInput().getChosenOutputLvl() == Input::OUTPUT_LVL_DETAILED){
					instance.output(outputCode);
				}
				std::chrono::_V2::system_clock::time_point end = std::chrono::high_resolution_clock::now();
				double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); 
				time_taken *= 1e-9; 
		
				std::cout << "Time taken by program is : " << std::fixed  << time_taken << std::setprecision(9); 
				std::cout << " sec" << std::endl; 
			}
			if (instance.getInput().getChosenOutputLvl() >= Input::OUTPUT_LVL_NORMAL){
				outputCode = getInBetweenString(nextFile, "/", ".");
				instance.output(outputCode + "_FINAL");
				instance.outputLogResults(outputCode);
			}
		}
		
		std::chrono::_V2::system_clock::time_point GLOBAL_END = std::chrono::high_resolution_clock::now();
		double GLOBAL_TIME = std::chrono::duration_cast<std::chrono::nanoseconds>(GLOBAL_END - GLOBAL_START).count(); 
		GLOBAL_TIME *= 1e-9; 

		std::cout << "Total time taken by program is : " << std::fixed  << GLOBAL_TIME << std::setprecision(9); 
		std::cout << " sec" << std::endl; 
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
