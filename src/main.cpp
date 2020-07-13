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
#include "ClockTime.h"

#include "cplexForm.h"
#include "YoussoufForm.h"
#include "subgradient.h"


using namespace lemon;

int main(int argc, char *argv[]) {
	try{
		ClockTime GLOBAL_TIME(ClockTime::getTimeNow());
		/********************************************************************/
		/* 						Get Parameter file 							*/
		/********************************************************************/
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
		
		/********************************************************************/
		/* 				For each file of demands, optimize it 				*/
		/********************************************************************/
		std::cout << "> Number of online demand files: " << input.getNbOnlineDemandFiles() << std::endl;
		for (int i = 0; i < input.getNbOnlineDemandFiles(); i++) {
			ClockTime OPTIMIZATION_TIME(ClockTime::getTimeNow());
			/********************************************************************/
			/* 						Create initial mapping 						*/
			/********************************************************************/
			std::cout << "--- READING INSTANCE... --- " << std::endl;
			Instance instance(input);
			std::cout << instance.getNbRoutedDemands() << " are present in the initial mapping." << std::endl;
			
			
			/********************************************************************/
			/* 					Define set of demands to be routed 				*/
			/********************************************************************/
			//instance.displayDetailedTopology();
			std::cout << "--- READING NEW ONLINE DEMANDS... --- " << std::endl;
			std::string nextFile = instance.getInput().getOnlineDemandFilesFromIndex(i);
			instance.generateDemandsFromFile(nextFile);
			//instance.generateRandomDemands(1);
			instance.displayNonRoutedDemands();
			std::cout << instance.getNbNonRoutedDemands() << " demands to be routed." << std::endl;
			
			
			/********************************************************************/
			/* 							Start optimizing 						*/
			/********************************************************************/
			int optimizationCounter = 0;
			std::string outputCode = getInBetweenString(nextFile, "/", ".") + "_" + std::to_string(optimizationCounter);
			//instance.output(outputCode);
			bool feasibility = true;
			bool lastIterationFeas = true;
			while(instance.getNextDemandToBeRoutedIndex() < instance.getNbDemands() && feasibility == true && (instance.getInput().getOptimizationTimeLimit() >= OPTIMIZATION_TIME.getTimeInSecFromStart())){
					
				/********************************************************************/
				/* 							Initialization	 						*/
				/********************************************************************/
				optimizationCounter++;
				outputCode = getInBetweenString(nextFile, "/", ".") + "_" + std::to_string(optimizationCounter);
				ClockTime ITERATION_TIME(ClockTime::getTimeNow());
				if ((instance.getInput().getOptimizationTimeLimit() - OPTIMIZATION_TIME.getTimeInSecFromStart()) < instance.getInput().getIterationTimeLimit()){
					instance.setTimeLimit(std::max(0, instance.getInput().getOptimizationTimeLimit() - (int)OPTIMIZATION_TIME.getTimeInSecFromStart()));
				}

				/********************************************************************/
				/* 								Solve	 							*/
				/********************************************************************/
				switch (instance.getInput().getChosenMethod()){
					case Input::METHOD_CPLEX:
					{
						FlowForm solver(instance);
						std::cout << "Status: " << solver.getCplex().getStatus() << std::endl;
						RSA::Status STATUS = solver.getStatus();
					
						if (STATUS == RSA::STATUS_FEASIBLE || STATUS == RSA::STATUS_OPTIMAL){
							solver.updateInstance(instance);
							int nextDemandToBeRouted = instance.getNextDemandToBeRoutedIndex() + solver.getNbDemandsToBeRouted();
							if (instance.getInput().isBlockingAllowed()){
								if (lastIterationFeas == false){
									nextDemandToBeRouted++;
								}
								instance.setNbDemandsAtOnce(input.getNbDemandsAtOnce());
							}
							else{
								if (lastIterationFeas == false){
									instance.setNbDemandsAtOnce(0);
								}
							}
							instance.setNextDemandToBeRoutedIndex(nextDemandToBeRouted);
							lastIterationFeas = true;
						}
						else{
							std::cout << "Decrease the number of demands to be treated." << std::endl;
							instance.decreaseNbDemandsAtOnce();
							lastIterationFeas = false;
							//instance.displayDetailedTopology();
						}

						if (instance.getInput().getNbDemandsAtOnce() <= 0){
							std::cout << "There is no room for an additional demand." << std::endl;
							if (instance.getInput().isBlockingAllowed()){
								lastIterationFeas = true;
								int nextDemandToBeRouted = instance.getNextDemandToBeRoutedIndex() + 1;
								instance.setNextDemandToBeRoutedIndex(nextDemandToBeRouted);
								instance.setNbDemandsAtOnce(input.getNbDemandsAtOnce());
							}
							else{
								feasibility = false;
							}
						}
						break;
					}
					case Input::METHOD_SUBGRADIENT:
					{
						Subgradient solver(instance);
						RSA::Status STATUS = solver.getStatus();
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
						//instance.output(outputCode);
						break;
					}
					case Input::METHOD_YOUSSOUF:
					{
						YoussoufForm solver(instance);	
							
						std::cout << "Status: " << solver.getCplex().getStatus() << std::endl;
						RSA::Status STATUS = solver.getStatus();
					
						if (STATUS == RSA::STATUS_FEASIBLE || STATUS == RSA::STATUS_OPTIMAL){
							solver.updateInstance(instance);
							int nextDemandToBeRouted = instance.getNextDemandToBeRoutedIndex() + solver.getNbDemandsToBeRouted();
							if (instance.getInput().isBlockingAllowed()){
								if (lastIterationFeas == false){
									nextDemandToBeRouted++;
								}
								instance.setNbDemandsAtOnce(input.getNbDemandsAtOnce());
							}
							else{
								if (lastIterationFeas == false){
									instance.setNbDemandsAtOnce(0);
								}
							}
							instance.setNextDemandToBeRoutedIndex(nextDemandToBeRouted);
							lastIterationFeas = true;
						}
						else{
							std::cout << "Decrease the number of demands to be treated." << std::endl;
							instance.decreaseNbDemandsAtOnce();
							lastIterationFeas = false;
							//instance.displayDetailedTopology();
						}

						if (instance.getInput().getNbDemandsAtOnce() <= 0){
							std::cout << "There is no room for an additional demand." << std::endl;
							if (instance.getInput().isBlockingAllowed()){
								lastIterationFeas = true;
								int nextDemandToBeRouted = instance.getNextDemandToBeRoutedIndex() + 1;
								instance.setNextDemandToBeRoutedIndex(nextDemandToBeRouted);
								instance.setNbDemandsAtOnce(input.getNbDemandsAtOnce());
							}
							else{
								feasibility = false;
							}
						}
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
					
				std::cout << "Time taken by iteration is : ";
				std::cout << std::fixed  << ITERATION_TIME.getTimeInSecFromStart() << std::setprecision(9); 
				std::cout << " sec" << std::endl; 
			}
			if (instance.getInput().getChosenOutputLvl() >= Input::OUTPUT_LVL_NORMAL){
				outputCode = getInBetweenString(nextFile, "/", ".");
				instance.output(outputCode + "_FINAL");
				instance.outputLogResults(outputCode);
			}
			instance.displayAllDemands();
			std::cout << "Time taken by optimization is : ";
			std::cout << std::fixed  << OPTIMIZATION_TIME.getTimeInSecFromStart() << std::setprecision(9); 
			std::cout << " sec" << std::endl; 
		}

		std::cout << "Total time taken by program is : ";
		std::cout << std::fixed  << GLOBAL_TIME.getTimeInSecFromStart() << std::setprecision(9); 
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
