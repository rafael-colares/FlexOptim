#include <ilcplex/ilocplex.h>
ILOSTLBEGIN
// cplex10Test.cpp : This file contains the 'main' function. Program execution begins and ends there.

#include <bits/stdc++.h> 
#include <chrono> 

#include <lemon/list_graph.h>
#include <lemon/concepts/graph.h>

#include "tools/clockTime.h"
#include "topology/instance.h"
#include "solver/solverFactory.h"

#include "lagrangian/solver/lagSolverFactory.h"


using namespace lemon;

void createFile(std::string parameterFile,std::string linkfile,std::string demandfolder,int nbdemands,
                int rl,int lagRelax,int relaxMethod,int nodeMethod,int solver, int lagFormulation,
                int heuristic,int projection,int warmstart,int alternativeStop,int directionMethod,
                double crowderParam,double carmeriniParam,double lagrangianLambda_zero,int nbIterationsWithoutImprovement,
                int maxNbIterations,int obj, std::string outputFolder){

    double lagrangianMultiplier_zero = 0.0;
    double time = 18000;

    std::ofstream fichier(parameterFile.c_str());
    if (!fichier.fail()) {

        fichier << "******* Input File Paths *******" << std::endl;;
        fichier << "topologyFile=" << linkfile <<""<< std::endl;
        fichier << "initialMappingDemandFile=" << std::endl;
        fichier << "initialMappingAssignmentFile=" << std::endl;
        fichier << "demandToBeRoutedFolder="<< demandfolder << "" << std::endl;
        fichier << std::endl;

        fichier << "******* GNPY parameters *******" << std::endl;
        fichier << "GNPY_activation=0 "<< std::endl;
        fichier << "GNPY_topologyFile=../oopt-gnpy/gnpy/example-data/spain_topo.json " << std::endl;
        fichier << "GNPY_equipmentFile=../oopt-gnpy/gnpy/example-data/spain_eqpt_config.json " << std::endl;
        fichier << std::endl;

        fichier << "******* Formulation parameters *******" << std::endl;
        fichier << "nbDemandsAtOnce=" << nbdemands << " " << std::endl;
        fichier << "formulation=0 " << std::endl;
        fichier << "userCuts=0 " << std::endl;
        fichier << "obj=" << obj<<"" << std::endl;
        fichier << "allowBlocking=0 " << std::endl;
        fichier << "hopPenalty=0 " << std::endl;
        fichier << "partitionPolicy=0 " << std::endl;
        fichier << "partitionLoad=4 " << std::endl;
        fichier << "partitionSlice=15 " << std::endl;
        fichier << std::endl;

        fichier << "******* Optimization parameters *******" << std::endl;
        fichier << "solver=" <<solver << " "<< std::endl;
        fichier << "method=" << nodeMethod << " " << std::endl;
        fichier << "preprocessingLevel=0 " << std::endl;
        fichier << "linearRelaxation="<< rl << " " << std::endl;
        fichier << "relaxMethod="<< relaxMethod << " " << std::endl;
        fichier << "lagrangianRelaxation="<< lagRelax << " " << std::endl;
        fichier << std::endl;

        fichier << "******* Execution parameters *******" << std::endl;
        fichier << "outputPath=../Parameters/Instances/Spain_N5/InitialMappingMet1Slice30/SimpleTest/Output " << std::endl;
        fichier << "outputLevel=2 " << std::endl;
        fichier << "nbSlicesInOutputFile=320 " << std::endl;
        fichier << "globalTimeLimit=70000 " << std::endl;
        fichier << "timeLimit="<< time << " " << std::endl;
        fichier << std::endl;

        fichier << "******* Fields below are reserved for team LIMOS ********" << std::endl;
        fichier << "lagrangianMultiplier_zero=" << lagrangianMultiplier_zero<<" " << std::endl;
        fichier << "lagrangianLambda_zero=" << lagrangianLambda_zero <<" "<< std::endl;
        fichier << "nbIterationsWithoutImprovement=" << nbIterationsWithoutImprovement<< " " << std::endl;
        fichier << "maxNbIterations=" << maxNbIterations << " "<< std::endl;
        fichier << std::endl;

        fichier << "lagFormulation=" << lagFormulation << " " << std::endl;
        fichier << "heuristic=" << heuristic << " " << std::endl;
        fichier << std::endl;

        fichier << "directionMethod=" << directionMethod << " " << std::endl;
        fichier << "crowderParam=" << crowderParam << " " << std::endl;
        fichier << "carmeriniParam=" << carmeriniParam << " " << std::endl;
        fichier << std::endl;

        fichier << "projection=" << projection << " " << std::endl;
        fichier << "alternativeStop=" << alternativeStop << " " << std::endl;
        fichier << "warmstart=" << warmstart << " " << std::endl;
        fichier << std::endl;

        fichier << "lagOutputPath=" << outputFolder << " " << std::endl;
    }
    fichier.close();
}


int main(int argc, char *argv[]) {
	/**************************************************************************************************************************/
	/* 						                             Get Parameter file 							                      */
	/**************************************************************************************************************************/
	std::string parameterFile;
	if (argc < 2){
		std::cerr << "A parameter file is required in the arguments. Please run the program as \n./exec parameterFile.par topology objective\n";
		throw std::invalid_argument( "did not receive an argument" );
	}
	else{
		parameterFile = argv[1];
	}
	std::cout << "PARAMETER FILE: " << parameterFile << std::endl;
    std::string topology;
    if (argc < 3){
        std::cerr << "A topology file is required in the arguments. Please run the program as \n./exec parameterFile.par topology objective\n";
		throw std::invalid_argument( "did not receive an argument" );
    }else{
        topology = argv[2];
    }
    int obj;
    if (argc < 4){
        std::cerr << "A objetive function is required in the arguments. Please run the program as \n./exec parameterFile.par topology objective\n";
		throw std::invalid_argument( "did not receive an argument" );
    }else{
        obj = std::atoi(argv[3]);
    }
    
    /*************************************************************************************************************************/
	/* 						                              Instances to test 							                     */
	/*************************************************************************************************************************/
    int m = 7; // instances

    std::string generalFolder              = "../Parameters/Instances/LagrangianTests/";
    std::string instances[7]               = {"10demands/","20demands/","30demands/","40demands/","50demands/","60demands/","70demands/"};

    std::string linkfile[7];
    for(int j=0;j<m;j++){
        linkfile[j]= generalFolder + topology + instances[j] + "Link.csv";
    }

    std::string demandfolders[7];
    for(int j=0;j<m;j++){
        demandfolders[j] = generalFolder + topology + instances[j] + "Demands1";
    }

    int numdemands[7] = {10,20,30,40,50,60,70};

    /************************  Parameters *************************/
    // We always use the shortest path heuristic.
    int heuristic = 0;

    // We always analyse the Flow formulation.
    int lagFormulation = 0;

    // Default configuration for the Subgradient
    int projection = 0;
    int warmstart = 0;
    int alternativeStop = 0;
    int directionMethod = 0;
    double crowderParam = 0.0;
    double carmeriniParam = 0.0;
    double lagrangianLambda_zero = 2.0;
    int nbIterationsWithoutImprovement = 10;
    int maxNbIterations = 5000;

    // Changing parameters
    int rl, relaxMethod, solver, nodeMethod, lagRelax; 

    /************************ File with the responses *************************/
    std::string outputFolder = generalFolder + topology + "outputs/";

    std::string nom_fichier = outputFolder + "obj_"+ std::to_string(obj) +"_general.csv";
    std::ofstream fichier(nom_fichier);
    std::string delimiter = ";";

    fichier << std::endl << std::endl;
    fichier << " Slices ; Demands ; MIP-VOL-UB ; MIP-VOL-LB ; MIP-VOL-GAP ; MIP-VOL-Tree-Size ; MIP-Time; " << std::endl;

    for(int j=0;j<m;j++){
        rl = 0; relaxMethod = 4; solver = 1; nodeMethod = 2; lagRelax = 0;
        createFile(parameterFile,linkfile[j],demandfolders[j],numdemands[j],
                                    rl,lagRelax,relaxMethod,nodeMethod,solver,lagFormulation,
                                    heuristic,projection,warmstart,alternativeStop,directionMethod,crowderParam,carmeriniParam,
                                    lagrangianLambda_zero,nbIterationsWithoutImprovement,maxNbIterations,obj,outputFolder);
        Input input(parameterFile);

        std::cout << "--- READING INSTANCE... --- " << std::endl;
        Instance instance(input);

        std::cout << "--- READING NEW ONLINE DEMANDS... --- " << std::endl;
        std::string nextFile = instance.getInput().getDemandToBeRoutedFilesFromIndex(0);
        instance.generateDemandsFromFile(nextFile);

        std::cout << "Solving with MIP-Vol" << std::endl;
        SolverFactory factory;
        AbstractSolver *solver = factory.createSolver(instance);
        solver->solve();

        fichier << instance.getMaxSlice() << delimiter;
        fichier << instance.getNbDemands() << delimiter;
        fichier << solver->getUpperBound() << delimiter;
        fichier << solver->getLowerBound() << delimiter;
        fichier << solver->getMipGap() << delimiter;
        fichier << solver->getTreeSize() << delimiter;
        fichier << solver->getDurationTime() << delimiter << std::endl;

        std::cout << "Mip-Vol completed" << std::endl;

        delete solver;
    }

    fichier << std::endl << std::endl;
    fichier << " Slices ; Demands ; MIP-UB ; MIP-LB ; MIP-GAP ; MIP-Tree-Size ; MIP-Time; " << std::endl;

    for(int j=0;j<m;j++){
        rl = 0; relaxMethod = 4; solver = 1; nodeMethod = 0; lagRelax = 0;
        createFile(parameterFile,linkfile[j],demandfolders[j],numdemands[j],
                                    rl,lagRelax,relaxMethod,nodeMethod,solver,lagFormulation,
                                    heuristic,projection,warmstart,alternativeStop,directionMethod,crowderParam,carmeriniParam,
                                    lagrangianLambda_zero,nbIterationsWithoutImprovement,maxNbIterations,obj,outputFolder);
        Input input_mip(parameterFile);

        std::cout << "--- READING INSTANCE... --- " << std::endl;
        Instance instance_mip(input_mip);

        std::cout << "--- READING NEW ONLINE DEMANDS... --- " << std::endl;
        std::string nextFile_mip = instance_mip.getInput().getDemandToBeRoutedFilesFromIndex(0);
        instance_mip.generateDemandsFromFile(nextFile_mip);

        std::cout << "Solving with MIP-Relax" << std::endl;
        SolverFactory factory_mip;
        AbstractSolver *solver_mip = factory_mip.createSolver(instance_mip);
        solver_mip->solve();

        fichier << instance_mip.getMaxSlice() << delimiter;
        fichier << instance_mip.getNbDemands() << delimiter;
        fichier << solver_mip->getUpperBound() << delimiter;
        fichier << solver_mip->getLowerBound() << delimiter;
        fichier << solver_mip->getMipGap() << delimiter;
        fichier << solver_mip->getTreeSize() << delimiter;
        fichier << solver_mip->getDurationTime() << delimiter << std::endl;

        std::cout << "Mip-Relax completed" << std::endl;

        delete solver_mip;
    }

    fichier << std::endl;
    fichier.close();
    
	return 0;
}

