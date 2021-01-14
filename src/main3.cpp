#include <ilcplex/ilocplex.h>
ILOSTLBEGIN
// cplex10Test.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <bits/stdc++.h> 
#include <chrono> 

#include <lemon/list_graph.h>
#include <lemon/concepts/graph.h>

#include "tools/clockTime.h"
#include "topology/instance.h"
#include "solver/solverFactory.h"

#include "lagrangian/solver/lagSolverFactory.h"

//#include "YoussoufForm.h"
//#include "subgradient.h"


using namespace lemon;

void createFile(std::string parameterFile,std::string linkfile,std::string demandfolder,int nbdemands,int rl,int lagMethod,int lagFormulation,int heuristic,int projection,int warmstart,int alternativeStop,int directionMethod,double crowderParam,double carmeriniParam,double lagrangianLambda_zero,int nbIterationsWithoutImprovement, int maxNbIterations){

    int obj = 1;
    double lagrangianMultiplier_zero = 0.0;

    std::ofstream fichier(parameterFile.c_str());
    if (!fichier.fail()) {
        fichier << "******* Input File Paths *******" << std::endl;;
        fichier << "topologyFile=" << linkfile <<" "<< std::endl;
        fichier << "initialMappingDemandFile= " << std::endl;
        fichier << "initialMappingAssignmentFile= " << std::endl;
        fichier << "demandToBeRoutedFolder="<< demandfolder << " " << std::endl;
        fichier << std::endl;
        fichier << "******* GNPY parameters *******" << std::endl;
        fichier << "GNPY_activation=0 "<< std::endl;
        fichier << "GNPY_topologyFile=../oopt-gnpy/gnpy/example-data/spain_topo.json " << std::endl;
        fichier << "GNPY_equipmentFile=../oopt-gnpy/gnpy/example-data/spain_eqpt_config.json " << std::endl;
        fichier << std::endl;
        fichier << "******* Formulation parameters *******" << std::endl;
        fichier << "nbDemandsAtOnce=" << nbdemands << " " << std::endl;
        fichier << "formulation=0 " << std::endl;
        fichier << "obj=" <<obj<<" " << std::endl;
        fichier << "allowBlocking=0 " << std::endl;
        fichier << "hopPenalty=0 " << std::endl;
        fichier << "partitionPolicy=0 " << std::endl;
        fichier << "partitionLoad=4 " << std::endl;
        fichier << "partitionSlice=15 " << std::endl;
        fichier << std::endl;
        fichier << "******* Optimization parameters *******" << std::endl;
        fichier << "solver=0 " << std::endl;
        fichier << "method=0 " << std::endl;
        fichier << "preprocessingLevel=2 " << std::endl;
        fichier << "linearRelaxation="<< rl << " " << std::endl;
        fichier << std::endl;
        fichier << "******* Execution parameters *******" << std::endl;
        fichier << "outputPath=../Parameters/Instances/Spain_N5/InitialMappingMet1Slice30/SimpleTest/Output " << std::endl;
        fichier << "outputLevel=2 " << std::endl;
        fichier << "nbSlicesInOutputFile=320 " << std::endl;
        fichier << "globalTimeLimit=70000 " << std::endl;
        fichier << "timeLimit=7200 " << std::endl;
        fichier << std::endl;
        fichier << "******* Fields below are reserved for team LIMOS ********" << std::endl;
        fichier << "lagrangianMultiplier_zero=" << lagrangianMultiplier_zero<<" " << std::endl;
        fichier << "lagrangianLambda_zero=" << lagrangianLambda_zero <<" "<< std::endl;
        fichier << "nbIterationsWithoutImprovement=" << nbIterationsWithoutImprovement<< " " << std::endl;
        fichier << "maxNbIterations=" << maxNbIterations << " "<< std::endl;
        fichier << std::endl;
        fichier << "lagMethod=" << lagMethod << " " << std::endl;
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
    }
    fichier.close();
}

int main(int argc, char *argv[]) {
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

    /********************************************************************/
	/* 						Instances to test 							*/
	/********************************************************************/
    int n = 10;
    int m = 4;
    std::string aux              = "../Parameters/MeusTestes/NSF/";
    std::string aux2[m]          = {"10demands/","20demands/","30demands/","40demands/"};//,"50demands/","60demands/"};
    std::string aux3[n]          = {"Demands1","Demands2","Demands3","Demands4","Demands5","Demands6","Demands7","Demands8","Demands9","Demands10"};
    std::string linkfile         = "../Parameters/MeusTestes/NSF/Link.csv";
    std::string demandfolders[m][n];
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            demandfolders[i][j] = aux + aux2[i] + aux3[j];
        }
    }
    int numdemands[m] = {10,20,30,40};//,50,60};

    for(int i=3;i<m;i++){
        /************************ File with the responses *************************/
        std::string nom_fichier = aux + aux2[i] + "general-teste2-faltou.csv";
        std::ofstream fichier(nom_fichier);
        
        fichier << "MIP-UB;MIP-LB;MIP-GAP;MIP-Time;;RELAX-OBJ;RELAX-Time;;";
        fichier << "SUB-FLOW-UB;SUB-FLOW-LB;SUB-FLOW-Iterations;SUB-FLOW-Stop Creterion;SUB-FLOW-Time;;";
        fichier << "SUB-OVERLAP-UB;SUB-OVERLAP-LB;SUB-OVERLAP-Iterations;SUB-OVERLAP-Stop Creterion;SUB-OVERLAP-Time;;";
        fichier << "VOL-FLOW-UB;VOL-FLOW-LB;VOL-FLOW-Iterations;VOL-FLOW-Stop Creterion;VOL-FLOW-Time;;";
        fichier << "VOL-OVERLAP-UB;VOL-OVERLAP-LB;VOL-OVERLAP-Iterations;VOL-OVERLAP-Stop Creterion;VOL-OVERLAP-Time;;";
        fichier <<std::endl; 

        /************************ File with the responses *************************/
        // We always use the shortest path heuristic
        int heuristic = 0;

        // Default configuration for the Subgradient
        int projection = 0;
        int warmstart = 0;
        int alternativeStop = 0;
        int directionMethod = 0;
        double crowderParam = 0.0;
        double carmeriniParam = 0.0;
        double lagrangianLambda_zero = 2.0;
        int nbIterationsWithoutImprovement = 10;
        int maxNbIterations = 150;

        // Changing parameters
        int rl;
        int lagMethod;
        int lagFormulation;

        for(int j=5;j<6;j++){
        
            rl = 0; lagMethod = 0; lagFormulation = 0;
            createFile(parameterFile,linkfile,demandfolders[i][j],numdemands[i],rl,lagMethod,lagFormulation,heuristic,projection,warmstart,alternativeStop,directionMethod,crowderParam,carmeriniParam,lagrangianLambda_zero,nbIterationsWithoutImprovement,maxNbIterations);
            Input input(parameterFile);

            std::cout << "--- READING INSTANCE... --- " << std::endl;
            Instance instance(input);

            std::cout << "--- READING NEW ONLINE DEMANDS... --- " << std::endl;
            std::string nextFile = instance.getInput().getDemandToBeRoutedFilesFromIndex(0);
            instance.generateDemandsFromFile(nextFile);

            /********************************************************************/
            /* 				        Solve - MIP SOLUTION	 					*/
            /********************************************************************/
            std::cout << "Solving with MIP-Cplex" << std::endl;
            ClockTime MIP_OPTIMIZATION_TIME(ClockTime::getTimeNow());
            SolverFactory factory;
            AbstractSolver *solver = factory.createSolver(instance);
            solver->solve();
            fichier << solver->getUpperBound() << ";" << solver->getLowerBound() << ";" << solver->getMipGap() << ";" << std::fixed  << MIP_OPTIMIZATION_TIME.getTimeInSecFromStart() << std::setprecision(9);
            std::cout << "Time taken by optimization is : ";
            std::cout << std::fixed  << MIP_OPTIMIZATION_TIME.getTimeInSecFromStart() << std::setprecision(9);
            std::cout << " sec" << std::endl; 
        
            /*********************************************************************/
            /*********************************************************************/
            
            rl = 1; lagMethod = 0; lagFormulation = 0;
            createFile(parameterFile,linkfile,demandfolders[i][j],numdemands[i],rl,lagMethod,lagFormulation,heuristic,projection,warmstart,alternativeStop,directionMethod,crowderParam,carmeriniParam,lagrangianLambda_zero,nbIterationsWithoutImprovement,maxNbIterations);
            Input input2(parameterFile);

            std::cout << "--- READING INSTANCE... --- " << std::endl;
            Instance instance2(input2);
            
            std::cout << "--- READING NEW ONLINE DEMANDS... --- " << std::endl;
            std::string nextFile2 = instance2.getInput().getDemandToBeRoutedFilesFromIndex(0);
            instance2.generateDemandsFromFile(nextFile2);

            /********************************************************************/
            /* 				        Solve - RELAXATION	 					*/
            /********************************************************************/
            std::cout << "Solving with RELAX-Cplex" << std::endl;
            ClockTime RELAX_OPTIMIZATION_TIME(ClockTime::getTimeNow());
            SolverFactory factory2;
            AbstractSolver *solver2 = factory2.createSolver(instance2);
            solver2->solve();
            fichier << ";;" << solver2->getUpperBound() << ";" << std::fixed  << RELAX_OPTIMIZATION_TIME.getTimeInSecFromStart() << std::setprecision(9);
            std::cout << "Time taken by optimization is : ";
            std::cout << std::fixed  << RELAX_OPTIMIZATION_TIME.getTimeInSecFromStart() << std::setprecision(9); 
            std::cout << " sec" << std::endl; 
        
            /*********************************************************************/
            /*********************************************************************/
            
            rl = 0; lagMethod = 0; lagFormulation = 0; maxNbIterations = 150;
            createFile(parameterFile,linkfile,demandfolders[i][j],numdemands[i],rl,lagMethod,lagFormulation,heuristic,projection,warmstart,alternativeStop,directionMethod,crowderParam,carmeriniParam,lagrangianLambda_zero,nbIterationsWithoutImprovement,maxNbIterations);
            Input input3(parameterFile);

            std::cout << "--- READING INSTANCE... --- " << std::endl;
            Instance instance3(input3);
            
            std::cout << "--- READING NEW ONLINE DEMANDS... --- " << std::endl;
            std::string nextFile3 = instance3.getInput().getDemandToBeRoutedFilesFromIndex(0);
            instance3.generateDemandsFromFile(nextFile3);
        
            /********************************************************************/
            /* 		      Solve - SUBGRADIENT LAGFLOW FORMULATION               */
            /********************************************************************/
            std::cout << "Solving with subgradient and flow formulation" << std::endl;
            ClockTime LAG1_OPTIMIZATION_TIME(ClockTime::getTimeNow());
            lagSolverFactory lagfactory;
            AbstractLagSolver *lagsolver = lagfactory.createSolver(instance3);
            lagsolver->run();
            fichier << ";;" << lagsolver->getUB() << ";" << lagsolver->getLB() << ";" << lagsolver->getIteration() << ";" << lagsolver->getStop()<< ";" << std::fixed  << LAG1_OPTIMIZATION_TIME.getTimeInSecFromStart() << std::setprecision(9);
            std::cout << "Time taken by optimization is : ";
            std::cout << std::fixed  << LAG1_OPTIMIZATION_TIME.getTimeInSecFromStart() << std::setprecision(9); 
            std::cout << " sec" << std::endl; 
        
            /*********************************************************************/
            /*********************************************************************/

            rl = 0; lagMethod = 0; lagFormulation = 1; maxNbIterations = 300;
            createFile(parameterFile,linkfile,demandfolders[i][j],numdemands[i],rl,lagMethod,lagFormulation,heuristic,projection,warmstart,alternativeStop,directionMethod,crowderParam,carmeriniParam,lagrangianLambda_zero,nbIterationsWithoutImprovement,maxNbIterations);
            Input input4(parameterFile);

            std::cout << "--- READING INSTANCE... --- " << std::endl;
            Instance instance4(input4);
		
            std::cout << "--- READING NEW ONLINE DEMANDS... --- " << std::endl;
            std::string nextFile4 = instance4.getInput().getDemandToBeRoutedFilesFromIndex(0);
            instance4.generateDemandsFromFile(nextFile4);

            /********************************************************************/
            /* 		Solve - SUBGRADIENT LAGNONOVERLAP FORMULATION      	 		*/
            /********************************************************************/
            std::cout << "Solving with subgradient non overlap formualtion" << std::endl;
            ClockTime LAG2_OPTIMIZATION_TIME(ClockTime::getTimeNow());
            lagSolverFactory lagfactory2;
            AbstractLagSolver *lagsolver2 = lagfactory2.createSolver(instance4);
            lagsolver2->run();
            fichier << ";;" << lagsolver2->getUB() << ";" << lagsolver2->getLB() << ";" << lagsolver2->getIteration() << ";" << lagsolver2->getStop()<< ";" << std::fixed  << LAG2_OPTIMIZATION_TIME.getTimeInSecFromStart() << std::setprecision(9);
            std::cout << "Time taken by optimization is : ";
            std::cout << std::fixed  << LAG2_OPTIMIZATION_TIME.getTimeInSecFromStart() << std::setprecision(9); 
            std::cout << " sec" << std::endl; 
        
            /*********************************************************************/
            /*********************************************************************/
        
            rl = 0; lagMethod = 1; lagFormulation = 0; maxNbIterations = 150;
            createFile(parameterFile,linkfile,demandfolders[i][j],numdemands[i],rl,lagMethod,lagFormulation,heuristic,projection,warmstart,alternativeStop,directionMethod,crowderParam,carmeriniParam,lagrangianLambda_zero,nbIterationsWithoutImprovement,maxNbIterations);
            Input input5(parameterFile);

            std::cout << "--- READING INSTANCE... --- " << std::endl;
            Instance instance5(input5);
            
            std::cout << "--- READING NEW ONLINE DEMANDS... --- " << std::endl;
            std::string nextFile5 = instance5.getInput().getDemandToBeRoutedFilesFromIndex(0);
            instance5.generateDemandsFromFile(nextFile5);
        
            /********************************************************************/
            /* 		      Solve - VOLUME LAGFLOW FORMULATION               */
            /********************************************************************/
            std::cout << "Solving with volume and flow formulation" << std::endl;
            ClockTime LAG3_OPTIMIZATION_TIME(ClockTime::getTimeNow());
            lagSolverFactory lagfactory3;
            AbstractLagSolver *lagsolver3 = lagfactory3.createSolver(instance5);
            lagsolver3->run();
            fichier << ";;" << lagsolver3->getUB() << ";" << lagsolver3->getLB() << ";" << lagsolver3->getIteration() << ";" << lagsolver3->getStop()<< ";" << std::fixed  << LAG3_OPTIMIZATION_TIME.getTimeInSecFromStart() << std::setprecision(9);
            std::cout << "Time taken by optimization is : ";
            std::cout << std::fixed  << LAG3_OPTIMIZATION_TIME.getTimeInSecFromStart() << std::setprecision(9); 
            std::cout << " sec" << std::endl; 
        
            /*********************************************************************/
            /*********************************************************************/

            rl = 0; lagMethod = 1; lagFormulation = 1; maxNbIterations = 300;
            createFile(parameterFile,linkfile,demandfolders[i][j],numdemands[i],rl,lagMethod,lagFormulation,heuristic,projection,warmstart,alternativeStop,directionMethod,crowderParam,carmeriniParam,lagrangianLambda_zero,nbIterationsWithoutImprovement,maxNbIterations);
            Input input6(parameterFile);

            std::cout << "--- READING INSTANCE... --- " << std::endl;
            Instance instance6(input6);
            
            std::cout << "--- READING NEW ONLINE DEMANDS... --- " << std::endl;
            std::string nextFile6 = instance6.getInput().getDemandToBeRoutedFilesFromIndex(0);
            instance6.generateDemandsFromFile(nextFile6);

            /********************************************************************/
            /* 		Solve - VOLUME LAGNONOVERLAP FORMULATION      	 		*/
            /********************************************************************/
            std::cout << "Solving with volume non overlap formualtion" << std::endl;
            ClockTime LAG4_OPTIMIZATION_TIME(ClockTime::getTimeNow());
            lagSolverFactory lagfactory4;
            AbstractLagSolver *lagsolver4 = lagfactory4.createSolver(instance6);
            lagsolver4->run();
            fichier << ";;" << lagsolver4->getUB() << ";" << lagsolver4->getLB() << ";" << lagsolver4->getIteration() << ";" << lagsolver4->getStop()<< ";" << std::fixed  << LAG4_OPTIMIZATION_TIME.getTimeInSecFromStart() << std::setprecision(9);
            std::cout << "Time taken by optimization is : ";
            std::cout << std::fixed  << LAG4_OPTIMIZATION_TIME.getTimeInSecFromStart() << std::setprecision(9); 
            std::cout << " sec" << std::endl; 

            fichier << std::endl;
        }
        fichier.close();
    }
	return 0;
}

