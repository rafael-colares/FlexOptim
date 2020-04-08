# FlexOptim
This is a repository for the online algorithms used for solving the RSA problem from FlexOptim project.


# Requirements
- This code was developped for linux. If you use Windows, please install Windows Subsystem Linux (WSL) by following the steps described here: https://docs.microsoft.com/en-us/windows/wsl/install-win10

- This code uses IBM ILOG CPLEX version 12.10 which can be installed by following the instructions here: https://developer.ibm.com/docloud/blog/2019/07/04/cplex-optimization-studio-for-students-and-academics/ 

- We also use the COIN-OR lemon library to model graphs. Lemon can be installed by followig the instructions described here:
https://lemon.cs.elte.hu/trac/lemon/wiki/InstallGuide

- Boost library is also used for simplifying the lecture from .csv files. Boost library can found here: https://www.boost.org/doc/libs/1_55_0/doc/html/bbv2/installation.html

# Compilation
The compilation of the project is done through a makefile. If you have all the requirements described above, just update the paths of CPLEX, lemon and boost in the makefile to the correspondent directories where you have installed each of these libraries.

Open a terminal in linux (or WSL), go to the project folder and type make. Some warnings are expected but harmless.

# Running
After compiling the project, one can execute it by typing ./exec onlineParameters

# Parameters
The file onlineParameters contains the parameters needed for the program to run properly.
- linkFile refers to the address of the file containing information on the physical topology links.
- demandFile refers to the address of the file containing information on the demands already routed in the initial mapping.
- assignmentFile refers to the address of the file containing information on the assignment of demands, that is, on which edge/slice each demand is routed.
- onlineDemandFile refers to the address of the file containing information on the demands that are to be routed through the online procedure.
- nbDemandsAtOnce states how many demands are treated in one optimization step.
- outputPath refers to the folder address where the output files will be sent by the end of the optimization procedure.
- nbSlicesInOutputFile states how many slices will be displayed in the output file
- chosenMethod refers to which method is applied for solving the online RSA. For using CPLEX, choose 1. For subgradient, choose 2.

Next parameters are reserved for team LIMOS.
- lagrangianMultiplier_zero refers to the initial value of the lagrangian multiplier used if subgradient method is chosen.
- lagrangianLambda_zero refers to the initial value of the lambda used for computing the step size if subgradient method is chosen.
- maxNbIterations states the maximal number of iterations subgradient method is allowed.
- nbIterationsWithoutImprovement states the maximal number of iterarions the subgradient method is allowed without improving the lower bound.
