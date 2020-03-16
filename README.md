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
Obs.: onlineParameters is the file where one can input the necessary parameters for the optimization such as the location of .csv files needed for the building of an initial mapping.


