#ifndef __input__h
#define __input__h

#include <iostream>
#include <fstream>
#include <string>

/*****************************************************************************************
 * This class contains all the information needed for the creation of an instance.
 * It stores input/output file paths, execution parameters (such as the chosen method 
 * used to solve the problem), and control parameters (such as max number of iterations).						
*****************************************************************************************/
class Input {

private:
	const std::string PARAMETER_FILE;	/**< Path to the file containing all the parameters. **/
	std::string linkFile;				/**< Path to the file containing information on the physical topology of the network.**/
	std::string demandFile;				/**< Path to the file containing information on the already routed demands. **/
	std::string assignmentFile;			/**< Path to the file containing information the assignment of demands (i.e., on which edge/slice each demand is routed).**/
	std::string onlineDemandFile;		/**< Path to the file containing information on the non-routed demands. **/
	std::string outputPath;				/**< Path to the folder where the output files will be sent by the end of the optimization procedure.**/
	
    int nbDemandsAtOnce;				/**< How many demands are treated in a single optimization.**/
	int nbSlicesInOutputFile;			/**< How many slices will be displayed in the output file. **/
	int chosenMethod;					/**< Refers to which method is applied for solving the online RSA. If 1, use CPLEX. If 2, use subgradient.**/

	double lagrangianMultiplier_zero;	/**< The initial value of the lagrangian multiplier used if subgradient method is chosen. **/
	double lagrangianLambda_zero;		/**< The initial value of the lambda used for computing the step size if subgradient method is chosen. **/
	int nbIterationsWithoutImprovement;	/**< The maximal number of iterations allowed in the subgradient method.**/
	int maxNbIterations;				/**< The maximal number of iterations allowed without improving the lower bound in the subgradient method.**/
public:
	/************************************************/
	/*				   Constructors					*/
	/************************************************/
	/** Default constructor initializes the object with the information contained in the parameterFile. **/
    Input(std::string parameterFilePath);
	/** Copy constructor. **/
    Input(const Input &i);

	/************************************************/
	/*				     Getters					*/
	/************************************************/
	/** Returns the path to the file containing all the parameters. **/
    const std::string getParameterFile() const { return PARAMETER_FILE; }

	/** Returns the path to the file containing information on the physical topology of the network.**/
    std::string getLinkFile() const { return linkFile; }
	
	/** Returns the path to the file containing information on the already routed demands. **/
    std::string getDemandFile() const { return demandFile; }

	/** Returns the path to the file containing information on the assignment of demands (i.e., on which edge/slice each demand is routed).**/
    std::string getAssignmentFile() const { return assignmentFile; }
	
	/** Returns the path to the file containing information on the non-routed demands. **/
    std::string getOnlineDemandFile() const { return onlineDemandFile; }

	/** Returns the path to the folder where the output files will be sent by the end of the optimization procedure.**/
    std::string getOutputPath() const { return outputPath; }
	
	/** Returns the number of demands to be treated in a single optimization. **/
    int getNbDemandsAtOnce() const {return nbDemandsAtOnce;}

	/** Returns the number of slices to be displayed in the output file. **/
    int getnbSlicesInOutputFile() const {return nbSlicesInOutputFile;}

	/** Returns the identifier of the method chosen for optimization. **/
    int getChosenMethod() const {return chosenMethod;}

	/** Returns the initial value of the lagrangian multiplier used if subgradient method is chosen. **/
	double getInitialLagrangianMultiplier() const { return lagrangianMultiplier_zero; }
	
	/** Returns the initial value of the lambda used for computing the step size if subgradient method is chosen. **/
	double getInitialLagrangianLambda() const { return lagrangianLambda_zero; }
	
	/** Returns the maximal number of iterations allowed in the subgradient method.**/
	int getNbIterationsWithoutImprovement() const { return nbIterationsWithoutImprovement; }

	/** Returns the maximal number of iterations allowed without improving the lower bound in the subgradient method.**/
	int getMaxNbIterations() const { return maxNbIterations; }

	/************************************************/
	/*				     Methods					*/
	/************************************************/
	/** Searches for a pattern in the Parameter File (PF) and returns its associated value. @warning It is surely not the most performant method but PF has a reasonable size and the code becomes much clearer. **/
    std::string getParameterValue(std::string pattern);

	/** Displays the main input file paths: link, demand and assignement. **/
    void displayMainParameters();
};

#endif