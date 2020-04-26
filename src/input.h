#ifndef __input__h
#define __input__h

#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string>
#include <vector>

/*****************************************************************************************
 * This class contains all the information needed for the creation of an instance.
 * It stores input/output file paths, execution parameters (such as the chosen method 
 * used to solve the problem), and control parameters (such as max number of iterations).						
*****************************************************************************************/
class Input {

public: 
	/** Enumerates the possible methods to be applied for solving the Online Routing and Spectrum Allocation problem.**/
	enum Method {						
		METHOD_CPLEX = 0,  /**< Solve it through a MIP using CPLEX. **/
		METHOD_SUBGRADIENT = 1 /**<  Solve it using the subgradient method.**/
	};

	/** Enumerates the possible levels of applying a preprocessing step fo reducing the graphs before optimization is called. **/
	enum PreprocessingLevel {
		PREPROCESSING_LVL_NO = 0,		/**< Only remove arcs that do not fit the demand load. **/
		PREPROCESSING_LVL_PARTIAL = 1,	/**< Previous steps + look for arcs that would induce length violation and arcs whose neighboors cannot forward the demand. **/
		PREPROCESSING_LVL_FULL = 2		/**< Previous steps recursively until no additional arc can be removed. **/
	};

	/** Enumerates the possible objectives to be optimized. of applying a preprocessing step fo reducing the graphs before optimization is called. **/
	enum ObjectiveMetric {
		OBJECTIVE_METRIC_1 = 0,		/**< Minimize the sum of (max used slice positions) over demands. **/
		OBJECTIVE_METRIC_1p = 1,	/**< Minimize the sum of (max used slice positions) over edges. **/
		OBJECTIVE_METRIC_2 = 2,		/**< Minimize the sum of (number of hops in paths) over demands. **/
		OBJECTIVE_METRIC_4 = 4,		/**< Minimize the path lengths. **/
		OBJECTIVE_METRIC_8 = 8		/**< Minimize the max used slice position overall. **/
	};
	
private:
	const std::string PARAMETER_FILE;	/**< Path to the file containing all the parameters. **/
	std::string linkFile;				/**< Path to the file containing information on the physical topology of the network.**/
	std::string demandFile;				/**< Path to the file containing information on the already routed demands. **/
	std::string assignmentFile;			/**< Path to the file containing information on the assignment of demands (i.e., on which edge/slice each demand is routed).**/
	std::string onlineDemandFolder;		/**< Path to the folder containing the files on the non-routed demands. **/
	std::vector< std::string > vecOnlineDemandFile;	/**< A vector storing the paths to the files containing information on the non-routed demands. **/
	std::string outputPath;				/**< Path to the folder where the output files will be sent by the end of the optimization procedure.**/
	
    int nbDemandsAtOnce;				/**< How many demands are treated in a single optimization.**/
	int nbSlicesInOutputFile;			/**< How many slices will be displayed in the output file. **/
	Method chosenMethod;				/**< Refers to which method is applied for solving the problem.**/
	PreprocessingLevel chosenPreprLvl;	/**< Refers to which level of preprocessing is applied before solving the problem.**/
	ObjectiveMetric chosenObj;			/**< Refers to which objective is optimized.**/

	double lagrangianMultiplier_zero;	/**< The initial value of the lagrangian multiplier used if subgradient method is chosen. **/
	double lagrangianLambda_zero;		/**< The initial value of the lambda used for computing the step size if subgradient method is chosen. **/
	int nbIterationsWithoutImprovement;	/**< The maximal number of iterations allowed in the subgradient method.**/
	int maxNbIterations;				/**< The maximal number of iterations allowed without improving the lower bound in the subgradient method.**/
public:
	/****************************************************************************************/
	/*									Constructors										*/
	/****************************************************************************************/
	/** Default constructor initializes the object with the information contained in the parameterFile. @param file The address of the parameter file (usually the address of file 'onlineParameters'). **/
    Input(std::string file);

	/** Copy constructor. @param i The input to be copied. **/
    Input(const Input &i);

	/****************************************************************************************/
	/*										Getters											*/
	/****************************************************************************************/
	/** Returns the path to the file containing all the parameters. **/
    const std::string getParameterFile() const { return PARAMETER_FILE; }

	/** Returns the path to the file containing information on the physical topology of the network.**/
    std::string getLinkFile() const { return linkFile; }
	
	/** Returns the path to the file containing information on the already routed demands. **/
    std::string getDemandFile() const { return demandFile; }

	/** Returns the path to the file containing information on the assignment of demands (i.e., on which edge/slice each demand is routed).**/
    std::string getAssignmentFile() const { return assignmentFile; }
	
	/** Returns the path to the folder containing the files on the non-routed demands. **/
    std::string getOnlineDemandFolder() const { return onlineDemandFolder; }

	/** Returns the number of online demand files to be treated. **/
	int getNbOnlineDemandFiles(){ return vecOnlineDemandFile.size(); }
	
	/** Returns the vector storing the paths to the files containing information on the non-routed demands. **/
    std::vector<std::string> getOnlineDemandFiles() const { return vecOnlineDemandFile; }

	/** Returns the path to the i-th file containing information on the non-routed demands. @param i The index of file. **/
    std::string getOnlineDemandFilesFromIndex(int i) const { return vecOnlineDemandFile[i]; }

	/** Returns the path to the folder where the output files will be sent by the end of the optimization procedure.**/
    std::string getOutputPath() const { return outputPath; }
	
	/** Returns the number of demands to be treated in a single optimization. **/
    int getNbDemandsAtOnce() const {return nbDemandsAtOnce;}

	/** Returns the number of slices to be displayed in the output file. **/
    int getnbSlicesInOutputFile() const {return nbSlicesInOutputFile;}

	/** Returns the identifier of the method chosen for optimization. **/
    Method getChosenMethod() const {return chosenMethod;}

	/** Returns the identifier of the method chosen for optimization. **/
    PreprocessingLevel getChosenPreprLvl() const {return chosenPreprLvl;}

	/** Returns the identifier of the objective chosen to be optimized. **/
    ObjectiveMetric getChosenObj() const {return chosenObj;}

	/** Returns the initial value of the lagrangian multiplier used if subgradient method is chosen. **/
	double getInitialLagrangianMultiplier() const { return lagrangianMultiplier_zero; }
	
	/** Returns the initial value of the lambda used for computing the step size if subgradient method is chosen. **/
	double getInitialLagrangianLambda() const { return lagrangianLambda_zero; }
	
	/** Returns the maximal number of iterations allowed in the subgradient method.**/
	int getNbIterationsWithoutImprovement() const { return nbIterationsWithoutImprovement; }

	/** Returns the maximal number of iterations allowed without improving the lower bound in the subgradient method.**/
	int getMaxNbIterations() const { return maxNbIterations; }

	/****************************************************************************************/
	/*										Methods											*/
	/****************************************************************************************/
	/** Searches for a pattern in the parameter file and returns its associated value. @param pattern The substring pattern to be searched. @note It is surely not the most performant method but the parameter file has a reasonable size and the code becomes much clearer. **/
    std::string getParameterValue(std::string pattern);

	/** Populates the vector of paths to the files containing information on the non-routed demands. **/
	void populateOnlineDemandFiles();

	/** Converts a string into an ObjectiveMetric. **/
	ObjectiveMetric to_ObjectiveMetric(std::string data);

	/** Displays the main input file paths: link, demand and assignement. **/
    void displayMainParameters();
};

#endif