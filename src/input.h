#ifndef __input__h
#define __input__h

#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string>
#include <vector>
#include <climits>

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
		METHOD_SUBGRADIENT = 1, /**<  Solve it using the subgradient method.**/
		METHOD_YOUSSOUF = 2  /**< Solve it through a MIP using CPLEX. **/
	};

	/** Enumerates the possible levels of applying a preprocessing step fo reducing the graphs before optimization is called. **/
	enum PreprocessingLevel {
		PREPROCESSING_LVL_NO = 0,		/**< Only remove arcs that do not fit the demand load. **/
		PREPROCESSING_LVL_PARTIAL = 1,	/**< Previous steps + look for arcs that would induce length violation and arcs whose neighboors cannot forward the demand. **/
		PREPROCESSING_LVL_FULL = 2		/**< Previous steps recursively until no additional arc can be removed. **/
	};

	/** Enumerates the possible output policies to be used. **/
	enum OutputLevel {
		OUTPUT_LVL_NO = 0,			/**< Do not create any output file. **/
		OUTPUT_LVL_NORMAL = 1,		/**< Generate output files corresponding to the last mapping before blocking. **/
		OUTPUT_LVL_DETAILED = 2		/**< Generate output files after every optimization procedure. **/
	};

	/** Enumerates the possible objectives to be optimized. **/
	enum ObjectiveMetric {
		OBJECTIVE_METRIC_0 = 0,		/**< Minimize nothing, just search for a feasible solution. **/
		OBJECTIVE_METRIC_1 = 1,		/**< Minimize the sum of (max used slice positions) over demands. **/
		OBJECTIVE_METRIC_1p = 11,	/**< Minimize the sum of (max used slice positions) over edges. **/
		OBJECTIVE_METRIC_2 = 2,		/**< Minimize the sum of (number of hops in paths) over demands. **/
		OBJECTIVE_METRIC_4 = 4,		/**< Minimize the path lengths. **/
		OBJECTIVE_METRIC_8 = 8		/**< Minimize the max used slice position overall. **/
	};
	
	/** Enumerates the possible spectrum partitioning policies to be applied. **/
	enum PartitionPolicy {
		PARTITION_POLICY_NO = 0,	/**< No spectrum partition. **/
		PARTITION_POLICY_SOFT = 1,	/**< A subset of demands is pushed to the Right region of the spectrum. **/
		PARTITION_POLICY_HARD = 2	/**< Divide the spectrum into two fixed regions (Left and Right). **/
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
	int partitionSlice;					/**< Refers to the max slice of Left spectrum region, if partitioning policy is set to Hard. **/
	int partitionLoad;					/**< Refers to the max load that can be routed on the Left spectrum region, if some partioning policy is set. **/
	int timeLimit;						/**< Refers to how much time (in seconds) can be spent during one optimization. **/
	int globalTimeLimit;				/**< Refers to how much time (in seconds) can be spent during the whole optmization. **/
	bool allowBlocking;					/**< If this option is inactive, optimization stops within first blocking. Otherwise, blocking is accepted (this only works in online case). **/

	Method chosenMethod;				/**< Refers to which method is applied for solving the problem.**/
	PreprocessingLevel chosenPreprLvl;	/**< Refers to which level of preprocessing is applied before solving the problem.**/
	ObjectiveMetric chosenObj;			/**< Refers to which objective is optimized.**/
	OutputLevel chosenOutputLvl;		/**< Refers to which output policy is adopted.**/
	PartitionPolicy chosenPartitionPolicy;	/**< Refers to which partition policy is adopted.**/

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
    int getNbDemandsAtOnce() const { return nbDemandsAtOnce; }

	/** Returns the number of slices to be displayed in the output file. **/
    int getnbSlicesInOutputFile() const { return nbSlicesInOutputFile; }

	/** Returns the max slice of Left spectrum region. **/
    int getPartitionSlice() const { return partitionSlice; }

	/** Returns the max load that should be pushed to the Left spectrum region. **/
    int getPartitionLoad() const { return partitionLoad; }

	/** Returns the timeLimit for a single optimization iteration. **/
    int getIterationTimeLimit() const { return timeLimit; }

	/** Returns the global time limit applied to the optimization. **/
    int getOptimizationTimeLimit() const { return globalTimeLimit; }

	/** Returns true if blocking is accepted and false, otherwise. **/
    bool isBlockingAllowed() const { return allowBlocking; }

	/** Returns the identifier of the method chosen for optimization. **/
    Method getChosenMethod() const { return chosenMethod; }

	/** Returns the identifier of the method chosen for optimization. **/
    PreprocessingLevel getChosenPreprLvl() const { return chosenPreprLvl; }

	/** Returns the identifier of the objective chosen to be optimized. **/
    ObjectiveMetric getChosenObj() const { return chosenObj; }

	/** Returns the identifier of the output policy adopted. **/
    OutputLevel getChosenOutputLvl() const { return chosenOutputLvl; }

	/** Returns the identifier of the partition policy adopted. **/
    PartitionPolicy getChosenPartitionPolicy() const { return chosenPartitionPolicy; }

	/** Returns the initial value of the lagrangian multiplier used if subgradient method is chosen. **/
	double getInitialLagrangianMultiplier() const { return lagrangianMultiplier_zero; }
	
	/** Returns the initial value of the lambda used for computing the step size if subgradient method is chosen. **/
	double getInitialLagrangianLambda() const { return lagrangianLambda_zero; }
	
	/** Returns the maximal number of iterations allowed in the subgradient method.**/
	int getNbIterationsWithoutImprovement() const { return nbIterationsWithoutImprovement; }

	/** Returns the maximal number of iterations allowed without improving the lower bound in the subgradient method.**/
	int getMaxNbIterations() const { return maxNbIterations; }

	/****************************************************************************************/
	/*										Setters											*/
	/****************************************************************************************/
	
	/** Changes the number of demands to be treated in a single optimization. @param val The new number of demands. **/
    void setNbDemandsAtOnce(const int val) { nbDemandsAtOnce = val; }

	/** Changes the time limit of one iteration. @param val The new time limit (in seconds). **/
    void setTimeLimit(const int val) { timeLimit = val; }

	/** Changes the global time limit. @param val The new time limit (in seconds). **/
    void setGlobalTimeLimit(const int val) { globalTimeLimit = val; }

	/****************************************************************************************/
	/*										Methods											*/
	/****************************************************************************************/
	/** Searches for a pattern in the parameter file and returns its associated value. @param pattern The substring pattern to be searched. @note It is surely not the most performant method but the parameter file has a reasonable size and the code becomes much clearer. **/
    std::string getParameterValue(std::string pattern);

	/** Populates the vector of paths to the files containing information on the non-routed demands. **/
	void populateOnlineDemandFiles();

	/** Converts a string into an ObjectiveMetric. **/
	ObjectiveMetric to_ObjectiveMetric(std::string data);

	/** Converts a string into a PartitionPolicy. **/
	PartitionPolicy to_PartitionPolicy(std::string data);

	/** Converts a string into time limit. \note By default, time limit is unlimited. **/
	int to_timeLimit(std::string data);

	/** Displays the main input file paths: link, demand and assignement. **/
    void displayMainParameters();
	
	/****************************************************************************************/
	/*										Destructor										*/
	/****************************************************************************************/

	/** Destructor. Clears the vector of strings. **/
	~Input();
};

#endif