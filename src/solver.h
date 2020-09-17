#ifndef __solver__h
#define __solver__h

#include "RSA.h"


//typedef IloArray<IloNumVarArray> IloNumVarMatrix;
typedef IloArray<IloBoolVarArray> IloBoolVarMatrix;


/**********************************************************************************************
 * This class specifies which MIP solver to be used The methods needed for solving 
 * the Routing and Spectrum Allocation problem via CPLEX. 
 * \note It uses the LEMON library to build the arc maps and CPLEX concert library to build 
 * expressions and constraints. 
 * *******************************************************************************************/
class Solver{

public:
	/** Enumerates the possible algorithm status according to the current model and solution in hand. **/
	enum Status {
        STATUS_UNKNOWN	= 0,        /**< The algorithm has no information about the solution of the model. **/					
		STATUS_FEASIBLE = 1,        /**< The algorithm found a feasible solution that may not necessarily be optimal. **/
		STATUS_OPTIMAL = 2,         /**< The algorithm found an optimal solution. **/
        STATUS_INFEASIBLE = 3,      /**< The algorithm proved the model infeasible; that is, it is not possible to find a feasible solution. **/
        STATUS_UNBOUNDED = 4,       /**< The algorithm proved the model unbounded. **/
        STATUS_INFEASIBLE_OR_UNBOUNDED = 5, /**< The model is infeasible or unbounded. **/
        STATUS_ERROR = 6            /**< An error occurred. **/
	};
protected:
	Status currentStatus;

public:
	/****************************************************************************************/
	/*										Constructor										*/
	/****************************************************************************************/

	/** Constructor. The RSA constructor is called and the arc map storing the index of the preprocessed graphs associated is built. @param inst The instance to be solved. **/
    Solver(const Status &s = STATUS_UNKNOWN);

	/****************************************************************************************/
	/*											Getters										*/
	/****************************************************************************************/
	/** Returns the status. **/
    virtual Status getStatus() = 0;

	/****************************************************************************************/
	/*											Setters										*/
	/****************************************************************************************/
	/** Changes the status. **/
    void setStatus(const Status &s){ currentStatus = s; }

};    
#endif