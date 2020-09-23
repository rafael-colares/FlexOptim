#ifndef __abstractSolver__h
#define __abstractSolver__h

#include "../formulation/formulationFactory.h"


//typedef IloArray<IloNumVarArray> IloNumVarMatrix;
typedef IloArray<IloBoolVarArray> IloBoolVarMatrix;


/**********************************************************************************************
 * This is an abstract class modelling a MIP solver. Derived classes need to specify how the 
 * RSA formulation is solved, i.e., which MIP solved is used.
 * *******************************************************************************************/
class AbstractSolver{

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
	AbstractFormulation *formulation;
	Status currentStatus;

public:
	/****************************************************************************************/
	/*										Constructor										*/
	/****************************************************************************************/

	/** Constructor. The RSA constructor is called and the arc map storing the index of the preprocessed graphs associated is built. @param inst The instance to be solved. **/
    AbstractSolver(const Instance &instance, const Status &s = STATUS_UNKNOWN);

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

	
	/****************************************************************************************/
	/*											Methods										*/
	/****************************************************************************************/
	virtual void solve() = 0;
	
	virtual void implementFormulation() = 0;

	virtual void updateRSA(Instance &instance) = 0;

};    
#endif