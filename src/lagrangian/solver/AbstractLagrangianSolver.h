#ifndef ABSTRACT_LAG_SOLVER_H
#define ABSTRACT_LAG_SOLVER_H

#include "../formulation/lagFormulationFactory.h"
#include "../heuristicFactory.h"

/**********************************************************************************************
 * This is an abstract class modelling a Lagrangian solver. Derived classes need to specify how the 
 * RSA formulation is solved, i.e., which lagrangian formulation is used.
 * *******************************************************************************************/

class AbstractLagSolver{
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
        
        const int MAX_NB_IT_WITHOUT_IMPROVEMENT;    /**< Maximum number of iterations without lower bound improvement. **/
        const int MAX_NB_IT;                        /**< Maximum number of performed iterations. **/
        const double INITIAL_STEPSIZE;
        const double MIN_STEPSIZE;

        Status currentStatus;
        AbstractLagFormulation *formulation;
        AbstractHeuristic *heuristic;

        int iteration;
        int itWithoutImprovement;

        double UB;
        double LB;
        double currentLagrCost;
        double currentRealCost;

        /** Stores the value of the step size used for updating the Lagrangian multipliers. **/
        double stepSize;

        /** Stores the value of lambda used for updating the step size. **/
        double lambda;

        std::string stop;

        std::ofstream fichier;
        std::ofstream fichier2;

    public:
        /****************************************************************************************/
        /*										Constructor										*/
        /****************************************************************************************/

        /** Constructor. The RSA constructor is called and the arc map storing the index of the preprocessed graphs associated is built. @param inst The instance to be solved. **/
        AbstractLagSolver(const Instance &instance, const Status &s = STATUS_UNKNOWN);

        /****************************************************************************************/
        /*											Getters										*/
        /****************************************************************************************/

        /** Returns the status. **/
        Status getStatus() const {return currentStatus;}

        std::string getStop() const { return stop;}

        int getIteration() const { return iteration; }
        int getItWithoutImprovement() const { return itWithoutImprovement; }

        double getUB() const { return UB; }
        double getLB() const { return LB; }

        double getStepSize() const { return stepSize; }
        double getLambda() const { return lambda; }

        /************************************************/
        /*					   Setters 		    		*/
        /************************************************/

        /** Changes the status. **/
        void setStatus(const Status &s){ currentStatus = s; }

        void setStop(std::string s){ stop = s;}

        void setIteration(int i) { iteration = i; }
        void setItWithoutImprovement(int i) { itWithoutImprovement = i; }
        void incIteration() { iteration++; }
        void incItWithoutImprovement() { itWithoutImprovement++; }

        void setUB(double val){ UB = val; }
        void setLB(double val){ LB = val; }
       
        /** Changes the value of lambda. @param val The new value of lambda. **/
        void setLambda(double val){ lambda = val; }

        void setStepSize(double val){ stepSize = val; }

        /************************************************/
        /*					   Update 		    		*/
        /************************************************/

        /* Updates the known lower bound. */
        virtual void updateLB(double bound) = 0;

        /* Updates the known upper bound. */
        void updateUB(double bound);

        /* Updates the step size with the rule: lambda*(UB - Z[u])/|slack| */
        virtual void updateStepSize() = 0;
        
        /* Updates the lambda used in the update of step size. Lambda is halved if LB has failed to increade in some fixed number of iterations */
        virtual void updateLambda() = 0;

        /************************************************/
        /*					   Methods 		    		*/
        /************************************************/

        /** Sets the initial lambda used for updating the step size. **/
        void initLambda();

        bool isGradientMoving();

        /** Sets all the initial parameters for the subgradient to run. **/
        virtual void initialization() = 0;

        /** Runs the subgradient method. **/
        virtual void run() = 0;

        /** Solves an iteration of the Subgradient Method. **/
        virtual void runIteration() = 0;

        /************************************************/
        /*					   Display 		    		*/
        /************************************************/

        void displayMainParameters();

        /************************************************/
        /*				   Destructors 		    		*/
        /************************************************/

        virtual ~AbstractLagSolver();


};

#endif