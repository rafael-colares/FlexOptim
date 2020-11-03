#ifndef LAG_SUBGRADIENT_H
#define LAG_SUBGRADIENT_H

#include "lagNonOverlapping.h"
#include "../topology/instance.h"

class lagSubgradient{

    private:
        const int MAX_NB_IT_WITHOUT_IMPROVEMENT;    /**< Maximum number of iterations without lower bound improvement. **/
        const int MAX_NB_IT;                        /**< Maximum number of performed iterations. **/
        const double INITIAL_STEPSIZE;
        const double MIN_STEPSIZE;

        
        lagNonOverlapping formulation;             /** Chosen Lagrangian formulation **/

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


    public:

        /************************************************/
	    /*				    Constructors 		   		*/
	    /************************************************/
        lagSubgradient(const Instance &inst);


        /************************************************/
        /*					   Getters 		    		*/
        /************************************************/
        int getIteration() const { return iteration; }
        int getItWithoutImprovement() const { return itWithoutImprovement; }

        double getUB() const { return UB; }
        double getLB() const { return LB; }

        double getStepSize() const { return stepSize; }
        double getLambda() const { return lambda; }

        /************************************************/
        /*					   Setters 		    		*/
        /************************************************/
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
        /*					   Methods 		    		*/
        /************************************************/

        /** Runs the subgradient method. **/
        void run();

        /** Sets all the initial parameters for the subgradient to run. **/
        void initialization();

        /** Sets the initial lambda used for updating the step size. **/
        void initLambda();

        /** Solves an iteration of the Subgradient Method. **/
        void runIteration();

        /* Updates the known lower bound. */
        void updateLB(double bound);

        /* Updates the known upper bound. */
        void updateUB(double bound);

        /* Updates the step size with the rule: lambda*(UB - Z[u])/|slack| */
        void updateStepSize();
        
        /* Updates the lambda used in the update of step size. Lambda is halved if LB has failed to increade in some fixed number of iterations */
        void updateLambda();

        bool isGradientMoving();

        void displayMainParameters();

};

#endif