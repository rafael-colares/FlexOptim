#ifndef LAG_SUBGRADIENT_H
#define LAG_SUBGRADIENT_H

#include "AbstractLagrangianSolver.h"
#include "../../topology/instance.h"

class lagSubgradient: public AbstractLagSolver{
    private:
        int globalItWithoutImprovement;

    public:

        /************************************************/
	    /*				    Constructors 		   		*/
	    /************************************************/
        lagSubgradient(const Instance &inst):AbstractLagSolver(inst){}

        /************************************************/
        /*					   Getters 		    		*/
        /************************************************/
        int getGlobalItWithoutImprovement() const { return globalItWithoutImprovement;}

        /************************************************/
        /*					   Setters 		    		*/
        /************************************************/
        void setGlobalItWithoutImprovement(int value) {globalItWithoutImprovement = value;}
        void incGlobalItWithoutImprovement() {globalItWithoutImprovement++;}

        /************************************************/
        /*					   Update 		    		*/
        /************************************************/
        
        /* Updates the lambda used in the update of step size. Lambda is halved if LB has failed to increade in some fixed number of iterations */
        void updateLambda();

        void updateLB(double);

        void updateStepSize();

        /************************************************/
        /*					   Methods 		    		*/
        /************************************************/

        /** Sets all the initial parameters for the subgradient to run. **/
        void initialization();

        /** Runs the subgradient method. **/
        void run();

        /** Solves an iteration of the Subgradient Method. **/
        void runIteration();

        ~lagSubgradient(){}

};

#endif