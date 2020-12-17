#ifndef LAG_VOLUME_H
#define LAG_VOLUME_H

#include "AbstractLagrangianSolver.h"

class lagVolume: public AbstractLagSolver{

    private:

        const int MAX_NB_IT_RED;                    /**< Maximum number of red iterations. **/
        const int MAX_NB_IT_YELLOW;                 /**< Maximum number of yellow iterations. **/
        const double UPD_STEPSIZE_GREEN_YELLOW;     /**< Update stepsize when green or consecutive yellow iterations. **/
        const double UPD_STEPSIZE_RED;              /**< Update stepsize when consecutive red iterations. **/
        const double MIN_MODULE_VALUE;              /**< If module smaller than this, stop. **/
        const double MIN_DIF_OBJ_VALUE;             /**< If difference between lower bound and primal cost smaller than this, stop. **/
        double MAX_ALPHA;                           /**< Initial maximum for alpha. **/

        bool greenIt;
        int nbRedIt;
        int nbYellowIt;

        double previous_lb;

    public:
        /************************************************/
	    /*				    Constructors 		   		*/
	    /************************************************/
        lagVolume(const Instance &inst):AbstractLagSolver(inst),MAX_NB_IT_RED(20),MAX_NB_IT_YELLOW(2),UPD_STEPSIZE_GREEN_YELLOW(1.1),UPD_STEPSIZE_RED(0.67),MIN_MODULE_VALUE(0.01),MIN_DIF_OBJ_VALUE(0.02),MAX_ALPHA(0.1){}

        /************************************************/
	    /*				    GETTERS      		   		*/
	    /************************************************/

        bool getGreenIt() const { return greenIt;}
        int getNbRedIt() const { return nbRedIt;}
        int getNbYellowIt() const {return nbYellowIt;}

        /************************************************/
	    /*				    SETTERS      		   		*/
	    /************************************************/

        void setGreenIt(bool b) {greenIt = b;}
        void setNbRedIt(int nb) {nbRedIt = nb;}
        void incNbRedIt() {nbRedIt++;}
        void setNbYellowIt(int nb) {nbYellowIt = nb;}
        void incNbYellowIt() {nbYellowIt++;}

        void initialization();

        void run();

        void runIteration();

        void updateLB(double);

        void updateParamAlpha();

        void updateLambda();

        void updateStepSize();

        double getAlphaValue();

};

#endif