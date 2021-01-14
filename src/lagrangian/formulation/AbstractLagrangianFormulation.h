#ifndef ABSTRACT_LAG_FORMULATION_H
#define ABSTRACT_LAG_FORMULATION_H

#include "../../formulation/rsa.h"

class AbstractLagFormulation: public RSA{

        protected:

                double currentLagrCost;
                double currentRealCost;

                std::vector< std::vector<double> > primal_linear_solution;

                /***************************** ASSIGNMENT MATRIX *******************************/

                std::vector< std::vector<bool> > assignmentMatrix_d;
        public:

                /***********************************************************************************/
                /*			            CONSTRUCTORS		                   */
                /***********************************************************************************/
                /** Constructor. Builds the Formulation.  @param instance The instance to be solved. **/
                AbstractLagFormulation(const Instance &instance);

                /* ******************************************************************************* /
                *                             INITIALIZATION METHODS                              
                ********************************************************************************* */
                
                /** Sets all initial parameters **/
                virtual void init() = 0;

                void initPrimalSolution();

                /* *******************************************************************************
                *                             RUNING METHODS
                ******************************************************************************* */
                
                /** Solves the RSA using the Subgradient Method. **/
                virtual void run() = 0;

                virtual bool checkFeasibility() = 0;

                /***********************************************************************************/
                /*		                       GETTERS 	     	    	                   */
                /***********************************************************************************/
                double getLagrCurrentCost() const { return currentLagrCost; }
                double getRealCurrentCost() const { return currentRealCost; }

                virtual double getSlackModule() = 0;

                virtual double getSlackModule_v2() = 0;

                virtual double get_prod_slack() = 0;

                virtual double getMeanSlackModule_v2() =0;

                virtual double getSlackDirectionProd() = 0;

                double getPrimalObjective();

                std::vector< std::vector<bool> > getVariables() const { return assignmentMatrix_d;}
                std::vector< std::vector<double> > getPrimalVariables() const { return primal_linear_solution;}

                std::shared_ptr<ListDigraph> getVecGraphD(int d) const { return vecGraph[d];}

                /***********************************************************************************/
                /*		                       SETTERS 		    	                   */
                /***********************************************************************************/
                void setCurrentLagrCost(double val){ currentLagrCost = val; }
                void incCurrentLagrCost(double val){ currentLagrCost += val; }
                void setCurrentRealCost(double val){ currentRealCost = val; }
                void incCurrentRealCost(double val){ currentRealCost += val; }

                /****************************************************************************************/
                /*				        UPDATE						*/
                /****************************************************************************************/

                /********************************* MULTIPLIERS ***********************************/
        
                /* Updates lagrangian multipliers with the rule: u[k+1] = u[k] + t[k]*violation */
                virtual void updateMultiplier(double) = 0;

                virtual void updateMultiplier_v2(double) = 0;

                virtual void updateStabilityCenter() = 0;

                /********************************* COSTS ***********************************/

                /** Updates the arc costs. @note cost = coeff + u_d*length **/
                virtual void updateCosts() = 0;

                /********************************* PRIMAL SOLUTION ***********************************/

                void updatePrimalSolution(double); 

                /****************************************************************************************/
                /*					Display						*/
                /****************************************************************************************/

                virtual void displaySlack(std::ostream & = std::cout) = 0;

                virtual void displayMultiplier(std::ostream & = std::cout) = 0;

                /* *******************************************************************************
                *                             DESTRUCTOR  
                ******************************************************************************* */
                virtual ~AbstractLagFormulation(){}

};

#endif


