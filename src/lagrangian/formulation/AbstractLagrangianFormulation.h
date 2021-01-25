#ifndef ABSTRACT_LAG_FORMULATION_H
#define ABSTRACT_LAG_FORMULATION_H

#include "../../formulation/rsa.h"
#include "../../tools/clockTime.h"

#include <lemon/bellman_ford.h>

/** This class implements a general Lagrangian Formulation considering the Flow formulation **/

class AbstractLagFormulation: public RSA{

        protected:

                ClockTime time;
                double constAuxGraphTime;

                double currentLagrCost;
                double currentRealCost;

                /********************************* MULTIPLIERS ***********************************/

                /** A vector storing the value of the Lagrangian multipliers associated with Length Constraints. **/
                std::vector<double> lagrangianMultiplierLength;

                /** A 2-dimensional vector storing the value of the Lagrangian multipliers associated with each Source/Target constraint. **/
                std::vector<std::vector<double>> lagrangianMultiplierSourceTarget;

                /** A 2-dimensional vector storing the value of the Lagrangian multipliers associated with each Flow constraint. **/
                std::vector<std::vector<double>> lagrangianMultiplierFlow;

                /** A 2-dimensional vector storing the value of the Lagrangian multipliers associated with each Non-Overlapping constraint. **/
                std::vector< std::vector<double> > lagrangianMultiplierOverlap;

                /** A 1-dimensional vector storing the value of the Lagrangian multipliers associated with Max used slice overall constrant. **/
                std::vector<double > lagrangianMultiplierMaxUsedSliceOverall;

                /** A 2-dimensional vector storing the value of the Lagrangian multipliers associated with Max used slice overall 2 constrant. **/
                std::vector< std::vector<double> > lagrangianMultiplierMaxUsedSliceOverall2;

                /** A 2-dimensional vector storing the value of the Lagrangian multipliers associated with Max used slice overall 3 constrant. **/
                std::vector< std::vector<double> > lagrangianMultiplierMaxUsedSliceOverall3;

                /********************************* STABILITY CENTER ***********************************/

                /** A vector storing the value of the Lagrangian stability center associated with Length Constraints. **/
                std::vector<double> lagrangianSCLength;

                std::vector<std::vector<double>> lagrangianSCSourceTarget;

                std::vector<std::vector<double>> lagrangianSCFlow;

                /** A 2-dimensional vector storing the value of the Lagrangian stability center associated with each Non-Overlapping constraint. **/
                std::vector< std::vector<double> > lagrangianSCOverlap;

                /** A 1-dimensional vector storing the value of the Lagrangian stability center associated with Max used slice overall constrant. **/
                std::vector<double > lagrangianSCMaxUsedSliceOverall;

                /** A 2-dimensional vector storing the value of the Lagrangian stability center associated with Max used slice overall 2 constrant. **/
                std::vector< std::vector<double> > lagrangianSCMaxUsedSliceOverall2;

                /** A 2-dimensional vector storing the value of the Lagrangian stability center associated with Max used slice overall 3 constrant. **/
                std::vector< std::vector<double> > lagrangianSCMaxUsedSliceOverall3;

                /********************************** SLACK ***************************************/

                /** Stores the value of the slack of lengths constraints (i.e., b - Dx). */
                std::vector<double> lengthSlack;

                /** Stores the value of the slack of Source/Target constraints (i.e., b - Dx). */
                std::vector<std::vector<double>> sourceTargetSlack;

                /** Stores the value of the slack of flow constraints (i.e., b - Dx). */
                std::vector<std::vector<double>> flowSlack;

                /** Stores the value of the slack of Non-Overlapping constraints (i.e., b - Dx). */
                std::vector< std::vector<double> > overlapSlack;

                /** Stores the value of the slack of Max used slice overall constraints (i.e., b - Dx). */
                std::vector<double> maxUsedSliceOverallSlack;

                /** Stores the value of the slack of Max used slice overall 2 constraints (i.e., b - Dx). */
                std::vector< std::vector<double> > maxUsedSliceOverallSlack2;

                /** Stores the value of the slack of Max used slice overall constraints (i.e., b - Dx). */
                std::vector< std::vector<double> > maxUsedSliceOverallSlack3;

                /************************** SLACK CONSIDERING PRIMAL VECTOR ***********************/
                std::vector<double> lengthSlack_v2;

                std::vector<std::vector<double>> sourceTargetSlack_v2;

                std::vector<std::vector<double>> flowSlack_v2;

                std::vector< std::vector<double> > overlapSlack_v2;

                /** Stores the value of the slack (primal variables) of Max used slice overall constraints (i.e., b - Dx). */
                std::vector<double> maxUsedSliceOverallSlack_v2;

                /** Stores the value of the slack (primal variables) of Max used slice overall 2 constraints (i.e., b - Dx). */
                std::vector< std::vector<double> > maxUsedSliceOverallSlack2_v2;

                /** Stores the value of the slack (primal variables) of Max used slice overall constraints (i.e., b - Dx). */
                std::vector< std::vector<double> > maxUsedSliceOverallSlack3_v2;

                /******************************** DIRECTION ***************************************/

                /** Stores the value of the direction of lengths constraints (i.e., b - Dx). */
                std::vector<double> lengthDirection;

                /** Stores the value of the direction of Source/Target constraints (i.e., b - Dx). */
                std::vector<std::vector<double>> sourceTargetDirection;

                /** Stores the value of the direction of flow constraints (i.e., b - Dx). */
                std::vector<std::vector<double>> flowDirection;

                /** Stores the value of the direction of overlap constraints (i.e., b - Dx). */
                std::vector< std::vector<double> > overlapDirection;

                /** Stores the value of the direction of Max used slice overall constraints (i.e., b - Dx). */
                std::vector<double> maxUsedSliceOverallDirection;

                /** Stores the value of the direction of Max used slice overall 2 constraints (i.e., b - Dx). */
                std::vector< std::vector<double> > maxUsedSliceOverall2Direction;

                /** Stores the value of the direction of Max used slice overall constraints (i.e., b - Dx). */
                std::vector< std::vector<double> > maxUsedSliceOverall3Direction;


                /***************************** PRIMAL APPROXIMATION *****************************/

                std::vector< std::vector<double> > primal_linear_solution;

                /***************************** ASSIGNMENT MATRIX *******************************/

                std::vector< std::vector<bool> > assignmentMatrix_d;

                double maxUsedSliceOverall;

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

                /********************************* MULTIPLIERS ***********************************/

                /** Sets the initial lagrangian multipliers values for the subgradient to run. **/
                virtual void initMultipliers() = 0;

                /** Sets the initial lagrangian multipliers associated with length constraints. **/
                void initializeLengthMultipliers(double);

                /** Sets the initial lagrangian multipliers associated with Source/Target constraints **/
                void initializeSourceTargetMultipliers(double);
                
                /** Sets the initial lagrangian multipliers associated with flow constraints **/
                void initializeFlowMultipliers(double);

                /** Sets the initial lagrangian multipliers associated with non-overlapping constraints. **/
                void initializeOverlapMultipliers(double);

                /** Sets the initial lagrangian multipliers associated with max used slice overall constraints. **/
                void initializeMaxUsedSliceOverallMultipliers(double);

                /** Sets the initial lagrangian multipliers associated with max used slice overall 2 constraints. **/
                void initializeMaxUsedSliceOverall2Multipliers(double);

                /** Sets the initial lagrangian multipliers associated with max used slice overall 3 constraints. **/
                void initializeMaxUsedSliceOverall3Multipliers(double);

                /********************************* STABILITY CENTER ***********************************/

                /** Initializes the slack of relaxed constraints. **/
                virtual void initStabilityCenter() = 0;

                /** Sets the initial lagrangian stability center associated with length constraints. **/
                void initializeLengthSC();

                /** Sets the initial lagrangian stability center associated with Source/Target constraints **/
                void initializeSourceTargetSC();
                
                /** Sets the initial lagrangian stability center associated with flow constraints **/
                void initializeFlowSC();

                /** Sets the initial lagrangian stability center associated with overlap constraints **/
                void initializeOverlapSC();

                /** Sets the initial lagrangian stability center associated with max used slice overall constraints **/
                void initializeMaxUsedSliceOverallSC();

                /** Sets the initial lagrangian stability center associated with max used slice overall 2 constraints **/
                void initializeMaxUsedSliceOverall2SC();

                /** Sets the initial lagrangian stability center associated with max used slice overall 3 constraints **/
                void initializeMaxUsedSliceOverall3SC();

                /********************************** SLACK ***************************************/

                /** Initializes the slack of relaxed constraints. **/
                virtual void initSlacks() = 0;

                /** Initializes the slack of length constraints. **/
                void initializeLengthSlacks();
                
                /** Initializes the slack of Source/Target constraints. **/
                void initializeSourceTargetSlacks();

                /** Initializes the slack of Flow constraints. **/
                void initializeFlowSlacks();

                /** Initializes the slack of non-overlap constraints. **/
                void initializeOverlapSlacks();

                /** Sets the initial lagrangian slack associated with max used slice overall constraints **/
                void initializeMaxUsedSliceOverallSlacks();

                /** Sets the initial lagrangian slack associated with max used slice overall 2 constraints **/
                void initializeMaxUsedSliceOverall2Slacks();

                /** Sets the initial lagrangian slack associated with max used slice overall 3 constraints **/
                void initializeMaxUsedSliceOverall3Slacks();

                /********************************** SLACK PRIMAL VARIABLES **********************************/

                /** Initialize slack considering primal variables **/
                virtual void initSlacks_v2() = 0;

                /** Initializes the slack of length constraints. **/
                void initializeLengthSlacks_v2();
                
                /** Initializes the slack of Source/Target constraints. **/
                void initializeSourceTargetSlacks_v2();

                /** Initializes the slack of Flow constraints. **/
                void initializeFlowSlacks_v2();

                /** Initializes the slack of non-overlap constraints. **/
                void initializeOverlapSlacks_v2();

                /** Sets the initial lagrangian slack (primal variables) associated with max used slice overall constraints **/
                void initializeMaxUsedSliceOverallSlacks_v2();

                /** Sets the initial lagrangian slack (primal variables) associated with max used slice overall 2 constraints **/
                void initializeMaxUsedSliceOverall2Slacks_v2();

                /** Sets the initial lagrangian slack (primal variables) associated with max used slice overall 3 constraints **/
                void initializeMaxUsedSliceOverall3Slacks_v2();

                /********************************** DIRECTION ***************************************/

                virtual void initDirection() = 0;

                void initializeLengthDirection();

                void initializeSourceTargetDirection();

                void initializeFlowDirection();

                void initializeOverlapDirection();

                void initializeMaxUsedSliceOverallDirection();

                void initializeMaxUsedSliceOverall2Direction();

                void initializeMaxUsedSliceOverall3Direction();

                /********************************** ASSIGNMENT MATRIX ************************************/

                virtual void initAssignmentMatrix();

                /********************************** PRIMAL SOLUTION ************************************/

                void initPrimalSolution();

                /* *******************************************************************************
                *                             RUNING METHODS
                ******************************************************************************* */
                
                /** Solves the RSA using the Subgradient Method. **/
                virtual void run() = 0;

                virtual void subtractConstantValuesFromLagrCost() = 0;

                virtual bool checkFeasibility() = 0;

                bool checkLengthFeasibility();

                bool checkSourceTargetFeasibility();

                bool checkFlowFeasibility();

                bool checkOverlapFeasibility();


                /***********************************************************************************/
                /*		                       GETTERS 	     	    	                   */
                /***********************************************************************************/
                double getLagrCurrentCost() const { return currentLagrCost; }
                double getRealCurrentCost() const { return currentRealCost; }

                double getConstAuxGraphTime() const { return constAuxGraphTime; }

                /********************************** MULTIPLIERS ***************************************/

                /** Returns the multiplier for the length constraint k **/
                double getLengthMultiplier_k(int k) const { return lagrangianMultiplierLength[k]; }

                /** Returns the multiplier for the source target constraint k,v **/
                double getSourceTargetMultiplier_k(int k, int v) const { return lagrangianMultiplierSourceTarget[k][v]; }

                /** Returns the multiplier for the flow constraint k,v **/
                double getFlowMultiplier_k(int k, int v) const { return lagrangianMultiplierFlow[k][v]; }

                /** Returns the multiplier for the overlap constraint e,s **/
                double getOverlapMultiplier_k(int e, int s) const { return lagrangianMultiplierOverlap[e][s]; }

                double getMaxUsedSliceOverallMultiplier_k(int k) const { return lagrangianMultiplierMaxUsedSliceOverall[k]; }
                double getMaxUsedSliceOverall2Multiplier_k(int e,int s) const { return lagrangianMultiplierMaxUsedSliceOverall2[e][s]; }
                double getMaxUsedSliceOverall3Multiplier_k(int v, int s) const { return lagrangianMultiplierMaxUsedSliceOverall3[v][s]; }

                /********************************** STABILITY CENTER ***************************************/

                /** Returns the stability center for the length constraint k **/
                double getLengthSC_k(int k) const { return lagrangianSCLength[k]; }

                /** Returns the stability center for the source target constraint k,v **/
                double getSourceTargetSC_k(int k, int v) const { return lagrangianSCSourceTarget[k][v]; }

                /** Returns the stability center for the flow constraint k,v **/
                double getFlowSC_k(int k, int v) const { return lagrangianSCFlow[k][v]; }

                /** Returns the stability center for the overlap constraint e,s **/
                double getOverlapSC_k( int e, int s) const { return lagrangianSCOverlap[e][s]; }

                double getMaxUsedSliceOverallSC_k(int k) const { return lagrangianSCMaxUsedSliceOverall[k]; }
                double getMaxUsedSliceOverall2SC_k(int e, int s) const { return lagrangianSCMaxUsedSliceOverall2[e][s]; }
                double getMaxUsedSliceOverall3SC_k(int v, int s) const { return lagrangianSCMaxUsedSliceOverall3[v][s]; }

                /********************************** SLACK ***************************************/

                /** Return the slack for the the length constraint k **/
                double getLengthSlack_k(int k) const { return lengthSlack[k];}

                /** Return the slack for the the source target constraint k,v **/
                double getSourceTargetSlack_k(int k, int v) const { return sourceTargetSlack[k][v];}

                /** Return the slack for the the flow constraint k,v **/
                double getFlowSlack_k(int k,int v) const { return flowSlack[k][v];}

                /** Return the slack for the the overlap constraint e,s **/
                double getOverlapSlack_k(int e, int s) const { return overlapSlack[e][s]; }

                double getMaxUsedSliceOverallSlack_k(int k) const { return maxUsedSliceOverallSlack[k]; }
                double getMaxUsedSliceOverall2Slack_k(int e, int s) const { return maxUsedSliceOverallSlack2[e][s]; }
                double getMaxUsedSliceOverall3Slack_k(int v, int s) const { return maxUsedSliceOverallSlack3[v][s]; }

                /********************************** SLACK PRIMAL VARIABLES ***************************************/
                
                /** Return the slack with primal variables for the the length constraint k **/
                double getLengthSlack_v2_k(int k) const { return lengthSlack_v2[k];}

                /** Return the slack with primal variables for the the source/target constraint k,v **/
                double getSourceTargetSlack_v2_k(int k, int v) const { return sourceTargetSlack_v2[k][v];}

                /** Return the slack with primal variables for the the flow constraint k,v **/
                double getFlowSlack_v2_k(int k,int v) const { return flowSlack_v2[k][v];}

                /** Return the slack with primal variables for the the non-overlap constraint e,s **/
                double getOverlapSlack_v2_k(int e, int s) const { return overlapSlack_v2[e][s]; }

                double getMaxUsedSliceOverallSlack_v2_k(int k) const { return maxUsedSliceOverallSlack_v2[k]; }
                double getMaxUsedSliceOverall2Slack_v2_k(int e, int s) const { return maxUsedSliceOverallSlack2_v2[e][s]; }
                double getMaxUsedSliceOverall3Slack_v2_k(int v, int s) const { return maxUsedSliceOverallSlack3_v2[v][s]; }

                /********************************** DIRECTION ***************************************/

                /** Return the slack with primal variables for the the length constraint k **/
                double getLengthDirection_k(int k) const { return lengthDirection[k];}

                /** Return the slack with primal variables for the the source/target constraint k,v **/
                double getSourceTargetDirection_k(int k, int v) const { return sourceTargetDirection[k][v];}

                /** Return the slack with primal variables for the the flow constraint k,v **/
                double getFlowDirection_k(int k,int v) const { return flowDirection[k][v];}

                /** Return the slack with primal variables for the the non-overlap constraint e,s **/
                double getOverlapDirection_k(int e,int s) const { return overlapDirection[e][s]; }

                double getMaxUsedSliceOverallDirection_k(int k) const { return maxUsedSliceOverallDirection[k]; }
                double getMaxUsedSliceOverall2Direction_k(int e, int s) const { return maxUsedSliceOverall2Direction[e][s]; }
                double getMaxUsedSliceOverall3Direction_k(int v, int s) const { return maxUsedSliceOverall3Direction[v][s]; }

                /********************************** MODULES ***************************************/

                /* Returns |slack|^2 */
                virtual double getSlackModule() = 0;

                virtual double getSlackModule_v2() = 0;

                virtual double get_prod_slack() = 0;

                virtual double getMeanSlackModule_v2() =0;

                virtual double getSlackDirectionProd() = 0;

                /********************************** PRIMAL SOLUTION ***************************************/

                double getPrimalObjective();
                std::vector< std::vector<double> > getPrimalVariables() const { return primal_linear_solution;}

                /********************************** CURRENT SOLUTION ***************************************/

                std::vector< std::vector<bool> > getVariables() const { return assignmentMatrix_d;}

                /********************************** GRAPH D ***************************************/

                std::shared_ptr<ListDigraph> getVecGraphD(int d) const { return vecGraph[d];}

                /***********************************************************************************/
                /*		                       SETTERS 		    	                   */
                /***********************************************************************************/
                void setCurrentLagrCost(double val){ currentLagrCost = val; }
                void incCurrentLagrCost(double val){ currentLagrCost += val; }
                void setCurrentRealCost(double val){ currentRealCost = val; }
                void incCurrentRealCost(double val){ currentRealCost += val; }

                void setConstAuxGraphTime(double value) { constAuxGraphTime = value;}

                /********************************** MULTIPLIERS ***************************************/

                void setLengthMultiplier_k (int k, double val) { lagrangianMultiplierLength[k] = val; }
                void setSourceTargetMultiplier_k (int k, int v, double val) { lagrangianMultiplierSourceTarget[k][v] = val; }
                void setFlowMultiplier_k (int k, int v, double val) { lagrangianMultiplierFlow[k][v] = val; }
                void setOverlapMultiplier_k ( int e, int s, double val) { lagrangianMultiplierOverlap[e][s] = val; }

                void setMaxUsedSliceOverallMultiplier_k(int k, double val) { lagrangianMultiplierMaxUsedSliceOverall[k] = val;}
                void setMaxUsedSliceOverall2Multiplier_k(int e,int s, double val) { lagrangianMultiplierMaxUsedSliceOverall2[e][s] = val;}
                void setMaxUsedSliceOverall3Multiplier_k(int v,int s, double val) { lagrangianMultiplierMaxUsedSliceOverall3[v][s] = val;}

                /********************************** STABILITY CENTER ***************************************/

                void setLengthSC_k (int k, double val) { lagrangianSCLength[k] = val; }
                void setSourceTargetSC_k (int k, int v, double val) { lagrangianSCSourceTarget[k][v] = val; }
                void setFlowSC_k (int k, int v, double val) { lagrangianSCFlow[k][v] = val; }
                void setOverlapSC_k ( int e, int s, double val) { lagrangianSCOverlap[e][s] = val; }

                void setMaxUsedSliceOverallSC_k(int k, double val) { lagrangianSCMaxUsedSliceOverall[k] = val;}
                void setMaxUsedSliceOverall2SC_k(int e,int s, double val) { lagrangianSCMaxUsedSliceOverall2[e][s] = val;}
                void setMaxUsedSliceOverall3SC_k(int v,int s, double val) { lagrangianSCMaxUsedSliceOverall3[v][s] = val;}

                /********************************** SLACK ***************************************/

                void setLengthSlack_k (int k, double val) { lengthSlack[k] = val; }
                void setSourceTargetSlack_k (int k, int v, double val) { sourceTargetSlack[k][v] = val; }
                void setFlowSlack_k (int k, int v, double val) { flowSlack[k][v] = val; }
                void setOverlapSlack_k( int e, int s, double val){ overlapSlack[e][s] = val; }

                void setMaxUsedSliceOverallSlack_k(int k, double val) { maxUsedSliceOverallSlack[k] = val;}
                void setMaxUsedSliceOverall2Slack_k(int e,int s, double val) { maxUsedSliceOverallSlack2[e][s] = val;}
                void setMaxUsedSliceOverall3Slack_k(int v,int s, double val) { maxUsedSliceOverallSlack3[v][s] = val;}

                /********************************** SLACK PRIMAL VARIABLES ***************************************/
                
                void setLengthSlack_v2_k (int k, double val) { lengthSlack_v2[k] = val; }
                void setSourceTargetSlack_v2_k (int k, int v, double val) { sourceTargetSlack_v2[k][v] = val; }
                void setFlowSlack_v2_k (int k, int v, double val) { flowSlack_v2[k][v] = val; }
                void setOverlapSlack_v2_k( int e, int s, double val){ overlapSlack_v2[e][s] = val; }

                void setMaxUsedSliceOverallSlack_v2_k(int k, double val) { maxUsedSliceOverallSlack_v2[k] = val;}
                void setMaxUsedSliceOverall2Slack_v2_k(int e,int s, double val) { maxUsedSliceOverallSlack2_v2[e][s] = val;}
                void setMaxUsedSliceOverall3Slack_v2_k(int v,int s, double val) { maxUsedSliceOverallSlack3_v2[v][s] = val;}

                /********************************** DIRECTION ***************************************/

                void setLengthDirection_k (int k, double val) { lengthDirection[k] = val; }
                void setSourceTargetDirection_k (int k, int v, double val) { sourceTargetDirection[k][v] = val; }
                void setFlowDirection_k (int k, int v, double val) { flowDirection[k][v] = val; }
                void setOverlapDirection_k( int e, int s, double val){ overlapDirection[e][s] = val; }

                void setMaxUsedSliceOverallDirection_k(int k, double val) { maxUsedSliceOverallDirection[k] = val;}
                void setMaxUsedSliceOverall2Direction_k(int e,int s, double val) { maxUsedSliceOverall2Direction[e][s] = val;}
                void setMaxUsedSliceOverall3Direction_k(int v,int s, double val) { maxUsedSliceOverall3Direction[v][s] = val;}

                /****************************************************************************************/
                /*				        UPDATE						*/
                /****************************************************************************************/

                /********************************* MULTIPLIERS ***********************************/

                /* Updates lagrangian multipliers with the rule: u[k+1] = u[k] + t[k]*violation */
                virtual void updateMultiplier(double) = 0;

                void updateLengthMultiplier(double);

                void updateSourceTargetMultiplier(double);

                void updateFlowMultiplier(double);

                void updateOverlapMultiplier(double);

                void updateMaxUsedSliceOverallMultiplier(double);

                void updateMaxUsedSliceOverall2Multiplier(double);

                void updateMaxUsedSliceOverall3Multiplier(double);

                /******** MULTIPLIER CONSIDERING THE STABILITY CENTER ************/

                virtual void updateMultiplier_v2(double) = 0;

                void updateLengthMultiplier_v2(double);

                void updateSourceTargetMultiplier_v2(double);

                void updateFlowMultiplier_v2(double);

                void updateOverlapMultiplier_v2(double);

                void updateMaxUsedSliceOverallMultiplier_v2(double);

                void updateMaxUsedSliceOverall2Multiplier_v2(double);

                void updateMaxUsedSliceOverall3Multiplier_v2(double);

                /********************************* STABILITY CENTER ***********************************/

                virtual void updateStabilityCenter() = 0;

                void updateLengthSC();

                void updateSourceTargetSC();

                void updateFlowSC();

                void updateOverlapSC();

                void updateMaxUsedSliceOverallSC();

                void updateMaxUsedSliceOverall2SC();

                void updateMaxUsedSliceOverall3SC();

                /********************************* SLACK ***********************************/

                virtual void updateSlack() = 0;

                /** Updates the slack of a lehgth constraint using the assigment matrix **/
                void updateLengthSlack();

                /** Updates the slack of a source/target constraint using the assigment matrix **/
                void updateSourceTargetSlack();

                /** Updates the slack of a flow constraint using the assigment matrix **/
                void updateFlowSlack();

                /** Updates the slack of a non-overlap constraint using the assigment matrix **/
                void updateOverlapSlack();

                void updateMaxUsedSliceOverallSlack();

                void updateMaxUsedSliceOverall2Slack();

                void updateMaxUsedSliceOverall3Slack();

                /********************************* SLACK PRIMAL VARIABLES ***********************************/

                virtual void updateSlack_v2() = 0;

                void updateLengthSlack_v2();

                void updateSourceTargetSlack_v2();

                void updateFlowSlack_v2();

                void updateOverlapSlack_v2();

                void updateMaxUsedSliceOverallSlack_v2();

                void updateMaxUsedSliceOverall2Slack_v2();

                void updateMaxUsedSliceOverall3Slack_v2();

                /********************************* DIRECTION ***********************************/

                virtual void updateDirection() = 0;

                void updateLengthDirection(double);

                void updateSourceTargetDirection(double);

                void updateFlowDirection(double);

                void updateOverlapDirection(double);

                void updateMaxUsedSliceOverallDirection(double);

                void updateMaxUsedSliceOverall2Direction(double);

                void updateMaxUsedSliceOverall3Direction(double);
        
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
                virtual ~AbstractLagFormulation(){
                        primal_linear_solution.clear();
                        assignmentMatrix_d.clear();
                }

};

#endif


