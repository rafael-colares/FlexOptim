#ifndef LAG_FLOW_H
#define LAG_FLOW_H

#include "AbstractLagrangianFormulation.h"

using namespace lemon;

class lagFlow :public AbstractLagFormulation{

    private:

        /********************************* MULTIPLIERS ***********************************/

        /** A vector storing the value of the Lagrangian multipliers associated with Length Constraints. **/
        std::vector<double> lagrangianMultiplierLength;

        /** A 2-dimensional vector storing the value of the Lagrangian multipliers associated with each Non-Overlapping constraint. **/
        std::vector< std::vector<double> > lagrangianMultiplierOverlap;

        /********************************* STABILITY CENTER ***********************************/

        /** A vector storing the value of the Lagrangian multipliers associated with Length Constraints. **/
        std::vector<double> lagrangianSCLength;

        /** A 2-dimensional vector storing the value of the Lagrangian multipliers associated with each Non-Overlapping constraint. **/
        std::vector< std::vector<double> > lagrangianSCOverlap;

        /********************************** SLACK ***************************************/

        /** Stores the value of the slack of lengths constraints (i.e., b - Dx). */
        std::vector<double> lengthSlack;

        /** Stores the value of the slack of Non-Overlapping constraints (i.e., b - Dx). */
        std::vector< std::vector<double> > overlapSlack;

        /******* SLACK CONSIDERING PRIMAL VARIABLES ***********/
        
        std::vector<double> lengthSlack_v2;

        std::vector< std::vector<double> > overlapSlack_v2;

        /*****************************  DIRECTION *******************************/

        std::vector<double> lengthDirection;

        std::vector< std::vector<double> > overlapDirection;

        /*****************************  COSTS *******************************/

        /* Refers to the cost of an arc during iteration k of subgradient. cost = c_{ij} + u_k*length_{ij} */
        std::vector< std::shared_ptr<ArcCost> > cost; 

        /***************************** ASSIGNMENT MATRIX *******************************/
        
        //std::vector< std::vector<bool> > assignmentMatrix_d;

    public:

        /* *******************************************************************************
        *                             INITIALIZATION METHODS
        ******************************************************************************* */
        lagFlow(const Instance &instance):AbstractLagFormulation(instance){}

        /** Sets all initial parameters **/
        void init();

        /********************************* MULTIPLIERS ***********************************/

        /** Sets the initial lagrangian multipliers values for the subgradient to run. **/
        void initMultipliers();
        
        /** Sets the initial lagrangian multipliers associated with length constraints. **/
        void initializeLengthMultipliers();
        
        /** Sets the initial lagrangian multipliers associated with non-overlapping constraints. **/
        void initializeOverlapMultipliers();

        /********************************* STABILITY CENTER ***********************************/

        void initStabilityCenter();

        void initializeLengthSC();

        void initializeOverlapSC();

        /********************************** SLACK ***************************************/

        /** Initializes the slack of relaxed constraints. **/
        void initSlacks();

        /** Initializes the slack of length constraints. **/
        void initializeLengthSlacks();
        
        /** Initializes the slack of non-overlap constraints. **/
        void initializeOverlapSlacks();

        /*********** SLACK CONSIDERING PRIMAL VARIABLES ***********/    

        void initSlacks_v2();

        void initializeLengthSlacks_v2();
        
        void initializeOverlapSlacks_v2();

        /******************************** DIRECTION *********************************/

        void initDirection();

        /******************************** COSTS *********************************/

        /** Initializes the costs in the objective function. **/
        void initCosts();

        /******************************** ASSIGMENT MATRIX *********************************/

        /** Initializes the assignement matrix. **/
        void initAssignmentMatrix();

        /* *******************************************************************************
        *                             RUNING METHODS
        ******************************************************************************* */

        /** Solves the RSA using the Subgradient Method. **/
        void run();

        void subtractConstantValuesFromLagrCost();

        /** Checks if all slacks are non-negative. **/
        bool checkFeasibility();
        
        /* Checks if all length slacks are non-negative. */
        bool checkLengthFeasibility();

        /* Checks if all overlap1 slacks are non-negative. */
        bool checkOverlapFeasibility();

        /****************************************************************************************/
        /*										Getters 										*/
        /****************************************************************************************/

        double getLengthMultiplier_k(int k) const { return lagrangianMultiplierLength[k]; }

        double getOverlapMultiplier_k( int e, int s) const { return lagrangianMultiplierOverlap[e][s]; }

        double getLengthSC_k(int k) const { return lagrangianSCLength[k]; }

        double getOverlapSC_k( int e, int s) const { return lagrangianSCOverlap[e][s]; }
      
        double getOverlapSlack_k(int e, int s) const { return overlapSlack[e][s]; }

        double getLengthSlack_k(int k) const { return lengthSlack[k]; }

        double getOverlapSlack_v2_k(int e, int s) const { return overlapSlack_v2[e][s]; }

        double getLengthSlack_v2_k(int k) const { return lengthSlack_v2[k]; }

        double getOverlapDirection_k(int e, int s) const { return overlapDirection[e][s]; }

        double getLengthDirection_k(int k) const { return lengthDirection[k]; }

        double getRealCostFromPath(int d, Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &SOURCE, const ListDigraph::Node &TARGET);

        double getSlackModule();

        double getSlackModule_v2();

        double getDirectionModule();

        double getSlackDirectionProd();

        double getSlackDirectionProdProjected(Input::ProjectionType);

        double getSlackDirectionProdNormal();

        double get_prod_slack();

        double getMeanSlackModule_v2();

        /* Returns the physical length of the path. */
        double getPathLength(int d, Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &s, const ListDigraph::Node &t);
        
        /* Returns the actual cost of the path according to the metric used. */
        double getPathCost(int d, Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &s, const ListDigraph::Node &t);

        /****************************************************************************************/
        /*										Setters											*/
        /****************************************************************************************/

        void setLengthMultiplier_k (int k, double val) { lagrangianMultiplierLength[k] = val; }
        void setOverlapMultiplier_k ( int e, int s, double val) { lagrangianMultiplierOverlap[e][s] = val; }
       
        void setLengthSC_k (int k, double val) { lagrangianSCLength[k] = val; }
        void setOverlapSC_k ( int e, int s, double val) { lagrangianSCOverlap[e][s] = val; }
       
        /** Changes the value of the slack from a length constraint. @param d The demand index. @param val The new slack value. **/
        void setLengthSlack_k(int d, double val){ lengthSlack[d] = val; }
        void setOverlapSlack_k( int e, int s, double val){ overlapSlack[e][s] = val; }

        void setLengthSlack_v2_k(int d, double val){ lengthSlack_v2[d] = val; }
        void setOverlapSlack_v2_k( int e, int s, double val){ overlapSlack_v2[e][s] = val; }

        void setLengthDirection_k(int d, double val){ lengthDirection[d] = val; }
        void setOverlapDirection_k( int e, int s, double val){ overlapDirection[e][s] = val; }

        /** Changes the cost of an arc in a graph. @param a The arc. @param d The graph index. @param val The new cost value. **/
        void setArcCost(const ListDigraph::Arc &a, int d, double val) { (*cost[d])[a] = val; }

        /** Increments the cost of an arc in a graph. @param a The arc. @param d The graph index. @param val The value to be added to cost. **/
        void incArcCost(const ListDigraph::Arc &a, int d, double val) { (*cost[d])[a] += val; }
    

        /****************************************************************************************/
        /*										Update											*/
        /****************************************************************************************/
        
        /********************************* SLACK ***********************************/

        void updateSlack();

        /** Updates the slack of a lehgth constraint using the assigment matrix **/
        void updateLengthSlack();
        
        void updateOverlapSlack();

        /******** SLACK CONSIDERING THE PRIMAL VECTOR ************/

        void updateSlack_v2();

        void updateLengthSlack_v2();
        
        void updateOverlapSlack_v2();

        /*********************************** DIRECTION********************************/

        void updateDirection();

        void updateLengthDirection(double);

        void updateOverlapDirection(double);


        /********************************* MULTIPLIERS ***********************************/

        /* Updates lagrangian multipliers with the rule: u[k+1] = u[k] + t[k]*violation */
        void updateMultiplier(double);

        /******** MULTIPLIER CONSIDERING THE STABILITY CENTER ************/

        void updateMultiplier_v2(double);

        void updateStabilityCenter();

        /********************************* COSTS ***********************************/

        /** Updates the arc costs. @note cost = coeff + u_d*length **/
        void updateCosts();

        /********************************* ASSIGNMENT MATRIX ***********************************/

        /** Updates the assignment of a demand based on the a given path. @param d The d-th demand. @param path The path on which the demand is routed. @param SOURCE The path's source. @param TARGET The path's target. **/
        void updateAssignment_k(int d, Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &SOURCE, const ListDigraph::Node &TARGET);

        
        /****************************************************************************************/
        /*										Display											*/
        /****************************************************************************************/
    
        /* Verifies if optimality condition has been achieved and update STOP flag. */
        void updateStop(bool &STOP);

        /* Tests if CSP is feasible by searching for a shortest path with arc costs based on their physical length. */
        bool testFeasibility(const ListDigraph::Node &s, const ListDigraph::Node &t);
        
        /* Assigns length as the main cost of the arcs. */
        void setLengthCost();
        
        /* Stores the solution found in the arcMap onPath. */
        void updateOnPath();
        
        std::string getPathString(int d, Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &s, const ListDigraph::Node &t);
          
        void displayMultiplier(std::ostream & = std::cout);

        void displaySlack(std::ostream & = std::cout);

        /* *******************************************************************************
        *                             DESTRUCTOR
        ******************************************************************************* */

        ~lagFlow();
        
};    

#endif