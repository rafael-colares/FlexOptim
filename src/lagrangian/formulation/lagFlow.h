#ifndef LAG_FLOW_H
#define LAG_FLOW_H

#include "AbstractLagrangianFormulation.h"

using namespace lemon;

class lagFlow :public AbstractLagFormulation{

    private:

        /*****************************  COSTS *******************************/

        /* Refers to the cost of an arc during iteration k of subgradient. cost = c_{ij} + u_k*length_{ij} */
        std::vector< std::shared_ptr<ArcCost> > cost; 

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

        /********************************* STABILITY CENTER ***********************************/

        void initStabilityCenter();

        /********************************** SLACK ***************************************/

        /** Initializes the slack of relaxed constraints. **/
        void initSlacks();

        /*********** SLACK CONSIDERING PRIMAL VARIABLES ***********/    

        void initSlacks_v2();

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

        void solveProblemMaxUsedSliceOverall();

        /** Checks if all slacks are non-negative. **/
        bool checkFeasibility();

        /****************************************************************************************/
        /*										Getters 										*/
        /****************************************************************************************/

        double getRealCostFromPath(int d, Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &SOURCE, const ListDigraph::Node &TARGET);

        /* Returns the physical length of the path. */
        double getPathLength(int d, Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &s, const ListDigraph::Node &t);
        
        /* Returns the actual cost of the path according to the metric used. */
        double getPathCost(int d, Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &s, const ListDigraph::Node &t);

        /******************************** MODULES *********************************/
        
        double getSlackModule();

        double getSlackModule_v2();

        double getDirectionModule();

        double getSlackDirectionProd();

        double getSlackDirectionProdProjected(Input::ProjectionType);

        double getSlackDirectionProdNormal();

        double get_prod_slack();

        double getMeanSlackModule_v2();

        /****************************************************************************************/
        /*										Setters											*/
        /****************************************************************************************/

        /** Changes the cost of an arc in a graph. @param a The arc. @param d The graph index. @param val The new cost value. **/
        void setArcCost(const ListDigraph::Arc &a, int d, double val) { (*cost[d])[a] = val; }

        /** Increments the cost of an arc in a graph. @param a The arc. @param d The graph index. @param val The value to be added to cost. **/
        void incArcCost(const ListDigraph::Arc &a, int d, double val) { (*cost[d])[a] += val; }

        /****************************************************************************************/
        /*										Update											*/
        /****************************************************************************************/

        /********************************* MULTIPLIERS ***********************************/

        /* Updates lagrangian multipliers with the rule: u[k+1] = u[k] + t[k]*violation */
        void updateMultiplier(double);

        /********************* MULTIPLIER CONSIDERING THE STABILITY CENTER ***************/

        void updateMultiplier_v2(double);

        /***************************** STABILITY CENTER ****************************/

        void updateStabilityCenter();
        
        /********************************* SLACK ***********************************/

        void updateSlack();

        /******************** SLACK CONSIDERING THE PRIMAL VECTOR *******************/

        void updateSlack_v2();

        /*********************************** DIRECTION********************************/

        void updateDirection();

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