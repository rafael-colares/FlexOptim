#ifndef ABSTRACT_HEURISTIC_H
#define ABSTRACT_HEURISTIC_H

#include "formulation/AbstractLagrangianFormulation.h"
#include <set>
//#include "formulation/lagFormulationFactory.h"

class AbstractLagFormulation;

class AbstractHeuristic{

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

        AbstractLagFormulation* formulation;

        double currentHeuristicCost;

        /* A vector for the heuristic solution */
        std::vector<std::vector<bool>> heuristicSolution;

        Status statusheuristic;

    public:

        /* *******************************************************************************
        *                              CONSTRUCTOR
        ******************************************************************************* */

        AbstractHeuristic(AbstractLagFormulation* form, const Status &s = STATUS_UNKNOWN);

        /* *******************************************************************************
        *                              GETTERS
        ******************************************************************************* */

        /** Returns the status. **/
        Status getStatus() const {return statusheuristic;}

        /* Get the current heuristic solution cost */
        double getCurrentHeuristicCost(){ return currentHeuristicCost;}

        /* *******************************************************************************
        *                              SETTERS
        ******************************************************************************* */

        /** Changes the status. **/
        void setStatus(const Status &s){ statusheuristic = s; }

        /* Set the current heuristic solution cost */
        void setCurrentHeuristicCost(int val){currentHeuristicCost = val;}

        void incCurrentHeuristicCost(int val){currentHeuristicCost += val;}

        /* *******************************************************************************
        *                                 METHODS
        ******************************************************************************* */

        virtual void init() = 0;

        void initSolution();

        virtual void run() = 0;

        void updateCostFromHeuristic();

        void printSolution();

        /* *******************************************************************************
        *                                DESTRUCTOR
        ******************************************************************************* */

        virtual ~AbstractHeuristic(){
            heuristicSolution.clear();
        }

};

#endif