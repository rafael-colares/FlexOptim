#ifndef LAG_SOLVER_FACTORY_H
#define LAG_SOLVER_FACTORY_H

// include all concrete solvers
#include "lagSubgradient.h"
#include "lagVolume.h"

class lagSolverFactory{
    public:
        inline AbstractLagSolver* createSolver(const Instance &instance){
            //Input::NodeMethod chosenSolver = instance.getInput().getChosenNodeMethod();
            Input::LagMethod chosenSolver = instance.getInput().getChosenLagMethod();
            switch (chosenSolver){
                //case Input::NODE_METHOD_SUBGRADIENT:{
                case Input::SUBGRADIENT:{
                    return new lagSubgradient(instance);
                    break;
                }
                //case Input::NODE_METHOD_VOLUME:{
                case Input::VOLUME:{
                    return new lagVolume(instance);
                    break;
                }
               
                default:{
                    std::cout << "ERROR: Invalid LAG_Solver." << std::endl;
                    exit(0);
                    break;
                }
            }
            return NULL;
        }
};

#endif