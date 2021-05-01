#include "abstractSolver.h"


/****************************************************************************************/
/*										Constructor										*/
/****************************************************************************************/

/* Constructor. The RSA constructor is called and the arc map storing the index of the preprocessed graphs associated is built. */
AbstractSolver::AbstractSolver(const Instance &instance, const Status &s) : currentStatus(s) {
    FormulationFactory factory;
    formulation = factory.createFormulation(instance);

    totalImpleTime = formulation->getTotalImpleTime();
	varImpleTime = formulation->getVarImpleTime();
	constImpleTime = formulation->getConstImpleTime();
	cutImpleTime = formulation->getCutImpleTime();
	objImpleTime = formulation->getObjImpleTime();
}

void AbstractSolver::updateRSA(Instance &instance){
    std::cout << "Update RSA" << std::endl;
    if (this->getStatus() == STATUS_OPTIMAL || this->getStatus() == STATUS_FEASIBLE){   
        std::cout << "Feasible" << std::endl;
        formulation->updatePath(this->getSolution());
        formulation->updateInstance(instance);
        instance.setWasBlocked(false);
        formulation->displayPaths();

    }
    else{
        std::cout << "Decrease the number of demands to be treated." << std::endl;
        instance.decreaseNbDemandsAtOnce();
        instance.setWasBlocked(true);
    }
}
