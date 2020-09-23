#include "abstractSolver.h"


/****************************************************************************************/
/*										Constructor										*/
/****************************************************************************************/

/* Constructor. The RSA constructor is called and the arc map storing the index of the preprocessed graphs associated is built. */
AbstractSolver::AbstractSolver(const Instance &instance, const Status &s) : currentStatus(s) {
    FormulationFactory factory;
    formulation = factory.createFormulation(instance);
}