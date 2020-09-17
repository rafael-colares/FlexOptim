#ifndef __FlowForm__h
#define __FlowForm__h

#include "RSA.h"
#include "FormulationComponents.h"

typedef std::vector<Variable> VarArray;
typedef std::vector<VarArray> VarMatrix;

/*********************************************************************************************
* This class implements the Online Routing and Spectrum Allocation through a Flow MIP 
* formulation using CPLEX.	
*********************************************************************************************/
class FlowForm : public RSA{

private:
    VarMatrix x;				    /**< The matrix of assignement variables used in the MIP. x[d][a]=1 if the d-th demand is routed through the arc from index a. **/
    VarArray maxSlicePerLink;	    /**< The array of variables used in the MIP for verifying the max used slice position for each link in the topology network. maxSlicePerLink[i]=p if p is the max used slice position from the link with id i. **/
	Variable maxSliceOverall;		/**< The max used slice position throughout all the network. **/
    int nbVar;                      /**< The total number of variables. **/

    std::vector<Constraint> constraintSet;
    std::vector<ObjectiveFunction> objectiveSet;
public:

	/****************************************************************************************/
	/*										Constructors									*/
	/****************************************************************************************/

	/** Constructor. Builds the Online RSA mixed-integer program and solves it using CPLEX.  @param instance The instance to be solved. **/
    FlowForm(const Instance &instance);

	/****************************************************************************************/
	/*										Getters											*/
	/****************************************************************************************/
    
	/** Returns the total number of variables. **/
    int getNbVar() const { return nbVar; }

	/** Returns the array of variables. **/
	VarArray getVariables();

	/** Returns the array of constraints. **/
	std::vector<Constraint> getConstraints(){ return constraintSet; }

	/** Returns the objective function expression. **/
    ObjectiveFunction getObjFunction(int i){ return objectiveSet[i]; }

	/** Returns the objective function expression. **/
    std::vector<ObjectiveFunction> getObjectiveSet(){ return objectiveSet; }

	/** Returns the objective function expression. **/
    int getNbObjectives(){ return objectiveSet.size(); }

	/** Returns the objective function expression for the given metric. **/
    Expression getObjFunctionFromMetric(Input::ObjectiveMetric chosenObjective);

	/** Returns the source constraint associated with a demand and a node. @param demand The demand. @param d The demand index. @param nodeLabel The node label. **/
    Constraint getSourceConstraint_d_n(const Demand & demand, int d, int nodeLabel);
	
	/** Returns the flow conservation constraint associated with a demand and a node. @param v The node. @param demand The demand. @param d The demand index. **/
    Constraint getFlowConservationConstraint_i_d(ListDigraph::Node &v, const Demand & demand, int d);

	/** Returns the target constraint associated with a demand. @param demand The demand. @param d The demand index. **/
    Constraint getTargetConstraint_d(const Demand & demand, int d);
	
    /** Returns the length constraint associated with a demand. @param demand The demand. @param d The demand index. **/
    Constraint getLengthConstraint(const Demand &demand, int d);

	/** Returns the non-overlapping constraint associated with an edge and a slice. @param linkLabel The arc label. @param slice The arc slice. **/
	Constraint getNonOverlappingConstraint(int linkLabel, int slice);
	
	/** Returns the Link's Max Used Slice Position constraint associated with a link and a demand. @param linkIndex The link's index. @param d The demand's index. **/
	Constraint getMaxUsedSlicePerLinkConstraints(int linkIndex, int d);

	/** Returns the Overall Max Used Slice Position constraints associated with a link and a demand. @param linkIndex The link's index. @param d The demand's index. **/
	Constraint getMaxUsedSliceOverallConstraints(int d);
	/** Returns the Overall Max Used Slice Position constraints associated with a link and a demand. @param linkIndex The link's index. **/
	Constraint getMaxUsedSliceOverallConstraints2(int linkIndex, int s);
	/** Returns the Overall Max Used Slice Position constraints associated with a link and a demand. @param linkIndex The link's index. **/
	Constraint getMaxUsedSliceOverallConstraints3(int linkIndex, int s);
	
	/****************************************************************************************/
	/*										Setters											*/
	/****************************************************************************************/

	/** Increments the total number of variables. **/
    void incNbVar() { nbVar++; }

	/** Defines the decision variables need in the MIP formulation. **/
    void setVariables();

	void setVariableValues(const std::vector<double> &value);

	/** Defines the objective function. @param chosenObjective The chosen objective metrics. **/
    void setObjective(std::vector<Input::ObjectiveMetric> chosenObjective);

	/** Defines Source constraints. At most one arc leaves each node and exactly one arc leaves the source. **/
    void setSourceConstraints();

	/** Defines Flow Conservation constraints. If an arc enters a node, then an arc must leave it. **/
    void setFlowConservationConstraints();
	
	/** Defines Target constraints. Exactly one arc enters the demand's target. **/
    void setTargetConstraints();
	
	/** Defines Length constraints. Demands must be routed within a length limit. **/
    void setLengthConstraints();

	/** Defines Non-Overlapping constraints. Demands must not overlap eachother's slices. **/
    void setNonOverlappingConstraints();
	
	/** Defines the Link's Max Used Slice Position constraints. The max used slice position on each link must be greater than every slice position used in the link. **/
	void setMaxUsedSlicePerLinkConstraints();

	/** Defines the Overall Max Used Slice Position constraints. The max used slice position overall must be greater than every other slice position used in the network. **/
	void setMaxUsedSliceOverallConstraints();
	/** Defines the Overall Max Used Slice Position constraints. The max used slice position overall must be greater than every other slice position used in the network. **/
	void setMaxUsedSliceOverallConstraints2();
	/** Defines the Overall Max Used Slice Position constraints. The max used slice position overall must be greater than every other slice position used in the network. **/
	void setMaxUsedSliceOverallConstraints3();

	/****************************************************************************************/
	/*										Methods											*/
	/****************************************************************************************/

	/** Recovers the obtained MIP solution and builds a path for each demand on its associated graph from RSA. **/
    void updatePath(const std::vector<double> &vals);

	/****************************************************************************************/
	/*										Display											*/
	/****************************************************************************************/

	/** Displays the obtained paths. **/
    void displayOnPath();

	/** Displays the value of each variable in the obtained solution. **/
    void displayVariableValues();

	/****************************************************************************************/
	/*										Destructor										*/
	/****************************************************************************************/

	/** Destructor. Clears the variable matrices, cplex model and environment. **/
	~FlowForm();

};


#endif