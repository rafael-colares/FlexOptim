#ifndef __subgradient__h
#define __subgradient__h

#include "RSA.h"

using namespace lemon;
// For the moment only k=1 is implemented
class Subgradient : public RSA {
    
private:
    const double EPSILON = 0.0001;
    const int SOURCE;
    const int TARGET;
    const int MAX_LENGTH;
    const int MAX_NB_IT_WITHOUT_IMPROVEMENT;
    const int MAX_NB_IT;

    int iteration;
    int itWithoutImprovement;

    double UB;
    double LB;
    double currentCost;

    bool isFeasible;
    bool isOptimal;

    /* stores the value of the Lagrangian multiplier during each iteration */
    std::vector<double> lagrangianMultiplier;

    /* stores the value of the slack of length constraint (i.e., b - Dx) during each iteration */
    std::vector<double> slack;

    /* stores the value of the step size used for updating the Lagrangian Multiplier during each iteration */
    std::vector<double> stepSize;

    /* stores the value of lambda used for updating the step size for each iteration */
    std::vector<double> lambda;

    /* refers to the cost of an arc during iteration k of subgradient. cost = c_{ij} + u_k*length_{ij} */
    ListDigraph::ArcMap<double> cost;  
             

public:
	/************************************************/
	/*				    Constructors 		   		*/
	/************************************************/
    Subgradient(const Instance &inst);

	/************************************************/
	/*					   Getters 		    		*/
	/************************************************/
    int getIteration() const { return iteration; }
    int getItWithoutImprovement() const { return itWithoutImprovement; }

    double getUB() const { return UB; }
    double getLB() const { return LB; }
    double getCurrentCost() const { return currentCost; }

    bool getIsFeasible() const { return isFeasible; }
    bool getIsOptimal() const { return isOptimal; }

    std::vector<double> getMultiplier() const { return lagrangianMultiplier; }
    std::vector<double> getStepSize() const { return stepSize; }
    std::vector<double> getSlack() const { return slack; }
    std::vector<double> getLambda() const { return lambda; }

    double getMultiplier_k(int k) const { return lagrangianMultiplier[k]; }
    double getStepSize_k(int k) const { return stepSize[k]; }
    double getSlack_k(int k) const { return slack[k]; }
    double getLambda_k(int k) const { return lambda[k]; }

    double getLastMultiplier() const { return lagrangianMultiplier[lagrangianMultiplier.size() - 1]; }
    double getLastStepSize() const { return stepSize[stepSize.size() - 1]; }
    double getLastSlack() const { return slack[slack.size() - 1]; }
    double getLastLambda() const { return lambda[lambda.size() - 1]; }
    
	/************************************************/
	/*					   Setters 		    		*/
	/************************************************/
    void setIteration(int i) { iteration = i; }
    void setItWithoutImprovement(int i) { itWithoutImprovement = i; }
    void incIteration() { iteration++; }
    void incItWithoutImprovement() { itWithoutImprovement++; }

    void setUB(double i){ UB = i; }
    void setLB(double i){ LB = i; }
    void setCurrentCost(double i){ currentCost = i; }

    void setIsFeasible(bool i) { isFeasible = i;}
    void setIsOptimal(bool i) { isOptimal = i; }

    void setMultiplier_k (int k, double i) { lagrangianMultiplier[k] = i; }
    void setStepSize_k(int k, double i){ stepSize[k] = i; }
    void setSlack_k(int k, double i){ slack[k] = i; }
    void setLambda_k(int k, double i){ lambda[k] = i; }

	/************************************************/
	/*					   Methods 		    		*/
	/************************************************/
    /* Sets the initial parameters for the subgradient to run. */
    void initialization();

    /* Call preprocessing functions. */
    void subgradientPreprocessing();

    /* Updates the arc costs according to the last lagrangian multiplier available. cost = c + u_k*length */
    void updateCosts();

    /* Solves the Constrained Shortest Path from node s to node t using the Subgradient Method. */
    void run(const ListDigraph::Node &s, const ListDigraph::Node &t);

    /* Updates the known lower bound. */
    void updateLB(double bound);

    /* Updates the known upper bound. */
    void updateUB(double bound);
    
    /* Updates lagrangian multiplier with the rule: u[k+1] = u[k] + t[k]*violation */
    void updateMultiplier();

    /* Updates the step size with the rule: lambda*(UB - Z[u])/|slack| */
    void updateStepSize();
    
    /* Updates the lambda used in the update of step size. Lambda is halved if LB has failed to increade in some fixed number of iterations */
    void updateLambda();
    
    /* Updates the slack of length constraint for a given path length */
    void updateSlack(double pathLength);
    
    /* Verifies if optimality condition has been achieved and update STOP flag. */
    void updateStop(bool &STOP);

    /* Tests if CSP is feasible by searching for a shortest path with arc costs based on their physical length. */
    bool testFeasibility(const ListDigraph::Node &s, const ListDigraph::Node &t);
    
    /* Assigns length as the main cost of the arcs. */
    void setLengthCost();
    
    /* Stores the path found in the arcMap onPath. */
    void updateOnPath(Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &s, const ListDigraph::Node &t);
    
    /* Returns the physical length of the path. */
    double getPathLength(Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &s, const ListDigraph::Node &t);
    
    /* Returns the actual cost of the path according to the metric used. */
    double getPathCost(Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &s, const ListDigraph::Node &t);
    
    /************************************************/
	/*					   Display 		    		*/
	/************************************************/
    std::string getPathString(Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &s, const ListDigraph::Node &t);
    void displayPath(Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &s, const ListDigraph::Node &t);
    void displayMainParameters(Dijkstra< ListDigraph, ListDigraph::ArcMap<double> > &path, const ListDigraph::Node &s, const ListDigraph::Node &t);
    void displayStepSize();
    void displayMultiplier();
    void displaySlack();
    void displayLambda();
};    
#endif