#ifndef __FormulationComponents__h
#define __FormulationComponents__h

#include <string>
#include <vector>
/********************************************************************************************
 * This class identifies a Variable in a MIP formulation. A demand is defined by its lower 
 * and upper bounds, and its type. 
 ********************************************************************************************/
class Variable
{    
    
	/** Enumerates the possible types of variable. **/
    enum Type {						
		TYPE_REAL = 0,    /**< Variable can admit real values. **/
		TYPE_BOOLEAN = 1, /**< Variable can only admit boolean values. **/
		TYPE_INTEGER = 2  /**< Variable can only admit integer values. **/
	};

private:
    double lb;      /**< The variable lower bound. **/
    double ub;      /**< The variable upper bound. **/
    Type type;      /**< The variable type. **/

public:
    /** Constructor. @param lowerBound The variable lower bound. @param upperBound The variable upper bound. @param varType The variable type. **/
    Variable(double lowerBound=0, double upperBound=1, Type varType=TYPE_BOOLEAN);

	/****************************************************************************************/
	/*										Getters											*/
	/****************************************************************************************/
	/** Returns the variable's lower bound. **/
    double getLb() const { return lb; }

	/** Returns the variable's upper bound. **/
    double getUb() const { return ub; }

	/** Returns the variable's type. **/
    Type getType() const { return type; }

	/****************************************************************************************/
	/*										Setters											*/
	/****************************************************************************************/
	/** Changes the variable's lower bound. @param lowerBound The new lower bound. **/
    void setLb(double lowerBound) { this->lb = lowerBound; }

	/** Changes the variable's upper bound. @param upperBound The new upper bound. **/
    void setUb(double upperBound) { this->ub = upperBound; }

	/** Changes the variable's type.  @param varType The new type. **/
    void setType(Type varType) { this->type = varType; }

};


/***************************************************************************************************
 * This class identifies a Term in an Expression. A term is defined by a Variable and a coefficent. 
 ***************************************************************************************************/
class Term
{    

private:
    double coeff;      /**< The term coefficient. **/
    Variable var;      /**< The term variable. **/

public:
    /** Constructor. @param variable The term variable. @param coefficient The term coeffiencient. **/
    Term(Variable variable, double coefficient);

	/****************************************************************************************/
	/*										Getters											*/
	/****************************************************************************************/
	/** Returns the term's coefficient. **/
    double getCoeff() const { return coeff; }

	/** Returns the term's variable. **/
    Variable getVar() const { return var; }

	/****************************************************************************************/
	/*										Setters											*/
	/****************************************************************************************/
	/** Changes the term's coefficient. @param coefficient The new coefficient. **/
    void setCoeff(double coefficient) { this->coeff = coefficient; }

	/** Changes the term's variable. @param v The new variable. **/
    void setVar(Variable v) { this->var = v; }
};


/********************************************************************************************
 * This class identifies an Expression in a MIP formulation. An Expression is defined by a 
 * vector of Term. 
 ********************************************************************************************/
class Expression
{    

private:
    std::vector<Term> termsArray;  /**< The array of terms. **/

public:
    /** Constructor. @param lowerBound The variable lower bound. @param upperBound The variable upper bound. @param varType The variable type. **/
    Expression();

    /** Copy constructor. @param e The expression to be copied. **/
    Expression(const Expression &e);

	/****************************************************************************************/
	/*										Getters											*/
	/****************************************************************************************/
	/** Returns the vector of terms. **/
    std::vector<Term> getTerms() const { return termsArray; }

	/** Returns the i-th term. @param pos Term position. **/
    Term getTerm_i(int pos) const { return termsArray[pos]; }

	/****************************************************************************************/
	/*										Setters											*/
	/****************************************************************************************/
	/** Changes the vector of terms. @param vec The new vector of terms. **/
    void setTerms(std::vector<Term> &vec) { this->termsArray = vec; }

	/** Changes the i-th term. @param pos Term position. @param term The new term. **/
    void setTerm_i(int pos, Term &term) { this->termsArray[pos] = term; }

	/****************************************************************************************/
	/*										Methods											*/
	/****************************************************************************************/
	/** Adds a new term to the expression. @param term The term to be added. **/
    void addTerm(Term &term) { termsArray.push_back(term); }

    /** Clear the expression. **/
    void clear() { termsArray.clear(); }
};



/********************************************************************************************
 * This class identifies an Expression in a MIP formulation. An Expression is defined by a 
 * vector of Term. 
 ********************************************************************************************/
class Constraint
{    

private:
    double lb;          /**< The constraint lower bound. **/
    Expression expr;    /**< The constraint expression. **/
    double ub;          /**< The constraint upper bound. **/

public:
    /** Constructor. @param lowerBound The variable lower bound. @param upperBound The variable upper bound. @param varType The variable type. **/
    Constraint(double lowerBound, Expression &e, double upperBound);

	/****************************************************************************************/
	/*										Getters											*/
	/****************************************************************************************/
	/** Returns the constraint's lower bound. **/
    double getLb() const { return lb; }

	/** Returns the constraint's upper bound. **/
    double getUb() const { return ub; }
    
	/** Returns the constraint's expression. **/
    Expression getExpression() const { return expr; }

	/****************************************************************************************/
	/*										Setters											*/
	/****************************************************************************************/
	/** Changes the constraint's lower bound. @param lowerBound The new lower bound. **/
    void setLb(double lowerBound) { this->lb = lowerBound; }

	/** Changes the constraint's upper bound. @param upperBound The new upper bound. **/
    void setUb(double upperBound) { this->ub = upperBound; }

	/** Changes the constraint's expression. @param e The new expression. **/
    void setExpression(Expression &e) {this->expr = e; }
    
    /** Clear the constriant's expression. **/
    void clear() { this->expr.clear(); }
};


#endif