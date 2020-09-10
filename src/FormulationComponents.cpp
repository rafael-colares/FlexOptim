#include "FormulationComponents.h"

Variable::Variable(double lowerBound, double upperBound, Type varType): lb(lowerBound), ub(upperBound), type(varType){}

Term::Term(Variable variable, double coefficient): coeff(coefficient), var(variable){}

Expression::Expression(const Expression &e): termsArray(e.getTerms()){}

Constraint::Constraint(double lowerBound, Expression &e, double upperBound): lb(lowerBound), expr(e), ub(upperBound){}