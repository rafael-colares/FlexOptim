#include "formulationComponents.h"

Variable::Variable(int identifier, double lowerBound, double upperBound, Type varType, double val, std::string varName): id(identifier), lb(lowerBound), ub(upperBound), type(varType), value(val), name(varName){}

Term::Term(Variable variable, double coefficient): coeff(coefficient), var(variable){}

Expression::Expression(const Expression &e): termsArray(e.getTerms()){}

Constraint::Constraint(double lowerBound, Expression &e, double upperBound, std::string constName): lb(lowerBound), expr(e), ub(upperBound){}

ObjectiveFunction::ObjectiveFunction(Expression &e, Direction d): expr(e), direction(d){}

void Expression::addTerm(Term &term) {
    for (unsigned int i = 0; i < termsArray.size(); i++){
        if (termsArray[i].getVar().getId() == term.getVar().getId()){
            double oldCoeff = termsArray[i].getCoeff();
            double newCoeff = oldCoeff + term.getCoeff();
            termsArray[i].setCoeff(newCoeff);
            return;
        }
    }
    termsArray.push_back(term); 
}