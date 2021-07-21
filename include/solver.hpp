#ifndef SOLVER_H
#define SOLVER_H

#include <common.hpp>
#include <vector>
#include <string>
#include <cassert>
#include <memory>

class Solver{
public:
    Solver();
    ~Solver() { };

    enum SolverState{NOT_INITIALIZED=0, PROBLEM_SET, PROBLEM_SOLVED, UNDEFINED};
    SolverState solverState;
    Settings settings;
    casadi::Function integratorDynamics;
    casadi::Function accelerationConsisitencyConstraint;
    casadi::MX initialStateParameters;
    casadi::MX referenceStateParameters;
    casadi::MX X, A, U, T;
    casadi::MX minCartHorizonPos, maxCartHorizonPos, minU, maxU;
    casadi::DM Xsol, Usol;

    casadi::Opti opti;
    std::unique_ptr<casadi::OptiSol> solution;
    casadi::DM linSpacePoints;
    casadi::Function getIntegratorDynamics();
    casadi::Function getAccelerationConsistencyConstraintFunction();
    void setupOpti();
    bool setupProblem(const Settings& _settings);
    void setParametersValue(const State& initialState, const References& references);

    bool solve(const State& initialState, const References& references);
    void fillSolution();
};

#endif