#include "solver.hpp"

int main(){
    Settings settings;
    settings.phaseLength = 10;
    settings._costWeights.control = 1;
    settings.solverVerbosity = 0;
    settings.ipoptLinearSolver = "mumps";
    
    cart_pole_model model;
    model.l = 0.5;
    model.m1 = 1;
    model.m2 = 0.3;
    model.g = 9.81;

    State initialState;
    initialState.state.zeros(4);
    References references;
    references.desiredState = {1, 3.14, 0, 0};
    
    Solver solver;
    solver.setupProblem(settings);
    bool ok = solver.solve(initialState, references);
    
    if(ok) solver.fillSolution();
    else ok = solver.solve(initialState, references);
    if(ok) solver.fillSolution();

    return 0;
}