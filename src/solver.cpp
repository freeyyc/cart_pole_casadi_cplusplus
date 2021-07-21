#include "solver.hpp"

Solver::Solver()
{
    solverState = SolverState::NOT_INITIALIZED;
    solution = nullptr;

    minCartHorizonPos = opti.parameter();
    maxCartHorizonPos = opti.parameter();
    minU = opti.parameter();
    maxU = opti.parameter();
    T = opti.parameter();
    initialStateParameters = opti.parameter(4);
    referenceStateParameters = opti.parameter(4);
    opti.set_value(minCartHorizonPos, 0);
    opti.set_value(maxCartHorizonPos, 2);
    opti.set_value(minU, 0);
    opti.set_value(maxU, 20);
    opti.set_value(T, 2);

 }

 casadi::Function Solver::getIntegratorDynamics()
 {
     casadi::MX x = casadi::MX::sym("x", 4);
     casadi::MX a = casadi::MX::sym("a", 2);
     casadi::MX dT = casadi::MX::sym("dt");

     casadi::MX p = x(casadi::Slice(0,2));
     casadi::MX v = x(casadi::Slice(2,4));
     casadi::MX rhs = casadi::MX::vertcat({p + dT * (v + 0.5 * dT * a),
                                           v + dT * a});
    
    return casadi::Function("Integrator", {x,a,dT}, {rhs});
 }

 casadi::Function Solver::getAccelerationConsistencyConstraintFunction()
 {
    casadi::MX X = casadi::MX::sym("x", 4);
    casadi::MX U = casadi::MX::sym("u", 1);
    casadi::MX A = casadi::MX::sym("a", 2);

    casadi::MX currentPosition = X(casadi::Slice(0,2));
    casadi::MX currentVelocity = X(casadi::Slice(2,4));

    cart_pole_model model;
    casadi::MX temp1 = model.l * model.m2 * sin(currentPosition(1)) * pow(currentVelocity(1), 2);
    casadi::MX temp2 = model.m2 * model.g * cos(currentPosition(1)) * sin(currentPosition(1));
    casadi::MX temp3 = model.m1 + model.m2 * (1 - pow(cos(currentPosition(1)), 2));
    casadi::MX constraint_cart = A(0) - (temp1 + U + temp2) / temp3;

    temp1 = model.l * model.m2 * cos(currentPosition(1)) * sin(currentPosition(1)) * pow(currentVelocity(1), 2);
    temp2 = (model.m1 + model.m2) * model.g * sin(currentPosition(1));
    temp3 = model.l * model.m1 + model.l * model.m2 * (1 - pow(cos(currentPosition(1)), 2));
    casadi::MX constraint_pole = A(1) + (temp1 + U * cos(currentPosition(1)) + temp2) / temp3;

    casadi::MX constarint = casadi::MX::vertcat({constraint_cart, constraint_pole});
    return casadi::Function("accelerationConsistency", {X, U, A}, {constarint});

 }

 void Solver::setupOpti()
 {
    using S1 = casadi::Slice;
    casadi_int phaseLength = static_cast<casadi_int> (settings.phaseLength);
    casadi_int N = phaseLength;//single phase

    X = opti.variable(4, N + 1);
    A = opti.variable(2, N);
    U = opti.variable(1, N);
    
    casadi::Function accelerationConstraint = accelerationConsisitencyConstraint;

    opti.subject_to(X(S1(), 0) == initialStateParameters);
    casadi::MX dT = T / phaseLength;
    for(casadi_int k = 0; k < N; ++k){
        opti.subject_to(X(S1(), k + 1) == casadi::MX::vertcat(integratorDynamics({X(S1(), k), A(S1(), k), dT})));
        opti.subject_to(minCartHorizonPos <= X(0, k) <= maxCartHorizonPos);
        if(k == N-1) opti.subject_to(minCartHorizonPos <= X(0, k+1) <= maxCartHorizonPos);
    } 
    for(casadi_int k = 0; k < N; ++k){
        opti.subject_to(minU <= U(0, k) <= maxU);
        opti.subject_to(casadi::MX::vertcat(accelerationConstraint({X(S1(), k), U(S1(), k), A(S1(), k)})) == casadi::MX::zeros(2, 1));
    }
   
    costWeights w = settings._costWeights;
    casadi::MX costFunction = w.control * casadi::MX::sumsqr(U(0, S1()));
    opti.minimize(costFunction);
 }

bool Solver::setupProblem(const Settings& _settings)
{
    settings = _settings;
    integratorDynamics = getIntegratorDynamics();
    accelerationConsisitencyConstraint = getAccelerationConsistencyConstraintFunction();

    casadi_int npoints = static_cast<casadi_int> (settings.phaseLength);
    linSpacePoints = casadi::DM::linspace(0, 1, npoints + 1);
    Xsol.resize(4, npoints + 1);
    Usol.resize(1, npoints);

    setupOpti();

    casadi::Dict casadiOptions;
    casadi::Dict ipoptOptions;

    casadiOptions["expand"] = true;//Replace MX with SX expressions in problem formulation, speed up
    unsigned long solverVerbosity = settings.solverVerbosity;
    if (solverVerbosity) {
        casadi_int ipoptVerbosity = static_cast<long long>(solverVerbosity - 1);
        ipoptOptions["print_level"] = ipoptVerbosity;
        casadiOptions["print_time"] = true;
        casadiOptions["bound_consistency"] = false;
    } else {
        ipoptOptions["print_level"] = 0;
        casadiOptions["print_time"] = false;
        //casadiOptions["bound_consistency"] = false;
        //ipoptOptions["fixed_variable_treatment"] = "make_constraint";
    }
    ipoptOptions["linear_solver"] = settings.ipoptLinearSolver;

    opti.solver("ipopt", casadiOptions, ipoptOptions);

    solverState = SolverState::PROBLEM_SET;

    return true;
}

void Solver::setParametersValue(const State& initialState, const References& references)
{
    opti.set_value(initialStateParameters, initialState.state);
    opti.set_value(referenceStateParameters, references.desiredState);
}

bool Solver::solve(const State& initialState, const References& references)
{
    if(solverState == SolverState::NOT_INITIALIZED){
        std::cerr << "problem not initialized" << std::endl;
        return false;
    }
    
    if(solverState == SolverState::PROBLEM_SOLVED){
        std::cout << "solve success" << std::endl;
        assert(solution);
        // opti.set_initial(solution -> value_variables());
        // opti.set_initial(opti.lam_g(), solution -> value(opti.lam_g()));
    }
    else{
        setParametersValue(initialState, references);
        casadi_int npoints = static_cast<casadi_int> (settings.phaseLength);
        casadi::DM initPos = casadi::DM::zeros(2,1);
        casadi::DM initVel = casadi::DM::zeros(2,1);
        casadi::DM refPos = casadi::DM::zeros(2,1);
        initPos = initialState.state(casadi::Slice(0,2));
        initVel = initialState.state(casadi::Slice(2,4));
        refPos = references.desiredState(casadi::Slice(0,2));
        opti.set_initial(X(casadi::Slice(), 0), casadi::DM::vertcat({initPos, initVel}));
        casadi::DM interpolatedPosition(2, 1);
        for(casadi_int k = 1; k < npoints; ++k){
            //initial guess
            interpolatedPosition = initPos + linSpacePoints(k) * (refPos - initPos);
            opti.set_initial(X(casadi::Slice(0,2), k), interpolatedPosition);
            opti.set_initial(X(casadi::Slice(2,4), k), 0);
            opti.set_initial(U(casadi::Slice(), k), 0);
        }

        
    }
    solverState = SolverState::PROBLEM_SET;
    solution = nullptr;
    
    try{
        solution = std::make_unique<casadi::OptiSol>(opti.solve());
    }catch(std::exception &e){
        opti.debug().show_infeasibilities(1e-5);
        std::cerr << "error while solving the optimization" << std::endl;
        std::cerr << "Details:\n " << e.what() << std::endl;
        return false;
    }
    solverState = SolverState::PROBLEM_SOLVED;
    return true;
}

void Solver::fillSolution(){
    Xsol = solution -> value(X);
    Usol = solution -> value(U);
}