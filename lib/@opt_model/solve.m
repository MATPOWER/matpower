function [x, f, eflag, output, lambda] = solve(om, opt)
%SOLVE  Solve the optimization model.
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = OM.SOLVE()
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = OM.SOLVE(OPT)
%
%   Solves the optimization model using one of the following, depending
%   on the problem type: QPS_MATPOWER, MIQPS_MATPOWER, NLPS_MATPOWER.
%
%   Inputs:
%       OPT : optional options structure with the following fields,
%           all of which are also optional (default values shown in
%           parentheses)
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           alg ('DEFAULT') : determines which solver to use, list of relevant
%                   problem types are listed in parens next to each
%               'DEFAULT' : automatic, depending on problem type, uses the
%                       the first available of:
%                   LP - Gurobi, CPLEX, MOSEK, Opt Tbx (if MATLAB), GLPK, BPMPD,
%                           MIPS
%                   QP - Gurobi, CPLEX, MOSEK, Opt Tbx (if MATLAB), BPMPD, MIPS
%                   MILP - Gurobi, CPLEX, MOSEK, Opt Tbx (intlingprog), GLPK
%                   MIQP - Gurobi, CPLEX, MOSEK
%                   NLP - MIPS
%                   MINLP - Artelys Knitro (not yet implemented)
%               'MIPS'    : (LP, QP, NLP) MIPS, MATPOWER Interior Point Solver
%                        pure MATLAB implementation of a primal-dual
%                        interior point method, if mips_opt.step_control = 1
%                        it uses MIPS-sc, a step controlled variant of MIPS
%               'BPMPD'   : (LP, QP) BPMPD_MEX
%               'CLP'     : (LP, QP) CLP
%               'CPLEX'   : (LP, QP, MILP, MIQP) CPLEX
%               'FMINCON' : (NLP) FMINCON, MATLAB Optimization Toolbox
%               'GLPK'    : (LP, MILP) GLPK
%               'GUROBI'  : (LP, QP, MILP, MIQP) Gurobi
%               'IPOPT'   : (LP, QP, NLP) IPOPT, requires MEX interface to IPOPT solver
%                           https://github.com/coin-or/Ipopt
%               'KNITRO'  : (NLP, MINLP) Artelys Knitro, requires Artelys Knitro solver
%                           https://www.artelys.com/solvers/knitro/
%               'MOSEK'   : (LP, QP, MILP, MIQP) MOSEK
%               'OT'      : (LP, QP, MILP) MATLAB Optimization Toolbox,
%                           LINPROG, QUADPROG or INTLINPROG
%           bp_opt      - options vector for BP
%           clp_opt     - options vector for CLP
%           cplex_opt   - options struct for CPLEX
%           fmincon_opt - options struct for FMINCON
%           glpk_opt    - options struct for GLPK
%           grb_opt     - options struct for GUROBI
%           intlinprog_opt - options struct for INTLINPROG
%           ipopt_opt   - options struct for IPOPT
%           knitro_opt  - options struct for Artelys Knitro
%           linprog_opt - options struct for LINPROG
%           mips_opt    - options struct for MIPS
%           mosek_opt   - options struct for MOSEK
%           quadprog_opt - options struct for QUADPROG
%           x0 (empty)  - initial value of x, overrides value stored in model
%           skip_prices (0) - flag that specifies whether or not to skip the
%               price computation stage for mixed integer problems, in which
%               the problem is re-solved for only the continuous variables,
%               with all others being constrained to their solved values
%           price_stage_warn_tol (1e-7) - tolerance on the objective fcn
%               value and primal variable relative match required to avoid
%               mis-match warning message if mixed integer price computation
%               stage is not skipped
%
%   Outputs:
%       X : solution vector
%       F : final objective function value
%       EXITFLAG : exit flag
%           1 = converged
%           0 or negative values = solver specific failure codes
%       OUTPUT : output struct with the following fields:
%           alg - algorithm code of solver used
%           (others) - algorithm specific fields
%               iterations - number of iterations performed
%               hist - struct array with trajectories of the following:
%                       feascond, gradcond, compcond, costcond, gamma,
%                       stepsize, obj, alphap, alphad
%               message - exit message
%       LAMBDA : struct containing the Langrange and Kuhn-Tucker
%           multipliers on the constraints, with fields:
%           eqnonlin - nonlinear equality constraints
%           ineqnonlin - nonlinear inequality constraints
%           mu_l - lower (left-hand) limit on linear constraints
%           mu_u - upper (right-hand) limit on linear constraints
%           lower - lower bound on optimization variables
%           upper - upper bound on optimization variables
%
%   See also OPT_MODEL, QPS_MATPOWER, MIQPS_MATPOWER, NLPS_MATPOWER

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 2
    opt = struct();
end

%% linear constraints
[A, l, u] = om.params_lin_constraint();

pt = om.problem_type();
if strcmp(pt, 'MINLP')          %% MINLP - mixed integer non-linear program
    error('@opt_model/solve: not yet implemented for ''MINLP'' problems.')
elseif strcmp(pt, 'NLP')        %% NLP - non-linear program
    %% optimization vars, bounds, types
    [x0, xmin, xmax] = om.params_var();
    if isfield(opt, 'x0')
        x0 = opt.x0;
    end

    %% run solver
    f_fcn = @(x)opf_costfcn(x, om);
    gh_fcn = @(x)opf_consfcn(x, om);
    hess_fcn = @(x, lambda, cost_mult)opf_hessfcn(x, lambda, cost_mult, om);
    [x, f, eflag, output, lambda] = ...
        nlps_matpower(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt);
    success = (eflag > 0);
else                            %% LP, QP, MILP or MIQP
    %% get cost parameters
    [HH, CC, C0] = om.params_quad_cost();

    if strcmp(pt(1:2), 'MI')    %% MILP, MIQP - mixed integer linear/quadratic program
        %% optimization vars, bounds, types
        [x0, xmin, xmax, vtype] = om.params_var();
        if isfield(opt, 'x0')
            x0 = opt.x0;
        end

        %% run solver
        [x, f, eflag, output, lambda] = ...
            miqps_matpower(HH, CC, A, l, u, xmin, xmax, x0, vtype, opt);
    else                        %% LP, QP - linear/quadratic program
        %% optimization vars, bounds, types
        [x0, xmin, xmax] = om.params_var();
        if isfield(opt, 'x0')
            x0 = opt.x0;
        end

        %% run solver
        [x, f, eflag, output, lambda] = ...
            qps_matpower(HH, CC, A, l, u, xmin, xmax, x0, opt);
    end
    f = f + C0;
    success = (eflag == 1);
end
