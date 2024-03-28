function [x, f, eflag, output, lambda] = miqps_master(H, c, A, l, u, xmin, xmax, x0, vtype, opt)
% miqps_master - Mixed Integer Quadratic Program Solver wrapper function.
% ::
%
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       MIQPS_MASTER(H, C, A, L, U, XMIN, XMAX, X0, VTYPE, OPT)
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = MIQPS_MASTER(PROBLEM)
%   A common wrapper function for various MILP/MIQP solvers. 
%   Solves the following MILP/MIQP (mixed integer linear programming/
%   mixed integer quadratic programming) problem:
%
%       min 1/2 X'*H*X + C'*X
%        X
%
%   subject to
%
%       L <= A*X <= U       (linear constraints)
%       XMIN <= X <= XMAX   (variable bounds)
%       X(i) is integer, for i in I (integer variable constraints)
%       X(b) is binary, for b in B (binary variable constraints)
%
%   Inputs (all optional except H, C, A and L):
%       H : matrix (possibly sparse) of quadratic cost coefficients
%       C : vector of linear cost coefficients
%       A, L, U : define the optional linear constraints. Default
%           values for the elements of L and U are -Inf and Inf,
%           respectively.
%       XMIN, XMAX : optional lower and upper bounds on the
%           X variables, defaults are -Inf and Inf, respectively.
%       X0 : optional starting value of optimization vector X
%       VTYPE : character string of length NX (number of elements in X),
%               or 1 (value applies to all variables in x),
%               allowed values are 'C' (continuous), 'B' (binary),
%               'I' (integer), 'S' (semi-continuous), or 'N' (semi-integer).
%               (MOSEK, GLPK, OT allow only 'C', 'B', or 'I')
%       OPT : optional options structure with the following fields,
%           all of which are also optional (default values shown in
%           parentheses)
%           alg ('DEFAULT') : determines which solver to use, can be either
%                   a string (new-style) or a numerical alg code (old-style)
%               'DEFAULT' : (or 0) automatic, first available of Gurobi,
%                       CPLEX, MOSEK, Opt Tbx (MILPs only), GLPK (MILPs only)
%               'CPLEX'   : (or 500) CPLEX
%               'GLPK'    : GLPK, (MILP problems only, i.e. empty H matrix)
%               'GUROBI'  : (or 700) Gurobi
%               'MOSEK'   : (or 600) MOSEK
%               'OT'      : (or 300) Optimization Toolbox, INTLINPROG
%                           (MILP problems only, i.e. empty H matrix)
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           skip_prices (0) - flag that specifies whether or not to
%               skip the price computation stage, in which the problem
%               is re-solved for only the continuous variables, with all
%               others being constrained to their solved values
%           price_stage_warn_tol (1e-7) - tolerance on the objective fcn
%               value and primal variable relative match required to avoid
%               mis-match warning message
%           cplex_opt - options struct for CPLEX
%           glpk_opt    - options struct for GLPK
%           grb_opt   - options struct for GUROBI
%           intlinprog_opt - options struct for INTLINPROG
%           linprog_opt - options struct for LINPROG
%           mosek_opt - options struct for MOSEK
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: H, c, A, l, u, xmin, xmax, x0, vtype, opt
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
%       LAMBDA : struct containing the Langrange and Kuhn-Tucker
%           multipliers on the constraints, with fields:
%           mu_l - lower (left-hand) limit on linear constraints
%           mu_u - upper (right-hand) limit on linear constraints
%           lower - lower bound on optimization variables
%           upper - upper bound on optimization variables
%
%   Note the calling syntax is almost identical to that of QUADPROG
%   from MathWorks' Optimization Toolbox. The main difference is that
%   the linear constraints are specified with A, L, U instead of
%   A, B, Aeq, Beq.
%
%   Calling syntax options:
%       [x, f, exitflag, output, lambda] = ...
%           miqps_master(H, c, A, l, u, xmin, xmax, x0, vtype, opt)
%
%       x = miqps_master(H, c, A, l, u)
%       x = miqps_master(H, c, A, l, u, xmin, xmax)
%       x = miqps_master(H, c, A, l, u, xmin, xmax, x0)
%       x = miqps_master(H, c, A, l, u, xmin, xmax, x0, vtype)
%       x = miqps_master(H, c, A, l, u, xmin, xmax, x0, vtype, opt)
%       x = miqps_master(problem), where problem is a struct with fields:
%                       H, c, A, l, u, xmin, xmax, x0, vtype, opt
%                       all fields except 'c', 'A' and 'l' or 'u' are optional
%       x = miqps_master(...)
%       [x, f] = miqps_master(...)
%       [x, f, exitflag] = miqps_master(...)
%       [x, f, exitflag, output] = miqps_master(...)
%       [x, f, exitflag, output, lambda] = miqps_master(...)
%
%   Example: (problem from from %% from MOSEK 6.0 Guided Tour, section  7.13.1
%             https://docs.mosek.com/6.0/toolbox/node009.html)
%       c = [-2; -3];
%       A = sparse([195 273; 4 40]);
%       u = [1365; 140];
%       xmax = [4; Inf];
%       vtype = 'I';
%       opt = struct('verbose', 2);
%       p = struct('c', c, 'A', A, 'u', u, 'xmax', xmax, 'vtype', vtype, 'opt', opt);
%       [x, f, s, out, lam] = miqps_master(p);

%   MP-Opt-Model
%   Copyright (c) 2010-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%%----- input argument handling  -----
%% gather inputs
if nargin == 1 && isstruct(H)       %% problem struct
    p = H;
    if isfield(p, 'opt'),   opt = p.opt;    else,   opt = [];   end
    if isfield(p, 'vtype'), vtype = p.vtype;else,   vtype = []; end
    if isfield(p, 'x0'),    x0 = p.x0;      else,   x0 = [];    end
    if isfield(p, 'xmax'),  xmax = p.xmax;  else,   xmax = [];  end
    if isfield(p, 'xmin'),  xmin = p.xmin;  else,   xmin = [];  end
    if isfield(p, 'u'),     u = p.u;        else,   u = [];     end
    if isfield(p, 'l'),     l = p.l;        else,   l = [];     end
    if isfield(p, 'A'),     A = p.A;        else,   A = [];     end
    if isfield(p, 'c'),     c = p.c;        else,   c = [];     end
    if isfield(p, 'H'),     H = p.H;        else,   H = [];     end
else                                %% individual args
    if nargin < 10
        opt = [];
        if nargin < 9
            vtype = [];
            if nargin < 8
                x0 = [];
                if nargin < 7
                    xmax = [];
                    if nargin < 6
                        xmin = [];
                    end
                end
            end
        end
    end
end

%% default options
if ~isempty(opt) && isfield(opt, 'alg') && ~isempty(opt.alg)
    alg = opt.alg;
    %% convert integer codes to string values
    if ~ischar(alg)
        switch alg
            case 0
                alg = 'DEFAULT';
            case 300
                alg = 'OT';
            case 500
                alg = 'CPLEX';
            case 600
                alg = 'MOSEK';
            case 700
                alg = 'GUROBI';
            otherwise
                error('miqps_master: %d is not a valid algorithm code', alg);
        end
    end
else
    alg = 'DEFAULT';
end
if ~isempty(opt) && isfield(opt, 'verbose') && ~isempty(opt.verbose)
    verbose = opt.verbose;
else
    verbose = 0;
end
if strcmp(alg, 'DEFAULT')
    if have_feature('gurobi')       %% use Gurobi by default, if available
        alg = 'GUROBI';
    elseif have_feature('cplex')    %% if not, then CPLEX, if available
        alg = 'CPLEX';
    elseif have_feature('mosek')    %% if not, then MOSEK, if available
        alg = 'MOSEK';
    elseif isempty(H) || ~any(any(H))   %% if not, and linear objective
        if have_feature('intlinprog')   %% then Optimization Tbx, if available
            alg = 'OT';
        elseif have_feature('glpk')     %% if not, and then GLPK, if available
            alg = 'GLPK';
        end
    else
        error('miqps_master: no solvers available - requires CPLEX, Gurobi, MOSEK, INTLINPROG or GLPK');
    end
end

%%----- call the appropriate solver  -----
switch alg
    case 'CPLEX'
        [x, f, eflag, output, lambda] = ...
            miqps_cplex(H, c, A, l, u, xmin, xmax, x0, vtype, opt);
    case 'GLPK'
        [x, f, eflag, output, lambda] = ...
            miqps_glpk(H, c, A, l, u, xmin, xmax, x0, vtype, opt);
    case 'GUROBI'
        [x, f, eflag, output, lambda] = ...
            miqps_gurobi(H, c, A, l, u, xmin, xmax, x0, vtype, opt);
    case 'MOSEK'
        [x, f, eflag, output, lambda] = ...
            miqps_mosek(H, c, A, l, u, xmin, xmax, x0, vtype, opt);
    case 'OT'
        [x, f, eflag, output, lambda] = ...
            miqps_ot(H, c, A, l, u, xmin, xmax, x0, vtype, opt);
    otherwise
        fcn = ['miqps_' lower(alg)];
        if exist([fcn '.m']) == 2
            [x, f, eflag, output, lambda] = ...
                feval(fcn, H, c, A, l, u, xmin, xmax, x0, vtype, opt);
        else
            error('miqps_master: ''%s'' is not a valid algorithm code', alg);
        end
end
if ~isfield(output, 'alg') || isempty(output.alg)
    output.alg = alg;
end
