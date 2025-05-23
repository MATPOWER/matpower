function [x, f, eflag, output, lambda] = qcqps_master(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, opt)
% qcqps_master - Quadratically Constrained Quadratic Program Solver wrapper function.
% ::
%
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       QCQPS_MASTER(H, C, Q, B, LQ, UQ, A, L, U, XMIN, XMAX, X0, OPT)
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = QCQPS_MASTER(PROBLEM)
%   A common wrapper function for various QCQP solvers.
%   Solves the following QCQP (quadratically constrained quadratic programming)
%   problem:
%
%       min 1/2 X'*H*X + C'*X
%        X
%
%   subject to
%
%       LQ(i) <= 1/2 X'*Q{i}*X + B(i,:)*X <= UQ(i), i = 1,2,...,nq
%                           (quadratic constraints)
%       L <= A*X <= U       (linear constraints)
%       XMIN <= X <= XMAX   (variable bounds)
%
%   Inputs (all optional except H, C, Q, B, LQ, and UQ):
%       H : matrix (possibly sparse) of quadratic cost coefficients
%       C : vector of linear cost coefficients
%       Q : nq x 1 cell array of sparse quadratic matrices for quadratic constraints
%       B : matrix (possibly sparse) of linear term of quadratic constraints
%       LQ, UQ: define the lower an upper bounds on the quadratic constraints
%       A, L, U : define the optional linear constraints. Default
%           values for the elements of L and U are -Inf and Inf,
%           respectively.
%       XMIN, XMAX : optional lower and upper bounds on the
%           X variables, defaults are -Inf and Inf, respectively.
%       X0 : optional starting value of optimization vector X
%       OPT : optional options structure with the following fields,
%           all of which are also optional (default values shown in
%           parentheses)
%           alg ('DEFAULT') : determines which solver to use
%               'DEFAULT' :  automatic, first available of IPOPT, Artelys
%                   Knitro, FMINCON, MIPS
%               'FMINCON' : FMINCON, MATLAB Optimization Toolbox
%               'GUROBI'  : Gurobi, requires Gurobi solver
%                           https://www.gurobi.com
%               'IPOPT'   : IPOPT, requires MEX interface to IPOPT
%                           solver, https://github.com/coin-or/Ipopt
%               'KNITRO'  : Artelys Knitro, requires Artelys Knitro solver
%                           https://www.artelys.com/solvers/knitro/
%               'KNITRO_NLP' : Artelys Knitro, via NLPS_MASTER, requires
%                           Artelys Knitro solver
%                           https://www.artelys.com/solvers/knitro/
%               'MIPS'    : MIPS, MATPOWER Interior Point Solver
%                        pure MATLAB implementation of a primal-dual
%                        interior point method, if mips_opt.step_control = 1
%                        it uses MIPS-sc, a step controlled variant of MIPS
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           grb_opt     - options struct for GUROBI
%           knitro_opt  - options struct for Artelys Knitro
%           ipopt_opt   - options struct for IPOPT
%           mips_opt    - options struct for MIPS
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, opt
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
%           mu_lq - lower (left-hand) limit on quadratic constraints
%           mu_uq - upper (right-hand) limit on quadratic constraints
%           lower - lower bound on optimization variables
%           upper - upper bound on optimization variables
%
%   Calling syntax options:
%       [x, f, exitflag, output, lambda] = ...
%           qcqps_master(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, opt)
%
%       x = qcqps_master(H, c, Q, B, lq, uq)
%       x = qcqps_master(H, c, Q, B, lq, uq, A, l, u)
%       x = qcqps_master(H, c, Q, B, lq, uq, A, l, u, xmin, xmax)
%       x = qcqps_master(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0)
%       x = qcqps_master(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, opt)
%       x = qcqps_master(problem), where problem is a struct with fields:
%                       H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, opt
%                       all fields except 'c', are optional, and problem with
%                       linear costs must include constraints
%       x = qcqps_master(...)
%       [x, f] = qcqps_master(...)
%       [x, f, exitflag] = qcqps_master(...)
%       [x, f, exitflag, output] = qcqps_master(...)
%       [x, f, exitflag, output, lambda] = qcqps_master(...)
%
%   Example: (problem from https://docs.gurobi.com/projects/examples/en/current/examples/matlab/qcp.html)
%       H = [];
%       c = [-1;0;0];
%       Q = { sparse([2 0 0; 0 2 0; 0 0 -2]), ...
%             sparse([2 0 0; 0 0 -2; 0 -2 0]) };
%       B = zeros(2,3);
%       lq = [-Inf;-Inf];
%       uq = [0; 0];
%       A = [1 1 1];
%       l = 1;
%       u = 1;
%       xmin = zeros(3,1);
%       xmax = Inf(3,1);
%       x0 = zeros(3,1);
%       opt = struct('verbose', 2);
%       [x, f, s, out, lambda] = ...
%           qcqps_master(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, opt);

%   MP-Opt-Model
%   Copyright (c) 2019-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model..
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%%----- input argument handling  -----
%% gather inputs
if nargin == 1 && isstruct(H)       %% problem struct
    p = H;
    if isfield(p, 'opt'),   opt = p.opt;    else,   opt = [];   end
    if isfield(p, 'x0'),    x0 = p.x0;      else,   x0 = [];    end
    if isfield(p, 'xmax'),  xmax = p.xmax;  else,   xmax = [];  end
    if isfield(p, 'xmin'),  xmin = p.xmin;  else,   xmin = [];  end
    if isfield(p, 'u'),     u = p.u;        else,   u = [];     end
    if isfield(p, 'l'),     l = p.l;        else,   l = [];     end
    if isfield(p, 'A'),     A = p.A;        else,   A = [];     end
    if isfield(p, 'uq'),    uq = p.uq;      else,   uq = [];    end
    if isfield(p, 'lq'),    lq = p.lq;      else,   lq = [];    end
    if isfield(p, 'B'),     B = p.B;        else,   B = [];     end
    if isfield(p, 'Q'),     Q = p.Q;        else,   Q = {};     end
    if isfield(p, 'c'),     c = p.c;        else,   c = [];     end
    if isfield(p, 'H'),     H = p.H;        else,   H = [];     end
else                                %% individual args
    if nargin < 13
        opt = [];
        if nargin < 12
            x0 = [];
            if nargin < 11
                xmax = [];
                if nargin < 10
                    xmin = [];
                    if nargin < 7
                        A = [];
                        l = [];
                        u = [];
                    end
                end
            end
        end
    end
end

%% default options
if ~isempty(opt) && isfield(opt, 'alg') && ~isempty(opt.alg)
    alg = opt.alg;
else
    alg = 'DEFAULT';
end
if strcmp(alg, 'DEFAULT')
    if have_feature('ipopt')
        alg = 'IPOPT';
    elseif have_feature('knitromatlab') %% if not, then Artelys Knitro, if available
        alg = 'KNITRO';
    elseif have_feature('fmincon') && have_feature('matlab')    %% if not, then Opt Tbx, if available in MATLAB
        alg = 'FMINCON';
    else                                %% otherwise MIPS
        alg = 'MIPS';
    end
end

%%----- call the appropriate solver  -----
switch alg
    case 'GUROBI'
        [x, f, eflag, output, lambda] = ...
            qcqps_gurobi(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, opt);
    case 'KNITRO'
        [x, f, eflag, output, lambda] = ...
            qcqps_knitro(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, opt);
%    case {'MOSEK'}
%        [x, f, eflag, output, lambda] = ...
%            qcqps_mosek(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, opt);
    otherwise
        %% strip trailing '_NLP' from alg, if present
        opt.alg = regexprep(alg, '_NLP$', '');

        [x, f, eflag, output, lambda] = ...
            qcqps_nlps(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, opt);
end
if ~isfield(output, 'alg') || isempty(output.alg)
    output.alg = alg;
end
