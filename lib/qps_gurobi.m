function [x, f, eflag, output, lambda] = qps_gurobi(H, c, A, l, u, xmin, xmax, x0, opt)
%QPS_GUROBI  Quadratic Program Solver based on GUROBI.
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       QPS_GUROBI(H, C, A, L, U, XMIN, XMAX, X0, OPT)
%   A wrapper function providing a MATPOWER standardized interface for using
%   GUROBI_MEX to solve the following QP (quadratic programming)
%   problem:
%
%       min 1/2 X'*H*X + C'*X
%        X
%
%   subject to
%
%       L <= A*X <= U       (linear constraints)
%       XMIN <= X <= XMAX   (variable bounds)
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
%       OPT : optional options structure with the following fields,
%           all of which are also optional (default values shown in
%           parentheses)
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           grb_opt - options struct for GUROBI, value in
%               verbose overrides these options
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: H, c, A, l, u, xmin, xmax, x0, opt
%
%   Outputs:
%       X : solution vector
%       F : final objective function value
%       EXITFLAG : GUROBI_MEX exit flag
%           1 = converged
%           0 or negative values = negative of GUROBI_MEX exit flag
%           (see GUROBI_MEX documentation for details)
%       OUTPUT : GUROBI_MEX output struct
%           (see GUROBI_MEX documentation for details)
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
%           qps_gurobi(H, c, A, l, u, xmin, xmax, x0, opt)
%
%       x = qps_gurobi(H, c, A, l, u)
%       x = qps_gurobi(H, c, A, l, u, xmin, xmax)
%       x = qps_gurobi(H, c, A, l, u, xmin, xmax, x0)
%       x = qps_gurobi(H, c, A, l, u, xmin, xmax, x0, opt)
%       x = qps_gurobi(problem), where problem is a struct with fields:
%                       H, c, A, l, u, xmin, xmax, x0, opt
%                       all fields except 'c', 'A' and 'l' or 'u' are optional
%       x = qps_gurobi(...)
%       [x, f] = qps_gurobi(...)
%       [x, f, exitflag] = qps_gurobi(...)
%       [x, f, exitflag, output] = qps_gurobi(...)
%       [x, f, exitflag, output, lambda] = qps_gurobi(...)
%
%
%   Example: (problem from from http://www.jmu.edu/docs/sasdoc/sashtml/iml/chap8/sect12.htm)
%       H = [   1003.1  4.3     6.3     5.9;
%               4.3     2.2     2.1     3.9;
%               6.3     2.1     3.5     4.8;
%               5.9     3.9     4.8     10  ];
%       c = zeros(4,1);
%       A = [   1       1       1       1;
%               0.17    0.11    0.10    0.18    ];
%       l = [1; 0.10];
%       u = [1; Inf];
%       xmin = zeros(4,1);
%       x0 = [1; 0; 0; 1];
%       opt = struct('verbose', 2);
%       [x, f, s, out, lambda] = qps_gurobi(H, c, A, l, u, xmin, [], x0, opt);
%
%   See also GUROBI_MEX.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2010-2011 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

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
    if isfield(p, 'c'),     c = p.c;        else,   c = [];     end
    if isfield(p, 'H'),     H = p.H;        else,   H = [];     end
else                                %% individual args
    if nargin < 9
        opt = [];
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

%% define nx, set default values for missing optional inputs
if isempty(H) || ~any(any(H))
    if isempty(A) && isempty(xmin) && isempty(xmax)
        error('qps_gurobi: LP problem must include constraints or variable bounds');
    else
        if ~isempty(A)
            nx = size(A, 2);
        elseif ~isempty(xmin)
            nx = length(xmin);
        else    % if ~isempty(xmax)
            nx = length(xmax);
        end
    end
else
    nx = size(H, 1);
end
if isempty(c)
    c = zeros(nx, 1);
end
if  ~isempty(A) && (isempty(l) || all(l == -Inf)) && ...
                   (isempty(u) || all(u == Inf))
    A = sparse(0,nx);           %% no limits => no linear constraints
end
nA = size(A, 1);                %% number of original linear constraints
if isempty(u)                   %% By default, linear inequalities are ...
    u = Inf * ones(nA, 1);      %% ... unbounded above and ...
end
if isempty(l)
    l = -Inf * ones(nA, 1);     %% ... unbounded below.
end
if isempty(xmin)                %% By default, optimization variables are ...
    xmin = -Inf * ones(nx, 1);  %% ... unbounded below and ...
end
if isempty(xmax)
    xmax = Inf * ones(nx, 1);   %% ... unbounded above.
end
if isempty(x0)
    x0 = zeros(nx, 1);
end

%% default options
if ~isempty(opt) && isfield(opt, 'verbose') && ~isempty(opt.verbose)
    verbose = opt.verbose;
else
    verbose = 0;
end
% if ~isempty(opt) && isfield(opt, 'max_it') && ~isempty(opt.max_it)
%     max_it = opt.max_it;
% else
%     max_it = 0;
% end

%% set up options struct for Gurobi
if ~isempty(opt) && isfield(opt, 'grb_opt') && ~isempty(opt.grb_opt)
    g_opt = gurobi_options(opt.grb_opt);
else
    g_opt = gurobi_options;
end
g_opt.Display = min(verbose, 3);
if verbose
    g_opt.DisplayInterval = 1;
else
    g_opt.DisplayInterval = Inf;
end

if ~issparse(A)
    A = sparse(A);
end

%% split up linear constraints
ieq = find( abs(u-l) <= eps );          %% equality
igt = find( u >=  1e10 & l > -1e10 );   %% greater than, unbounded above
ilt = find( l <= -1e10 & u <  1e10 );   %% less than, unbounded below
ibx = find( (abs(u-l) > eps) & (u < 1e10) & (l > -1e10) );

%% grab some dimensions
nlt = length(ilt);      %% number of upper bounded linear inequalities
ngt = length(igt);      %% number of lower bounded linear inequalities
nbx = length(ibx);      %% number of doubly bounded linear inequalities
neq = length(ieq);      %% number of equalities
niq = nlt+ngt+2*nbx;    %% number of inequalities

AA  = [ A(ieq, :); A(ilt, :); -A(igt, :); A(ibx, :); -A(ibx, :) ];
bb  = [ u(ieq);    u(ilt);    -l(igt);    u(ibx);    -l(ibx)    ];
contypes = char([ double('=')*ones(1,neq) double('<')*ones(1,niq) ]);

%% call the solver
if isempty(H) || ~any(any(H))
    lpqp = 'LP';
else
    lpqp = 'QP';
    [rr, cc, vv] = find(H);
    g_opt.QP.qrow = int32(rr' - 1);
    g_opt.QP.qcol = int32(cc' - 1);
    g_opt.QP.qval = 0.5 * vv';
end
if verbose
    methods = {
        'primal simplex',
        'dual simplex',
        'interior point',
        'concurrent',
        'deterministic concurrent'
    };
    fprintf('Gurobi Version %s -- %s %s solver\n', ...
        '<unknown>', methods{g_opt.Method+1}, lpqp);
end
[x, f, eflag, output, lambda, rc] = ...
    gurobi_mex(c', 1, AA, bb, contypes, xmin, xmax, 'C', g_opt);
output.flag = eflag;
if eflag == 2
    eflag = 1;          %% optimal solution found
else
    eflag = -eflag;     %% failed somehow
end

%% check for empty results (in case optimization failed)
if isempty(x)
    x = NaN(nx, 1);
    lam.lower   = NaN(nx, 1);
    lam.upper   = NaN(nx, 1);
else
    lam.lower   = zeros(nx, 1);
    lam.upper   = zeros(nx, 1);
end
if isempty(f)
    f = NaN;
end
if isempty(lambda)
    lambda      = NaN(length(bb), 1);
end

kl = find(rc > 0);   %% lower bound binding
ku = find(rc < 0);   %% upper bound binding
lam.lower(kl)   =  rc(kl);
lam.upper(ku)   = -rc(ku);
lam.eqlin   = lambda(1:neq);
lam.ineqlin = lambda(neq+(1:niq));
mu_l        = zeros(nA, 1);
mu_u        = zeros(nA, 1);

%% repackage lambdas
kl = find(lam.eqlin > 0);   %% lower bound binding
ku = find(lam.eqlin < 0);   %% upper bound binding

mu_l(ieq(kl)) = lam.eqlin(kl);
mu_l(igt) = -lam.ineqlin(nlt+(1:ngt));
mu_l(ibx) = -lam.ineqlin(nlt+ngt+nbx+(1:nbx));

mu_u(ieq(ku)) = -lam.eqlin(ku);
mu_u(ilt) = -lam.ineqlin(1:nlt);
mu_u(ibx) = -lam.ineqlin(nlt+ngt+(1:nbx));

lambda = struct( ...
    'mu_l', mu_l, ...
    'mu_u', mu_u, ...
    'lower', lam.lower, ...
    'upper', lam.upper ...
);
