function [x, f, eflag, output, lambda] = qps_matpower(H, c, A, l, u, xmin, xmax, x0, opt)
%QPS_MATPOWER  Quadratic Program Solver for MATPOWER.
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       QPS_MATPOWER(H, C, A, L, U, XMIN, XMAX, X0, OPT)
%   A common wrapper function for various QP solvers. 
%   Solves the following QP (quadratic programming) problem:
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
%           alg ('DEFAULT') : determines which solver to use, can be either
%                   a (new-style) string or an (old-style) numerical alg code
%               'DEFAULT' : (or 0) automatic, first available of CPLEX,
%                       Gurobi, MOSEK, BPMPD, Opt Tbx, MIPS
%               'BPMPD'   : (or 100) BPMPD_MEX
%               'MIPS'    : (or 200) MIPS, MATLAB Interior Point Solver
%                        pure MATLAB implementation of a primal-dual
%                        interior point method, if mips_opt.step_control = 1
%                        (or alg=250) it uses MIPS-sc, a step controlled
%                        variant of MIPS
%               'OT'      : (or 300) Optimization Toolbox, QUADPROG or LINPROG
%               'IPOPT'   : (or 400) IPOPT
%               'CPLEX'   : (or 500) CPLEX
%               'MOSEK'   : (or 600) MOSEK
%               'GUROBI'  : (or 700) Gurobi
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           max_it (0) - maximum number of iterations allowed
%               0 = use algorithm default
%           bp_opt    - options vector for BP
%           cplex_opt - options struct for CPLEX
%           grb_opt   - options struct for GBUROBI_MEX
%           ipopt_opt - options struct for IPOPT
%           mips_opt  - options struct for QPS_MIPS
%           mosek_opt - options struct for MOSEK
%           ot_opt    - options struct for QUADPROG/LINPROG
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: H, c, A, l, u, xmin, xmax, x0, opt
%
%   Outputs:
%       X : solution vector
%       F : final objective function value
%       EXITFLAG : exit flag
%           1 = converged
%           0 or negative values = algorithm specific failure codes
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
%           qps_matpower(H, c, A, l, u, xmin, xmax, x0, opt)
%
%       x = qps_matpower(H, c, A, l, u)
%       x = qps_matpower(H, c, A, l, u, xmin, xmax)
%       x = qps_matpower(H, c, A, l, u, xmin, xmax, x0)
%       x = qps_matpower(H, c, A, l, u, xmin, xmax, x0, opt)
%       x = qps_matpower(problem), where problem is a struct with fields:
%                       H, c, A, l, u, xmin, xmax, x0, opt
%                       all fields except 'c', 'A' and 'l' or 'u' are optional
%       x = qps_matpower(...)
%       [x, f] = qps_matpower(...)
%       [x, f, exitflag] = qps_matpower(...)
%       [x, f, exitflag, output] = qps_matpower(...)
%       [x, f, exitflag, output, lambda] = qps_matpower(...)
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
%       [x, f, s, out, lambda] = qps_matpower(H, c, A, l, u, xmin, [], x0, opt);

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2010-2013 by Power System Engineering Research Center (PSERC)
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

%% default options
if ~isempty(opt) && isfield(opt, 'alg') && ~isempty(opt.alg)
    alg = opt.alg;
    %% convert integer codes to string values
    if ~ischar(alg)
        switch alg
            case 0
                alg = 'DEFAULT';
            case 100
                alg = 'BPMPD';
            case 200
                alg = 'MIPS';
                opt.mips_opt.step_control = 0;
            case 250
                alg = 'MIPS';
                opt.mips_opt.step_control = 1;
            case 300
                alg = 'OT';
            case 400
                alg = 'IPOPT';
            case 500
                alg = 'CPLEX';
            case 600
                alg = 'MOSEK';
            case 700
                alg = 'GUROBI';
            otherwise
                error('qps_matpower: %d is not a valid algorithm code', alg);
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
    if have_fcn('cplex')        %% use CPLEX by default, if available
        alg = 'CPLEX';
    elseif have_fcn('gurobi')   %% if not, then Gurobi, if available
        alg = 'GUROBI';
    elseif have_fcn('mosek')    %% if not, then MOSEK, if available
        alg = 'MOSEK';
    elseif have_fcn('bpmpd')    %% if not, then BPMPD_MEX, if available
        alg = 'BPMPD';
    elseif have_fcn('quadprog') %% if not, then Optimization Tbx, if available
        alg = 'OT';
    else                        %% otherwise MIPS
        alg = 'MIPS';
    end
end

%%----- call the appropriate solver  -----
switch alg
    case 'BPMPD'                    %% use BPMPD_MEX
        [x, f, eflag, output, lambda] = ...
            qps_bpmpd(H, c, A, l, u, xmin, xmax, x0, opt);

        if eflag == -99
            if verbose
                fprintf('         Retrying with QPS_MIPS solver ...\n\n');
            end
            %% save (incorrect) solution from BPMPD
            bpmpd = struct('x', x, 'f', f, 'eflag', eflag, ...
                            'output', output, 'lambda', lambda);
            opt.alg = 'MIPS';
            [x, f, eflag, output, lambda] = ...
                qps_matpower(H, c, A, l, u, xmin, xmax, x0, opt);
            output.bpmpd = bpmpd;
        end
    case 'MIPS'
        %% set up options
        if ~isempty(opt) && isfield(opt, 'mips_opt') && ~isempty(opt.mips_opt)
            mips_opt = opt.mips_opt;
        else
            mips_opt = [];
        end
        if ~isempty(opt) && isfield(opt, 'max_it') && ~isempty(opt.max_it)
            mips_opt.max_it = opt.max_it;
        end
        mips_opt.verbose = verbose;
        
        %% call solver
        [x, f, eflag, output, lambda] = ...
            qps_mips(H, c, A, l, u, xmin, xmax, x0, mips_opt);
    case 'OT'                    %% use QUADPROG or LINPROG from Opt Tbx ver 2.x+
        [x, f, eflag, output, lambda] = ...
            qps_ot(H, c, A, l, u, xmin, xmax, x0, opt);
    case 'IPOPT'
        [x, f, eflag, output, lambda] = ...
            qps_ipopt(H, c, A, l, u, xmin, xmax, x0, opt);
    case 'CPLEX'
        [x, f, eflag, output, lambda] = ...
            qps_cplex(H, c, A, l, u, xmin, xmax, x0, opt);
    case 'MOSEK'
        [x, f, eflag, output, lambda] = ...
            qps_mosek(H, c, A, l, u, xmin, xmax, x0, opt);
    case 'GUROBI'
        [x, f, eflag, output, lambda] = ...
            qps_gurobi(H, c, A, l, u, xmin, xmax, x0, opt);
    otherwise
        error('qps_matpower: %d is not a valid algorithm code', alg);
end
if ~isfield(output, 'alg') || isempty(output.alg)
    output.alg = alg;
end
