function [x, f, eflag, output, lambda] = miqcqps_master(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, vtype, opt)
% miqcqps_master - Mixed Integer Quadratically Constrained Quadratic Program Solver wrapper function.
% ::
%
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       MIQCQPS_MASTER(H, C, Q, B, LQ, UQ, A, L, U, XMIN, XMAX, X0, VTYPE, OPT)
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = MIQCQPS_MASTER(PROBLEM)
%   A common wrapper function for various MIQCQP solvers.
%   Solves the following MIQCQP (mixed integer quadratically constrained quadratic 
%   programming) problem:
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
%       X(i) is integer, for i in I (integer variable constraints)
%       X(b) is binary, for b in B (binary variable constraints)
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
%       VTYPE : character string of length NX (number of elements in X),
%               or 1 (value applies to all variables in x),
%               allowed values are 'C' (continuous), 'B' (binary), or
%               'I' (integer), 'S' (semi-continuous), or 'N' (semi-integer).
%               (MOSEK allows only 'C', 'B', or 'I')
%       OPT : optional options structure with the following fields,
%           all of which are also optional (default values shown in
%           parentheses)
%           alg ('DEFAULT') - determines which solver to use
%               'DEFAULT' : automatic, first available of Gurobi, ...
%                           (Gurobi is only option currently implemented)
%               'GUROBI'  : Gurobi, requires Gurobi solver
%                           https://www.gurobi.com
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           fix_integer (0) - fix integer variables at value in x0, if true
%           relax_integer (0) - relax integer constraints, if true
%           grb_opt - options struct for GUROBI
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, vtype, opt
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
%   Note: Currently, LAMBDA is computed only for problems with all continuous
%   variables.
%
%   Calling syntax options:
%       [x, f, exitflag, output, lambda] = ...
%           miqcqps_master(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, vtype, opt)
%
%       x = miqcqps_master(H, c, Q, B, lq, uq)
%       x = miqcqps_master(H, c, Q, B, lq, uq, A, l, u)
%       x = miqcqps_master(H, c, Q, B, lq, uq, A, l, u, xmin, xmax)
%       x = miqcqps_master(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0)
%       x = miqcqps_master(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, vtype)
%       x = miqcqps_master(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, vtype, opt)
%       x = miqcqps_master(problem), where problem is a struct with fields:
%                       H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, vtype, opt
%                       all fields except 'c', are optional, and problem with
%                       linear costs must include constraints
%       x = miqcqps_master(...)
%       [x, f] = miqcqps_master(...)
%       [x, f, exitflag] = miqcqps_master(...)
%       [x, f, exitflag, output] = miqcqps_master(...)
%       [x, f, exitflag, output, lambda] = miqcqps_master(...)
%
%   Example: (problem from OPTI Toolbox, see:
%             https://jonathancurrie.github.io/OPTI/examples/problem-types/miqcqp/)
%
%       H = [1  0; 0  1];
%       c = [-2 -2]';
%       Q = {2*speye(2)};
%       B  = [0 -2];
%       lq = [];
%       uq = 1;
%       A = [-1 1; 1 3];
%       u = [2 5]';
%       l = [];
%       xmin = zeros(2,1);
%       xmax = Inf(2,1);
%       x0 = [];
%       vtype = 'IC';
%       opt = struct('verbose', 2);
%       [x, f, s, out, lam] = ...
%           miqcqps_master(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, vtype, opt);

%   MP-Opt-Model
%   Copyright (c) 2019-2026, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%   and Ray Zimmerman, PSERC Cornell
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
    if isfield(p, 'uq'),    uq = p.uq;      else,   uq = [];    end
    if isfield(p, 'lq'),    lq = p.lq;      else,   lq = [];    end
    if isfield(p, 'B'),     B = p.B;        else,   B = [];     end
    if isfield(p, 'Q'),     Q = p.Q;        else,   Q = {};     end
    if isfield(p, 'c'),     c = p.c;        else,   c = [];     end
    if isfield(p, 'H'),     H = p.H;        else,   H = [];     end
else                                %% individual args
    if nargin < 14
        opt = [];
        if nargin < 13
            vtype = [];
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
end

%% default options
if ~isempty(opt) && isfield(opt, 'alg') && ~isempty(opt.alg)
    alg = opt.alg;
    % convert integer codes to string values
    if ~ischar(alg)
        switch alg            
            case 700
                alg = 'GUROBI';
            otherwise
                error('miqcqps_master: %d is not a valid algorithm code. the only available solver with MIQCQP interface is Gurobi (700)', alg);
        end
    end
else
    alg = 'DEFAULT';
end
if strcmp(alg, 'DEFAULT')
    if have_feature('gurobi')
        alg = 'GUROBI';
    end
end

%% handle relax_integer and fix_integer options
done = false;
fix_integer_presolve = true;
if ~isempty(vtype) && (isfield(opt, 'relax_integer') && opt.relax_integer || ...
        isfield(opt, 'fix_integer') && opt.fix_integer)
    nx = length(x0);
    if length(vtype) == 1   %% expand if necessary
        vtype = char(vtype * ones(1, nx));
    end
    j = (vtype == 'B' | vtype == 'I')';
    if isfield(opt, 'fix_integer') && opt.fix_integer
        x0(j) = round(x0(j));
        if fix_integer_presolve
            x = x0;
            % 1) Reorganize quadratic constraints
            if isempty(Q)
                QQ = [];
                BB = B;
                Bxj = 0;
            else
                [nq, ncolQ] = size(Q);   %% Check dimension of quadratic terms of quadratic constraints
                if ~iscell(Q) || ncolQ ~= 1
                    error('miqcqps_master: Q must be column vector cell array.')
                end
                sizeQi  = cell2mat(cellfun(@(x)(size(x)), Q, 'UniformOutput', false));
                if any(sizeQi(:) - sizeQi(1))
                    error('miqcqps_master: All matrices Q{i}, i=1,...,%d must be square of the same size.', nq)
                end

                QQ = cellfun(@(w)(w(~j,~j)), Q, 'UniformOutput', false);
                BB = cell2mat(cellfun(@(w)( x(j)' * (1/2*(w(~j,j)'+w(j,~j))) ), Q, 'UniformOutput', false));
                Bxj = cell2mat(cellfun(@(w)( 1/2 * x(j)' * w(j,j) * x(j) ), Q, 'UniformOutput', false));
            end
            % 2) Reorganize linear constraints
            Axj = A(:,j) * x(j);            
            
            % 3) Reorganize objective function
            cc = c(~j);
            if isempty(H)
                HH = [];
            else
                HH = H(~j, ~j);
                cc = cc + 1/2*(H(~j,j)+H(j,~j)') * x(j);
            end
            [x(~j), f, eflag, output, lambda] = ...
                qcqps_master(HH, cc, QQ, BB, lq-Bxj, uq-Bxj, A(:, ~j), l-Axj, u-Axj, ...
                    xmin(~j), xmax(~j), x(~j), opt);
            f = f + c(j)' * x(j);
            if ~isempty(H)
                f = f + 1/2 * x(j)' * H(j,j) * x(j);
            end
            mu_lower = zeros(size(x));
            mu_upper = zeros(size(x));
            mu_lower(~j) = lambda.lower;
            mu_upper(~j) = lambda.upper;
            lambda.lower = mu_lower;
            lambda.upper = mu_upper;
            done = true;
        else
            %% fix integer variables
            if ~isempty(j)
                %% expand if necessary
                if length(xmin) == 1
                    xmin = xmin * ones(nx, 1);
                elseif isempty(xmin)
                    xmin = -Inf(nx, 1);
                end
                if length(xmax) == 1
                    xmax = xmax * ones(nx, 1);
                elseif isempty(xmax)
                    xmax = Inf(nx, 1);
                end
                xmin(j) = x0(j);
                xmax(j) = x0(j);
            end
        end
    end
    if ~done
        vtype(j) = 'C';
    end
end

miqp = 1;
if any(cell2mat(cellfun(@(x)(any(x(:))), Q, 'UniformOutput', false)))
    miqp = 0;
end

%%----- call the appropriate solver  -----
if ~done
    if miqp && ~strcmp(alg, 'GUROBI')
        [x, f, eflag, output, lambda] = ...
            miqps_master(H, c, A, l, u, xmin, xmax, x0, vtype, opt);
    else
        switch alg
            case 'GUROBI'
                [x, f, eflag, output, lambda] = ...
                    miqcqps_gurobi(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, vtype, opt);
            otherwise
                error('miqcqps_master: currently, the only available solver with MIQCQP interface is Gurobi')
        end
    end
end
if ~isfield(output, 'alg') || isempty(output.alg)
    output.alg = alg;
end
