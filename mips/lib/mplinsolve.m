function [x, info] = mplinsolve(A, b, solver, opt)
%MPLINSOLVE  Solves A * x = b using specified solver
%   X = MPLINSOLVE(A, B)
%   X = MPLINSOLVE(A, B, SOLVER)
%   X = MPLINSOLVE(A, B, SOLVER, OPT)
%   [X, INFO] = MPLINSOLVE(...)
%
%   Solves the linear system of equations A * x = b, using the selected
%   solver.
%
%   Inputs:
%       A      : sparse matrix
%       B      : full vector or matrix
%       SOLVER : ('') selected linear system solver
%           ''  - use default solver, currently this is
%                   always the built-in backslash operator
%           '\' - built-in backslash operator
%           'LU' - use explicit default LU decomposition and back substitution
%           'LU_AMD' - use LU decomposition with approximate minimum degree
%                     (AMD) reordering
%           'LU_GP' - use Gilbert-Peierls algorithm for LU Decomposition with
%                  (AMD) reordering (e.g. ideal for power flow Jacobian)
%           'PARDISO' - PARDISO
%       OPT    : struct of options, with the following fields
%                (currently used only by PARDISO, default shown in parens,
%                 see PARDISO documentation for details)
%           verbose (0) - true or false
%           mtype (11)  - matrix type (default is real and nonsymmetric)
%           solver (0)  - solver method (default is sparse direct)
%           iparm ([])  - n x 2 matrix of integer parameters
%               1st, 2nd columns are index, value of parameter respectively
%           dparm ([])  - n x 2 matrix of double parameters
%               1st, 2nd columns are index, value of parameter respectively

%   Note: SOLVER can also take value 'LU_GP0' or 'LU_AMD0', both of which
%         just use permutation matrices instead of permutation vectors,
%         which is slightly slower.
%         Included for compatibility with Matlab < 7.3.

%   MIPS
%   Copyright (c) 2015-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MIPS.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mips for more info.

if nargin < 4
    opt = [];
    if nargin < 3
        solver = ''
    end
end

info = [];

switch upper(solver)
    case {'', '\'}      %% use built-in \ operator
        x = A \ b;
    case {'LU'}         %% explicit LU (is this what's done internally by '\'?)
        [L, U, p, q] = lu(A, 'vector');
        x = zeros(size(A, 1), 1);
        x(q) = U \ ( L \ b(p));
    case {'LU_GP', 'LU_AMD', 'LU_GP0', 'LU_AMD0'}   %% LU with AMD reordering
        q = amd(A);                 %% permutation vector for AMD reordering
        n = size(A, 1);
        if have_lu_vec() && solver(end) ~= '0'
            if solver(4) == 'G'     %% Gilbert-Peierls LU
                %% set pivot threshold to 1 (partial pivoting)
                %% to ensure that Gilbert-Peierls is used
                %% (thanks to Jose Luis Marin, see also T. Davis)
                [L, U, p] = lu(A(q,q), 1.0, 'vector');
            else
                [L, U, p] = lu(A(q,q), 'vector');
            end
            x = zeros(n, 1);
            x(q) = U \ ( L \ b(q(p)) );
        else
            Q = sparse(q, 1:n, 1);  %% permutation matrix for AMD reordering
            if solver(4) == 'G'     %% Gilbert-Peierls LU
                %% set pivot threshold to 1 (partial pivoting)
                %% to ensure that Gilbert-Peierls is used
                %% (thanks to Jose Luis Marin, see also T. Davis)
                [L, U, P] = lu(Q'*A*Q, 1.0);
            else
                [L, U, P] = lu(Q'*A*Q);
            end
            x = Q * ( U \ ( L \ (P * Q' * b)) );
        end
    case {'PARDISO'}
        %% get number of threads from OpenMP env variable
        persistent num_threads;
        if isempty(num_threads)
            num_threads = str2num(getenv('OMP_NUM_THREADS'));
            if ~num_threads
                num_threads = 1;
            end
        end

        %% set default options
        verbose = false;
        mtype = 11;
        solver = 0;

        %% override if provided via opt
        if ~isempty(opt)
            if isfield(opt, 'verbose')
                verbose = opt.verbose;
            end
            if isfield(opt, 'mtype')
                mtype = opt.mtype;
            end
            if isfield(opt, 'solver')
                solver = opt.solver;
            end
        end

        %% begin setup and solve
        info = pardisoinit(mtype, solver);
        info.iparm(3) = num_threads;
        if ~isempty(opt)
            if isfield(opt, 'iparm') && ~isempty(opt.iparm)
                info.iparm(opt.iparm(:, 1)) = opt.iparm(:, 2);
            end
            if isfield(opt, 'dparm') && ~isempty(opt.dparm)
                info.iparm(opt.dparm(:, 1)) = opt.dparm(:, 2);
            end
        end
        info = pardisoreorder(A, info, verbose);
        info = pardisofactor(A, info, verbose);
        [x, info] = pardisosolve(A, b, info, verbose);
        pardisofree(info);
    otherwise
        warning('mplinsolve: ''%s'' is not a valid value for SOLVER, using default.', solver);
        x = A \ b;
end


function TorF = have_lu_vec()
% Checks whether or not LU supports lu(..., 'vector') syntax
persistent lu_vec;      %% cache the result for performance reasons
if isempty(lu_vec)
    lu_vec = 1;         %% assume it does, unless this is Matlab ver < 7.3
    v = ver('matlab');
    if length(v) > 1
        warning('The built-in VER command is behaving strangely, probably as a result of installing a 3rd party toolbox in a directory named ''matlab'' on your path. Check each element of the output of ver(''matlab'') to find the offending toolbox, then move the toolbox to a more appropriately named directory.');
        v = v(1);
    end
    if ~isempty(v) && isfield(v, 'Version') && ~isempty(v.Version)
        vstr = v.Version;
        if ~isempty(vstr) && vstr2num_(vstr) < 7.003
            lu_vec = 0;
        end
    end
end
TorF = lu_vec;


function num = vstr2num_(vstr)
% Converts version string to numerical value suitable for < or > comparisons
% E.g. '3.11.4' -->  3.011004
pat = '\.?(\d+)';
[s,e,tE,m,t] = regexp(vstr, pat);
b = 1;
num = 0;
for k = 1:length(t)
    num = num + b * str2num(t{k}{1});
    b = b / 1000;
end
