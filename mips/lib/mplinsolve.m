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

%   MIPS
%   Copyright (c) 2015-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MIPS.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 4
    opt = [];
    if nargin < 3
        solver = ''
    end
end

info = [];

switch upper(solver)
    case {'', '\'}
        x = A \ b;
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
