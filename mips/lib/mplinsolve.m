function [x, f] = mplinsolve(A, b, solver, opt)
% mplinsolve - Solve A * x = b using specified solver.
% ::
%
%   x = mplinsolve(A, b)
%   x = mplinsolve(A, b, solver)
%   x = mplinsolve(A, b, solver, opt)
%   x = mplinsolve(factors, b, ...)
%   factors = mplinsolve(A)
%   factors = mplinsolve(A, [], ...)
%   [x, factors] = mplinsolve(A, b, 'lu' ...)
%
% Solves the linear system of equations ``A * x == b`` for ``x``
% using the selected solver.
%
% Inputs:
%     A (double): sparse matrix
%     factors (struct) : contains (e.g. LU) factors of ``A`` from previous
%       call to mplinsolve, with a ``type`` field to identify which of the
%       following sets of additional fields are included:
%
%       ========  =================  ========================================
%       ``type``  add'l fields       description
%       ========  =================  ========================================
%          1      ``L, U, p, qa``    3 output Gilbert-Peierls, permutation vectors
%          2      ``L, U, P, Qa``    3 output Gilbert-Peierls, permutation matrices
%          3      ``L, U, p, q``     4 output UMFPACK, permutation vectors
%          4      ``L, U, P, Q``     4 output UMFPACK, permutation matrices
%          5      ``L, U, p, q, R``  same as 3, with row scaling
%          6      ``L, U, P, Q, R``  same as 4, with row scaling
%       ========  =================  ========================================
%     b (double) : RHS vector (or full matrix)
%     solver (char array) : *(optional, default* ``''`` *)* selected linear
%       system solver
%
%         - ``''``  -- use default solver, currently this is
%           always the built-in backslash operator
%         - ``'\'`` -- built-in backslash operator
%         - ``'LU'`` -- use explicit LU decomposition and back substitution
%
%           The following are also provided as short-cuts, with less
%           overhead and thus better performance on small systems, as
%           alternatives for ``'LU'`` with the following options:
%
%             ===========  ===============  ==============  ================
%             ``solver``   ``opt.lu.nout``  ``opt.lu.vec``  ``opt.lu.thresh``
%             ===========  ===============  ==============  ================
%             ``'LU3'``        3               1               1.0
%             ``'LU3a'``       3               1
%             ``'LU4'``        4               1
%             ``'LU5'``        5               1
%             ``'LU3m'``       3               0               1.0
%             ``'LU3am'``      3               0
%             ``'LU4m'``       4               0
%             ``'LU5m'``       5               0
%             ===========  ===============  ==============  ================
%         - ``'PARDISO'`` -- PARDISO
%     opt (struct) : options for certain solvers (e.g. :func:`lu` and PARDISO),
%       with fields:
%
%         - ``lu`` : struct of options to determine form of call to :func:`lu`
%           solver, with the following possible fields *(default value in
%           parens)*:
%
%             - ``nout`` *(4)* - number of output args for call to :func:`lu`,
%               UMFPACK is used for 4 or 5 output args, Gilbert-Peierls
%               algorithm with AMD ordering for 3 output args.
%             - ``vec`` *(1)*  - use permutation vectors instead of matrices
%               (permutation matrices used by default for MATLAB < 7.3)
%             - ``thresh``   - pivot threshold, see ``help lu`` for details
%         - ``pardiso`` : struct of PARDISO options *(default shown in parens)*,
%           see PARDISO documentation for details
%
%             - ``verbose`` *(0)* - true or false
%             - ``mtype`` *(11, i.e. real and nonsymmetric)*  - matrix type
%             - ``solver`` *(0, i.e. sparse direct)*  - solver method
%             - ``iparm`` *([])*  - :math:`n \by 2` matrix of integer parameters
%               1st, 2nd columns are index, value of parameter respectively
%             - ``dparm`` *([])*  - :math:`n \by 2` matrix of double parameters
%               1st, 2nd columns are index, value of parameter respectively
%         - ``tr`` *(0)* : if true, solve transposed system ``A' * x = b``
%
% Outputs:
%     x : solution vector (or matrix) satisfying ``A * x == b``
%     factors : see description under Inputs

%   MIPS
%   Copyright (c) 2015-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MIPS.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mips for more info.

if nargin < 4
    opt = [];
    if nargin < 3
        solver = '';
    end
end

%% solve transpose?
if isfield(opt, 'tr')
    tr = opt.tr;
else
    tr = 0;
end
have_f = isstruct(A);
have_b = nargin >= 2 && ~isempty(b);
is_LU = length(solver) >= 2 && all(solver(1:2) == 'LU');
if ~is_LU && ~have_b
    solver = 'LU';
    is_LU = 1;
end
is_LU = is_LU || have_f;

%% prepare LU factors, if necessary
if have_f
    f = A;      %% factors provided
elseif is_LU    %% perform factorization if it's LU & we have A (not factors)
    switch solver
        case 'LU3'      %% 3 output LU: Gilbert-Peierls alg, perm vec, 1.0 piv thresh
            q = amd(A);     %% permutation vector for AMD reordering
            if issparse(A)
                [L, U, p] = lu(A(q,q), 1.0, 'vector');
            else
                [L, U, p] = lu(A(q,q), 'vector');
            end
            f = struct('type', 1, 'L', L, 'U', U, 'p', p, 'qa', q);
        case 'LU3a'     %% 3 output LU: Gilbert-Peierls alg, permutation vectors
            q = amd(A);     %% permutation vector for AMD reordering
            [L, U, p] = lu(A(q,q), 'vector');
            f = struct('type', 1, 'L', L, 'U', U, 'p', p, 'qa', q);
        case 'LU4'      %% 4 output LU: UMFPACK, permutation vectors
            [L, U, p, q] = lu(A, 'vector');
            f = struct('type', 3, 'L', L, 'U', U, 'p', p, 'q', q);
        case 'LU5'      %% 5 output LU: UMFPACK w/row scaling, permutation vectors
            [L, U, p, q, R] = lu(A, 'vector');
            f = struct('type', 5, 'L', L, 'U', U, 'p', p, 'q', q, 'R', R);
        case 'LU3m'     %% 3 output LU: Gilbert-Peierls alg, perm mat, 1.0 piv thresh
            Q = sparse(amd(A), 1:size(A, 1), 1);    %% permutation matrix for AMD reordering
            if issparse(A)
                [L, U, P] = lu(Q'*A*Q, 1.0);
            else
                [L, U, P] = lu(Q'*A*Q);
            end
            f = struct('type', 2, 'L', L, 'U', U, 'P', P, 'Qa', Q);
        case 'LU3am'    %% 3 output LU: Gilbert-Peierls alg, permutation matrices
            Q = sparse(amd(A), 1:size(A, 1), 1);  %% permutation matrix for AMD reordering
            [L, U, P] = lu(Q'*A*Q);
            f = struct('type', 2, 'L', L, 'U', U, 'P', P, 'Qa', Q);
        case 'LU4m'     %% 4 output LU: UMFPACK, permutation matrices
            [L, U, P, Q] = lu(A);
            f = struct('type', 4, 'L', L, 'U', U, 'P', P, 'Q', Q);
        case 'LU5m'     %% 5 output LU: UMFPACK w/row scaling, permutation matrices
            [L, U, P, Q, R] = lu(A);
            f = struct('type', 6, 'L', L, 'U', U, 'P', P, 'Q', Q, 'R', R);
        case 'LU'       %% explicit LU, with options struct
            %% default options
            nout = 4;               %% 4 output args, UMFPACK
            if ~issparse(A)
                nout = 3;
            end
            vec = have_feature('lu_vec');   %% use permulation vectors, if available
            thresh = [];                    %% use default pivot threshold
            if isfield(opt, 'lu')
                opt_lu = opt.lu;
                if isfield(opt_lu, 'nout')
                    nout = opt_lu.nout;
                end
                if isfield(opt_lu, 'vec')
                    vec = opt_lu.vec;
                end
                if isfield(opt_lu, 'thresh')
                    thresh = opt_lu.thresh;
                end
            end
            %% call the appropriate form
            switch nout
                case 3      %% 3 output args: Gilbert-Peierls algorithm, with AMD reordering
                    q = amd(A);         %% permutation vector for AMD reordering
                    n = size(A, 1);
                    if vec
                        if isempty(thresh)
                            [L, U, p] = lu(A(q,q), 'vector');
                        else
                            [L, U, p] = lu(A(q,q), thresh, 'vector');
                        end
                        f = struct('type', 1, 'L', L, 'U', U, 'p', p, 'qa', q);
                    else
                        Q = sparse(q, 1:n, 1);  %% permutation matrix for AMD reordering
                        if isempty(thresh)
                            [L, U, P] = lu(Q'*A*Q);
                        else
                            [L, U, P] = lu(Q'*A*Q, thresh);
                        end
                        f = struct('type', 2, 'L', L, 'U', U, 'P', P, 'Qa', Q);
                    end
                case 4      %% 4 output args: UMFPACK
                    if vec
                        [L, U, p, q] = lu(A, 'vector');
                        f = struct('type', 3, 'L', L, 'U', U, 'p', p, 'q', q);
                    else
                        [L, U, P, Q] = lu(A);
                        f = struct('type', 4, 'L', L, 'U', U, 'P', P, 'Q', Q);
                    end
                case 5      %% 5 output args: UMFPACK w/row scaling
                    if vec
                        [L, U, p, q, R] = lu(A, 'vector');
                        f = struct('type', 5, 'L', L, 'U', U, 'p', p, 'q', q, 'R', R);
                    else
                        [L, U, P, Q, R] = lu(A);
                        f = struct('type', 6, 'L', L, 'U', U, 'P', P, 'Q', Q, 'R', R);
                    end
            end
    end
end

%% solve system
if have_b
    if is_LU        %% use LU factors
        if tr                       %% solve transposed system
            switch f.type
                case 1      %% 3 output LU: Gilbert-Peierls alg, perm vectors
                    x = zeros(size(f.L, 1), 1);
                    x(f.qa(f.p)) = f.L' \ ( f.U' \ b(f.qa) );
                case 2      %% 3 output LU: Gilbert-Peierls alg, perm matrices
                    x = f.Qa * f.P' * ( f.L' \ (f.U' \ (f.Qa' * b)) );
                case 3      %% 4 output LU: UMFPACK, perm vectors
                    x = zeros(size(f.L, 1), 1);
                    x(f.p) = f.L' \ ( f.U' \ b(f.q) );
                case 4      %% 4 output LU: UMFPACK, perm matrices
                    x = f.P' * ( f.L' \ (f.U' \ (f.Q' * b)) );
                case 5      %% 5 output LU: UMFPACK w/row scaling, perm vectors
                    x = f.R(f.p, :) \ (f.L' \ ( f.U' \ b(f.q) ));
                case 6      %% 5 output LU: UMFPACK w/row scaling, perm vectors
                    x = f.R \ (f.P' * ( f.L' \ (f.U' \ (f.Q' * b))) );
            end
        else                        %% solve non-transposed system
            switch f.type
                case 1      %% 3 output LU: Gilbert-Peierls alg, perm vectors
                    x = zeros(size(f.L, 1), 1);
                    x(f.qa) = f.U \ ( f.L \ b(f.qa(f.p)) );
                case 2      %% 3 output LU: Gilbert-Peierls alg, perm matrices
                    x = f.Qa * ( f.U \ (f.L \ (f.P * f.Qa' * b)) );
                case 3      %% 4 output LU: UMFPACK, perm vectors
                    x = zeros(size(f.L, 1), 1);
                    x(f.q) = f.U \ ( f.L \ b(f.p) );
                case 4      %% 4 output LU: UMFPACK, perm matrices
                    x = f.Q * ( f.U \ (f.L \ (f.P * b)) );
                case 5      %% 5 output LU: UMFPACK w/row scaling, perm vectors
                    x = zeros(size(f.L, 1), 1);
                    x(f.q) = f.U \ ( f.L \ (f.R(:, f.p) \ b));
                case 6      %% 5 output LU: UMFPACK w/row scaling, perm vectors
                    x = f.Q * ( f.U \ (f.L \ (f.P * (f.R \ b))) );
            end
        end
    else        %% not LU
        if tr                       %% use transposed system
            A = A';
        end

        switch solver
            case {'', '\'}  %% use built-in \ operator
                x = A \ b;
            case {'PARDISO'}
                %% set default options
                verbose = false;
                mtype = 11;
                pardiso_solver = 0;

                %% override if provided via opt
                if ~isempty(opt) && isfield(opt, 'pardiso')
                    if isfield(opt.pardiso, 'verbose') && opt.pardiso.verbose
                        verbose = true;
                    end
                    if isfield(opt.pardiso, 'mtype')
                        mtype = opt.pardiso.mtype;
                    end
                    if isfield(opt.pardiso, 'solver')
                        pardiso_solver = opt.pardiso.solver;
                    end
                end

                %% begin setup and solve
                v6 = have_feature('pardiso_object');
                if v6               %% PARDISO v6+
                    id = 1;
                    p = pardiso(id, mtype, pardiso_solver);
                    if verbose
                        p.verbose();
                    end
                else                %% PARDISO v5
                    p = pardisoinit(mtype, pardiso_solver);
                end
                if ~isempty(opt) && isfield(opt, 'pardiso')
                    if isfield(opt.pardiso, 'iparm') && ~isempty(opt.pardiso.iparm)
                        p.iparm(opt.pardiso.iparm(:, 1)) = opt.pardiso.iparm(:, 2);
                    end
                    if isfield(opt.pardiso, 'dparm') && ~isempty(opt.pardiso.dparm)
                        p.dparm(opt.pardiso.dparm(:, 1)) = opt.pardiso.dparm(:, 2);
                    end
                end
                if v6 || abs(mtype) == 2 || mtype == 6  %% need non-zero diagonal
                    nx = size(A, 1);
                    if abs(mtype) == 2 || mtype == 6    %% symmetric
                        myeps = 1e-14;
                        A = tril(A);
                    else                                %% non-symmetric
                        myeps = 1e-8;
                    end
                    A = A + myeps * speye(nx, nx);
                end
                if v6
                    p.factorize(id, A);
                    x = p.solve(id, A, b);
                    p.free(id);
                    p.clear();
                else
                    p = pardisoreorder(A, p, verbose);
                    p = pardisofactor(A, p, verbose);
                    [x, p] = pardisosolve(A, b, p, verbose);
                    pardisofree(p);
                end
            otherwise
                warning('mplinsolve: ''%s'' is not a valid value for SOLVER, using default.', solver);
                x = A \ b;
        end
    end
end

%% return only factors
if ~have_b
    x = f;
end
