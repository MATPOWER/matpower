function [h, g, dh, dg] = qcqp_nlp_consfcn(x, QQ, BB, dd)
% qcqp_nlp_consfcn - Evaluates quadratic constraints and their Jacobian for NLP solver.
% ::
%
%   [H, G] = QCQP_NLP_CONSFCN(X, QQ, BB, DD)
%   [H, G, DH, DG] = QCQP_NLP_CONSFCN(X, QQ, BB, DD)
%
%   Constraint evaluation function for quadratic constraints, suitable
%   for use with MIPS, FMINCON, etc. Computes constraint vectors and their
%   gradients for a set of quadratic constraints of the form:
%
%       1/2 Xblk' * blkQe * Xblk + Be * X == de
%       1/2 Xblk' * blkQi * Xblk + Bi * X <= di
%
%   where Xblk is formed by creating a block diagonal matrix with X repeated
%   along the block diagonal.
%
%   Inputs:
%     X   : optimization vector
%     QQ  : struct with (possibly sparse) quadratic matrices for
%           equality/inequality constraints with the following fields:
%         blkQi : block diagonal matrix formed from the nqi x 1 cell array
%                 of sparse quadratic matrices for inequaliy constraints
%         blkQe : block diagonal matrix formed from the nqe x 1 cell array
%                 of sparse quadratic matrices for equaliy constraints
%     BB : struct with the matrices (possibly sparse) of linear parameters
%          of equality/inequality constraints with the following fields:
%         Bi : matrix with linear parameters for inequality constraints
%         Be : matrix with linear parameters for equality constraints
%     DD : struct with the vector of constant terms of equality/inequality
%          constraints with the following fields:
%         di : vector with constant terms for inequality constraints
%         de : vector with constant terms for equality constraints
%
%   Outputs:
%     H  : vector of inequality constraint values
%     G  : vector of equality constraint values
%     DH : (optional) inequality constraint gradients, column j is
%          gradient of H(j)
%     DG : (optional) equality constraint gradients
%
%   Examples:
%       [h, g] = qcqp_nlp_consfcn(x, QQ, BB, dd);
%       [h, g, dh, dg] = qcqp_nlp_consfcn(x, QQ, BB, dd);
%
% See also qcqp_nlp_costfcn, qcqp_nlp_hessfcn, qcqps_master.

%   MP-Opt-Model
%   Copyright (c) 2019-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% gather parameters for equality/inequality quadratic constraints
blkQe = QQ.blkQe; blkQi = QQ.blkQi;
Be = BB.Be; Bi = BB.Bi;
de = dd.de; di = dd.di;

%% assess constraints and (possibly) their derivatives
if nargout > 1       %% constraints
    if isempty(blkQi) || isempty(x)         %% inequalities
        h = [];
    else
        nin = size(Bi,1);
        xxineq = mat2cell(repmat(sparse(x), nin, 1), length(x)*ones(nin,1));
        blkxineq = blkdiag(xxineq{:});
        h = 1/2 * diag(blkxineq' * blkQi * blkxineq) + Bi * x - di;
    end

    if isempty(blkQe) || isempty(x)         %% equalities
        g = [];
    else
        neq = size(Be,1);
        xxeq = mat2cell(repmat(sparse(x), neq, 1), length(x)*ones(neq,1));
        blkxeq = blkdiag(xxeq{:});
        g = 1/2 * diag(blkxeq' * blkQe * blkxeq) + Be * x - de;
    end
end

if nargout > 2       %% derivatives
    if isempty(blkQi) || isempty(x)         %% inequalities
        dh = [];
    else
        blkprod = blkxineq' * blkQi;
        [col1, row1] = size(blkxineq);
        [ ~  , col2] = size(blkQi);
        N = col1 / length(x);
        id_start = [(0:N-1)*(row1/N)+1 ; (0:N-1)*(col2/N)+1 ];
        id_end   = [   (1:N)*(row1/N)  ;    (1:N)*(col2/N)  ];
        id_start = mat2cell(id_start(:), 2*ones(N,1));
        id_end   = mat2cell(id_end(:), 2*ones(N,1));
        Qix = cellfun(@(x,y)(blkprod(x(1):y(1), x(2):y(2))), id_start, id_end, 'UniformOutput', false);
        Qix = cell2mat(Qix);
        dh = (Qix + Bi)';
    end

    if isempty(blkQe) || isempty(x)                       %% equalities
        dg = [];
    else
        blkprod = blkxeq' * blkQe;
        [col1, row1] = size(blkxeq);
        [ ~  , col2] = size(blkQe);
        N = col1 / length(x);
        id_start = [(0:N-1)*(row1/N)+1 ; (0:N-1)*(col2/N)+1 ];
        id_end   = [   (1:N)*(row1/N)  ;    (1:N)*(col2/N)  ];
        id_start = mat2cell(id_start(:), 2*ones(N,1));
        id_end   = mat2cell(id_end(:), 2*ones(N,1));
        Qex = cellfun(@(x,y)(blkprod(x(1):y(1), x(2):y(2))), id_start, id_end, 'UniformOutput', false);
        Qex = cell2mat(Qex);
        dg = (Qex + Be)';
    end
end
