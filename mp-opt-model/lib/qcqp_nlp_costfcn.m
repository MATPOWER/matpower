function [f, df, d2f] = qcqp_nlp_costfcn(x, H, c)
% qcqp_nlp_costfcn - Evaluates QCQP objective function, gradient and Hessian for NLP solver.
% ::
%
%   [F, DF, D2F] = QCQP_NLP_COSTFCN(X, H, C)
%
%   Objective function evaluation routine, suitable for use with MIPS,
%   FMINCON, etc. Computes objective function value, gradient and Hessian
%   of the quadratic function:
%
%       F(X) = 1/2 X'*H*X + C'*X
%
%   Inputs:
%     X : optimization vector
%     H : matrix (possibly sparse) of quadratic cost coefficients
%     C : vector of linear cost coefficients
%
%   Outputs:
%     F   : value of objective function
%     DF  : (optional) gradient of objective function (column vector)
%     D2F : (optional) Hessian of objective function (sparse matrix)
%
%   Examples:
%       f = qcqp_nlp_costfcn(x, H, c);
%       [f, df] = qcqp_nlp_costfcn(x, H, c);
%       [f, df, d2f] = qcqp_nlp_costfcn(x, H, c);
%
% See also qcqp_nlp_consfcn, qcqp_nlp_hessfcn, qcqps_master.

%   MP-Opt-Model
%   Copyright (c) 2019-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%%  ----- check inputs -----
if isempty(H)
    nx = length(c);
    H = sparse(nx,nx);
end

%%  ----- evaluate objective function -----
if nargout == 3
    f = 1/2 * x' * H * x + c' * x;
    df = H * x + c;
    d2f = H;
elseif nargout == 2
    f = 1/2 * x' * H * x + c' * x;
    df = H * x + c;
else
    f = 1/2 * x' * H * x + c' * x;
end
