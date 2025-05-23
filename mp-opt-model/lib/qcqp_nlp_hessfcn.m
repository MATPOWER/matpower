function Lxx = qcqp_nlp_hessfcn(x, lambda, H, QQi, QQe, cost_mult)
% qcqp_nlp_hessfcn - Evaluates QCQP Hessian of Lagrangian for NLP solver.
% ::
%
%   LXX = QCQP_NLP_HESSFCN(X, LAMBDA, H, QQI, QQE)
%   LXX = QCQP_NLP_HESSFCN(X, LAMBDA, H, QQI, QQE, COST_MULT)
%
%   Hessian evaluation function, suitable for use with MIPS, FMINCON, etc.
%
%   Inputs:
%     X : optimization vector
%     LAMBDA (struct)
%       .eqnonlin : Lagrange multipliers on equality constraints
%       .ineqnonlin : Kuhn-Tucker multipliers on inequality constraints
%     H   : matrix (possibly sparse) of quadratic cost coefficients
%     QQi : sparse matrix formed by stacking vertically all the quadratic
%           matrices for inequality constraints.
%     QQe : sparse matrix formed by stacking vertically all the quadratic
%           matrices for equality constraints.
%     COST_MULT : (optional) Scale factor to be applied to the cost
%          (default = 1).
%
%   Outputs:
%     LXX : Hessian of the Lagrangian.
%
%   Examples:
%       Lxx = qcqp_nlp_hessfcn(x, lambda, H, QQi, QQe);
%       Lxx = qcqp_nlp_hessfcn(x, lambda, H, QQi, QQe, cost_mult);
%
% See also qcqp_nlp_costfcn, qcqp_nlp_consfcn, qcqps_master.

%   MP-Opt-Model
%   Copyright (c) 1996-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% ----- evaluate d2f -----
if isempty(H)
    d2f = 0;
else
    d2f = H;
end

if nargin < 6
    cost_mult = 1;
end
d2f = d2f * cost_mult;

%% ----- evaluate d2H and d2G -----
n = length(x);
%%----- evaluate Hessian of inequality constraints -----
if isempty(QQi) || isempty(x)
    d2H = 0;
else
    lambdaineq = num2cell(lambda.ineqnonlin');
    Lambdaineq = cellfun(@(x)(x*speye(n)), lambdaineq, 'UniformOutput', false);
    d2H = cell2mat(Lambdaineq) * QQi;
end

%%----- evaluate Hessian of equality constraints -----
if isempty(QQe) || isempty(x)
    d2G = 0;
else
    lambdaeq = num2cell(lambda.eqnonlin');
    Lambdaeq = cellfun(@(x)(x*speye(n)), lambdaeq, 'UniformOutput', false);
    d2G = cell2mat(Lambdaeq) * QQe;
end

%% ----- evaluate Hessian of the Lagrangian -----
Lxx = d2f + d2G + d2H;
