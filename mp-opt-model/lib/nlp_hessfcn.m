function Lxx = nlp_hessfcn(om, x, lambda, cost_mult, Hs)
% nlp_hessfcn - Evaluates Hessian of Lagrangian.
% ::
%
%   LXX = NLP_HESSFCN(OM, X, LAMBDA, COST_MULT)
%   LXX = NLP_HESSFCN(OM, X, LAMBDA, COST_MULT, HS)
%
%   Hessian evaluation function, suitable for use with MIPS, FMINCON, etc.
%
%   Inputs:
%     OM : Opt-Model object
%     X : optimization vector
%     LAMBDA (struct)
%       .eqnonlin : Lagrange multipliers on equality constraints
%       .ineqnonlin : Kuhn-Tucker multipliers on inequality constraints
%     COST_MULT : (optional) Scale factor to be applied to the cost
%          (default = 1).
%     HS : (optional) sparse matrix with tiny non-zero values specifying
%          the fixed sparsity structure that the resulting LXX should match
%
%   Outputs:
%     LXX : Hessian of the Lagrangian.
%
%   Examples:
%       Lxx = nlp_hessfcn(om, x, lambda, cost_mult);
%       Lxx = nlp_hessfcn(om, x, lambda, cost_mult, Hs);
%
% See also nlp_costfcn, nlp_consfcn.

%   MP-Opt-Model
%   Copyright (c) 1996-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% ----- evaluate d2f -----
[f, df, d2f] = nlp_costfcn(om, x);
d2f = d2f * cost_mult;

%%----- evaluate Hessian of equality constraints -----
d2G = om.eval_nln_constraint_hess(x, lambda.eqnonlin, 1);

%%----- evaluate Hessian of inequality constraints -----
d2H = om.eval_nln_constraint_hess(x, lambda.ineqnonlin, 0);

%%-----  do numerical check using (central) finite differences  -----
if 0
    nx = length(x);
    step = 1e-5;
    num_d2f = sparse(nx, nx);
    num_d2G = sparse(nx, nx);
    num_d2H = sparse(nx, nx);
    for i = 1:nx
        xp = x;
        xm = x;
        xp(i) = x(i) + step/2;
        xm(i) = x(i) - step/2;
        % evaluate cost & gradients
        [fp, dfp] = nlp_costfcn(om, xp);
        [fm, dfm] = nlp_costfcn(om, xm);
        % evaluate constraints & gradients
        [Hp, Gp, dHp, dGp] = nlp_consfcn(om, xp);
        [Hm, Gm, dHm, dGm] = nlp_consfcn(om, xm);
        num_d2f(:, i) = cost_mult * (dfp - dfm) / step;
        num_d2G(:, i) = (dGp - dGm) * lambda.eqnonlin   / step;
        num_d2H(:, i) = (dHp - dHm) * lambda.ineqnonlin / step;
    end
    d2f_err = full(max(max(abs(d2f - num_d2f))));
    d2G_err = full(max(max(abs(d2G - num_d2G))));
    d2H_err = full(max(max(abs(d2H - num_d2H))));
    if d2f_err > 1e-6
        fprintf('Max difference in d2f: %g\n', d2f_err);
    end
    if d2G_err > 1e-5
        fprintf('Max difference in d2G: %g\n', d2G_err);
    end
    if d2H_err > 1e-6
        fprintf('Max difference in d2H: %g\n', d2H_err);
    end
end

Lxx = d2f + d2G + d2H;

%% force specified sparsity structure
if nargin > 4
    %% add sparse structure (with tiny values) to current matrices to
    %% ensure that sparsity structure matches that supplied
    Lxx = Lxx + Hs;

%     %% check sparsity structure against that supplied
%     if nnz(Lxx) ~= nnz(Hs)
%         fprintf('=====> nnz(Lxx) is %d, expected %d <=====\n', nnz(Lxx), nnz(Hs));
%     else
%         [iHs, jHs] = find(Hs);
%         [iH, jH] = find(Lxx);
%         if any(iH ~= iHs) || any(jH ~= jHs)
%             fprintf('=====> structure of Lxx is not as expected <=====\n');
%         end
%     end
end
