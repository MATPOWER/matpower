function [h, g, dh, dg] = nlp_consfcn(mm, x, dhs, dgs)
% nlp_consfcn - Evaluates nonlinear constraints and their Jacobian for NLP solver.
% ::
%
%   [H, G] = NLP_CONSFCN(MM, X)
%   [H, G, DH, DG] = NLP_CONSFCN(MM, X)
%   [H, G, DH, DG] = NLP_CONSFCN(MM, X, DHS, DGS)
%
%   Constraint evaluation function nonlinear constraints, suitable
%   for use with MIPS, FMINCON, etc. Computes constraint vectors and their
%   gradients.
%
%   Inputs:
%     MM : MP-Opt-Model object
%     X : optimization vector
%     DHS : (optional) sparse matrix with tiny non-zero values specifying
%          the fixed sparsity structure that the resulting DH should match
%     DGS : (optional) sparse matrix with tiny non-zero values specifying
%          the fixed sparsity structure that the resulting DG should match
%
%   Outputs:
%     H  : vector of inequality constraint values
%     G  : vector of equality constraint values
%     DH : (optional) inequality constraint gradients, column j is
%          gradient of H(j)
%     DG : (optional) equality constraint gradients
%
%   Examples:
%       [h, g] = nlp_consfcn(mm, x);
%       [h, g, dh, dg] = nlp_consfcn(mm, x);
%       [...] = nlp_consfcn(mm, x, dhs, dgs);
%
% See also nlp_costfcn, nlp_hessfcn.

%   MP-Opt-Model
%   Copyright (c) 1996-2025, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if nargout == 2     %% contraints only
    g = mm.nle.eval(mm.var, x);         %% equalities
    h = mm.nli.eval(mm.var, x);         %% inequalities
else                %% constraints and derivatives
    [g, dg] = mm.nle.eval(mm.var, x);   %% equalities
    [h, dh] = mm.nli.eval(mm.var, x);   %% inequalities
    dg = dg';
    dh = dh';

    %% force specified sparsity structure
    if nargin > 2
        %% add sparse structure (with tiny values) to current matrices to
        %% ensure that sparsity structure matches that supplied
        dg = dg + dgs;
        dh = dh + dhs;

%         %% check sparsity structure against that supplied
%         if nnz(dg) ~= nnz(dgs)
%             fprintf('=====> nnz(dg) is %d, expected %d <=====\n', nnz(dg), nnz(dgs));
%         else
%             [idgs, jdgs] = find(dgs);
%             [idg, jdg] = find(dg);
%             if any(idg ~= idgs) || any(jdg ~= jdgs)
%                 fprintf('=====> structure of dg is not as expected <=====\n');
%             end
%         end
%         if nnz(dh) ~= nnz(dhs)
%             fprintf('=====> nnz(dh) is %d, expected %d <=====\n', nnz(dh), nnz(dhs));
%         else
%             [idhs, jdhs] = find(dhs);
%             [idh, jdh] = find(dh);
%             if any(idh ~= idhs) || any(jdh ~= jdhs)
%                 fprintf('=====> structure of dh is not as expected <=====\n');
%             end
%         end
    end
end
