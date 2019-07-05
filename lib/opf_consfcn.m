function [h, g, dh, dg] = opf_consfcn(x, om, Ybus, Yf, Yt, mpopt, il, dhs, dgs)
%OPF_CONSFCN  Evaluates nonlinear constraints and their Jacobian for OPF.
%   [H, G, DH, DG] = OPF_CONSFCN(X, OM, YBUS, YF, YT, MPOPT, IL)
%
%   Constraint evaluation function for AC optimal power flow, suitable
%   for use with MIPS or FMINCON. Computes constraint vectors and their
%   gradients.
%
%   Inputs:
%     X : optimization vector
%     OM : OPF model object
%     YBUS : bus admittance matrix
%     YF : admittance matrix for "from" end of constrained branches
%     YT : admittance matrix for "to" end of constrained branches
%     MPOPT : MATPOWER options struct
%     IL : (optional) vector of branch indices corresponding to
%          branches with flow limits (all others are assumed to be
%          unconstrained). The default is [1:nl] (all branches).
%          YF and YT contain only the rows corresponding to IL.
%     DHS : (optional) sparse matrix with tiny non-zero values specifying
%          the fixed sparsity structure that the resulting DH should match
%     DGS : (optional) sparse matrix with tiny non-zero values specifying
%          the fixed sparsity structure that the resulting DG should match
%
%   Outputs:
%     H  : vector of inequality constraint values (flow limits)
%          where the flow can be apparent power, real power, or
%          current, depending on the value of opf.flow_lim in MPOPT
%          (only for constrained lines), normally expressed as
%          (limit^2 - flow^2), except when opf.flow_lim == 'P',
%          in which case it is simply (limit - flow).
%     G  : vector of equality constraint values (power balances)
%     DH : (optional) inequality constraint gradients, column j is
%          gradient of H(j)
%     DG : (optional) equality constraint gradients
%
%   Examples:
%       [h, g] = opf_consfcn(x, om, Ybus, Yf, Yt, mpopt);
%       [h, g, dh, dg] = opf_consfcn(x, om, Ybus, Yf, Yt, mpopt);
%       [h, g, dh, dg] = opf_consfcn(x, om, Ybus, Yf, Yt, mpopt, il);
%
%   See also OPF_COSTFCN, OPF_HESSFCN.

%   MATPOWER
%   Copyright (c) 1996-2017, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargout == 2     %% contraints only
    g = om.eval_nln_constraint(x, 1);       %% equalities (power flow)
    h = om.eval_nln_constraint(x, 0);       %% inequalities (branch flow limits)
else                %% constraints and derivatives
    [g, dg] = om.eval_nln_constraint(x, 1); %% equalities (power flow)
    [h, dh] = om.eval_nln_constraint(x, 0); %% inequalities (branch flow limits)
    dg = dg';
    dh = dh';

    %% force specified sparsity structure
    if nargin > 7
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
