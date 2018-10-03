function [f, df, d2f] = opf_gen_cost_fcn(x, baseMVA, gencost, ig, mpopt)
%OPF_GEN_COST_FCN  Evaluates polynomial generator costs and derivatives.
%   [F, DF, D2F] = OPF_GEN_COST_FCN(X, BASEMVA, COST, MPOPT)
%
%   Evaluates the polynomial generator costs and derivatives.
%
%   Inputs:
%     X : single-element cell array with vector of (active or reactive)
%         dispatches (in per unit)
%     BASEMVA : system per unit base
%     GENCOST : standard gencost matrix corresponding to dispatch
%               (active or reactive) provided in X
%     IG : vector of generator indices of interest
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     F  : sum of generator polynomial costs
%     DF : (optional) gradient (column vector) of polynomial costs
%     D2F : (optional) Hessian of polynomial costs
%
%   Examples:
%       f = opf_gen_cost_fcn(x, baseMVA, gencost, ig, mpopt);
%       [f, df] = opf_gen_cost_fcn(x, baseMVA, gencost, ig, mpopt);
%       [f, df, d2f] = opf_gen_cost_fcn(x, baseMVA, gencost, ig, mpopt);

%   MATPOWER
%   Copyright (c) 1996-2017, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialize -----
if isempty(ig)
    [PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
    ig = find(gencost(:, MODEL) == POLYNOMIAL);     %% poly MW / MVAr costs
end
PQg = x{1};                 %% active or reactive dispatch in p.u.
ng = length(PQg);           %% number of dispatchable injections

%%----- evaluate objective function -----
xx = PQg(ig) * baseMVA;     %% active or reactive dispatch in MW/MVAr
f = sum( totcost(gencost(ig, :), xx) ); %% cost of poly P or Q

%%----- evaluate cost gradient -----
if nargout > 1
    %% polynomial cost of P and Q
    df = zeros(ng, 1);      %% w.r.t p.u. Pg / Qg
    df(ig) = baseMVA * polycost(gencost(ig, :), xx, 1);

    %% ---- evaluate cost Hessian -----
    if nargout > 2
        %% polynomial generator costs
        d2f = sparse(ig, ig, baseMVA^2 * polycost(gencost(ig, :), xx, 2), ng, ng);
    end
end
