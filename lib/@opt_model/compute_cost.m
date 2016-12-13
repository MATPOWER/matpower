function f = compute_cost(om, x, name, idx)
%COMPUTE_COST  Computes a user-defined cost.
%   F_U = COMPUTE_COST(OM, X)
%   F_U = COMPUTE_COST(OM, X, NAME)
%   F_U = COMPUTE_COST(OM, X, NAME, IDX)
%
%   Computes the value of a user defined cost, either for all user
%   defined costs or for a named set of costs. Requires calling
%   BUILD_COST_PARAMS first to build the full set of parameters.
%
%   Let X be the full set of optimization variables and F_U(X, CP) be the
%   user-defined cost at X, corresponding to the set of cost parameters in
%   the CP struct returned by GET_COST_PARAMS, where CP is a struct with the
%   following fields:
%       N      - nw x nx sparse matrix
%       Cw     - nw x 1 vector
%       H      - nw x nw sparse matrix (optional, all zeros by default)
%       dd, mm - nw x 1 vectors (optional, all ones by default)
%       rh, kk - nw x 1 vectors (optional, all zeros by default)
%
%   These parameters are used as follows to compute F_U(X, CP)
%
%       R  = N*x - rh
%
%               /  kk(i),  R(i) < -kk(i)
%       K(i) = <   0,     -kk(i) <= R(i) <= kk(i)
%               \ -kk(i),  R(i) > kk(i)
%
%       RR = R + K
%
%       U(i) =  /  0, -kk(i) <= R(i) <= kk(i)
%               \  1, otherwise
%
%       DDL(i) = /  1, dd(i) = 1
%                \  0, otherwise
%
%       DDQ(i) = /  1, dd(i) = 2
%                \  0, otherwise
%
%       Dl = diag(mm) * diag(U) * diag(DDL)
%       Dq = diag(mm) * diag(U) * diag(DDQ)
%
%       w = (Dl + Dq * diag(RR)) * RR
%
%       F_U(X, CP) = 1/2 * w'*H*w + Cw'*w
%
%   See also OPT_MODEL, ADD_COSTS, BUILD_COST_PARAMS, GET_COST_PARAMS.

%   MATPOWER
%   Copyright (c) 2008-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 3
    cp = get_cost_params(om);
elseif nargin < 4
    cp = get_cost_params(om, name);
else
    cp = get_cost_params(om, name, idx);
end

[N, Cw, H, dd, rh, kk, mm] = deal(cp.N, cp.Cw, cp.H, cp.dd, ...
                                    cp.rh, cp.kk, cp.mm);
nw = size(N, 1);
r = N * x - rh;                 %% Nx - rhat
iLT = find(r < -kk);            %% below dead zone
iEQ = find(r == 0 & kk == 0);   %% dead zone doesn't exist
iGT = find(r > kk);             %% above dead zone
iND = [iLT; iEQ; iGT];          %% rows that are Not in the Dead region
iL = find(dd == 1);             %% rows using linear function
iQ = find(dd == 2);             %% rows using quadratic function
LL = sparse(iL, iL, 1, nw, nw);
QQ = sparse(iQ, iQ, 1, nw, nw);
kbar = sparse(iND, iND, [   ones(length(iLT), 1);
                            zeros(length(iEQ), 1);
                            -ones(length(iGT), 1)], nw, nw) * kk;
rr = r + kbar;                  %% apply non-dead zone shift
M = sparse(iND, iND, mm(iND), nw, nw);  %% dead zone or scale
diagrr = sparse(1:nw, 1:nw, rr, nw, nw);

%% linear rows multiplied by rr(i), quadratic rows by rr(i)^2
w = M * (LL + QQ * diagrr) * rr;

f = full((w' * H * w) / 2 + Cw' * w);
