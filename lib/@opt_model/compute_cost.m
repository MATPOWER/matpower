function f = compute_cost(om, x, varargin)
%COMPUTE_COST  Computes a user-defined cost.
%
%   -----  DEPRECATED - use EVAL_LEGACY_COST instead  -----
%
%   F_U = OM.COMPUTE_COST(X)
%   F_U = OM.COMPUTE_COST(X, NAME)
%   F_U = OM.COMPUTE_COST(X, NAME, IDX)
%
%   Computes the value of a user defined cost, either for all user
%   defined costs or for a named set of costs. When specifying a named set,
%   if NAME refers to an indexed name but IDX is not provided, then
%   COMPUTE_COST will be called recursively for all corresponding values
%   of IDX and the costs summed.
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
%   See also OPT_MODEL, ADD_COSTS, PARAMS_LEGACY_COST.

%   MATPOWER
%   Copyright (c) 2008-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

f = om.eval_legacy_cost(x, varargin{:});
