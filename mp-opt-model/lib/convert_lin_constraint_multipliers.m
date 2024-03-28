function [mu_l, mu_u] = convert_lin_constraint_multipliers(lam, mu, ieq, igt, ilt)
% convert_lin_constraint_multipliers - Convert multipliers for eq/ineq pair to bounded.
% ::
%
%   [mu_l, mu_u] = convert_lin_constraint_multipliers(lam, mu, ieq, igt, ilt)
%
% Convert multipliers back for constraints converted by convert_lin_constraint.
% That is, we assume a set of constraints that have been converted from the
% form::
%
%   l <= A * x <= u       (1)
%
% to::
%
%   Ae * x = be           (2)
%   Ai * x <= bi          (3)
%
% where::
%
%   Ae = A(ieq, :);
%   be = u(ieq, 1);
%   Ai  = [ A(ilt, :); -A(igt, :) ];
%   bi  = [ u(ilt, 1); -l(igt, 1) ];
%
% Given the constraint multipliers for constraints of the form (2)-(3), this
% function converts them back to the equivalent multipliers for form (1).
%
% Inputs:
%   lam (double) : vector of multipliers for (2)
%   mu (double) : vector of multipliers for (3)
%   ieq (integer) : vector of indices of equality constraints
%   igt (integer) : vector of indices of lower bounded inequality constraints
%   ilt (integer) : vector of indices of upper bounded inequality constraints
%
% Outputs:
%   mu_l (double) : vector of multipliers on lower bound of (1)
%   mu_u (double) : vector of multipliers on upper bound of (1)
%
% See also convert_lin_constraint.

%   MP-Opt-Model
%   Copyright (c) 2010-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

m = max([ieq; igt; ilt]);   %% number of original constraints in (1)
nlt = length(ilt);          %% number of upper bounded inequalities

%% initialize multipliers
mu_l = zeros(m, 1);
mu_u = zeros(m, 1);

%% convert equality constraint multipliers
kl = find(lam < 0);
ku = find(lam >= 0);
mu_l(ieq(kl)) = -lam(kl);
mu_u(ieq(ku)) = lam(ku);

%% convert inequality constraint multipliers
mu_l(igt) = mu(nlt+1:end);
mu_u(ilt) = mu(1:nlt);
