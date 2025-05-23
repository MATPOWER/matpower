function [mu_l, mu_u] = convert_constraint_multipliers(lam, mu, ieq, igt, ilt)
% convert_constraint_multipliers - Convert multipliers for eq/ineq pair to bounded.
% ::
%
%   [mu_l, mu_u] = convert_constraint_multipliers(lam, mu, ieq, igt, ilt)
%
% Convert multipliers back for constraints converted by convert_lin_constraint
% or convert_quad_constraint. That is, we assume a set of constraints that have
% been converted from the form:
%
% .. math::
%    :label: eq_mu_conv_bounded
%
%    \l \le \g(\x) \le \u
%
% to 
%
% .. math::
%    :label: eq_mu_conv_eq
%
%    \g_e(\x) = \rvec{b}_e
%
% .. math::
%    :label: eq_mu_conv_ineq
%
%    \g_i(\x) = \rvec{b}_i
%
% The relationships, using code notation, can be summarized as::
%
%   g_e = g(ieq, :);
%   b_e = u(ieq, 1);
%   g_i  = [ g(ilt, :); -g(igt, :) ];
%   b_i  = [ u(ilt, 1); -l(igt, 1) ];
%
% Given the constraint multipliers for constraints of the form
% :eq:`eq_mu_conv_eq`-:eq:`eq_mu_conv_ineq`, this function converts them back to
% the equivalent multipliers for form :eq:`eq_mu_conv_bounded`.
%
% Inputs:
%   lam (double) : vector of multipliers for :eq:`eq_mu_conv_eq`
%   mu (double) : vector of multipliers for :eq:`eq_mu_conv_ineq`
%   ieq (integer) : vector of indices of equality constraints
%   igt (integer) : vector of indices of lower bounded inequality constraints
%   ilt (integer) : vector of indices of upper bounded inequality constraints
%
% Outputs:
%   mu_l (double) : vector of multipliers on lower bound of :eq:`eq_mu_conv_bounded`
%   mu_u (double) : vector of multipliers on upper bound of :eq:`eq_mu_conv_bounded`
%
% See also convert_lin_constraint, convert_quad_constraint.

%   MP-Opt-Model
%   Copyright (c) 2010-2025, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

m = max([ieq; igt; ilt]);   %% number of original constraints in (1)
nlt = length(ilt);          %% number of upper bounded inequalities

%% initialize multipliers
if all(isnan(lam)) && all(isnan(mu))
    mu_l = NaN(m, 1);
    mu_u = NaN(m, 1);
else
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
end
