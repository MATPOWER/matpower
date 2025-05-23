function [ieq, igt, ilt, Qe, Be, de, Qi, Bi, di] = convert_quad_constraint(Q, B, l, u)
% convert_quad_constraint - Convert quadratic constraints from bounded to equality/inequality pair.
% ::
%
%   [ieq, igt, ilt, Qe, Be, de, Qi, Bi, di] = convert_quad_constraint(Q, B, l, u)
%   [ieq, igt, ilt, Q, B, d] = convert_quad_constraint(Q, B, l, u)
%
% Convert a set of constraints of the form::
%
%       l(j) <= x'*Q{j}*x + B(j,:)*x <= u(j),  j = 1,2,...,nq
%
% to::
%
%       x'*Qe{j}*x + Be(j,:)*x  =  de(j),  j = 1,2,...,nqe
%       x'*Qi{j}*x + Bi(j,:)*x  <= di(j),  j = 1,2,...,nqi   , nq = nqe+nqi
%
% where::
%
%       Qe = Q(ieq)
%       Be = B(ieq,:)
%       de = u(ieq)
%       Qi = [Q(ilt); -Q(igt)]
%       Bi = [B(ilt,:); -B(igt,:)]
%       di = [u(ilt);  -l(igt)]
%
% Alternatively, the returned cell arrays, matrices, and RHS vectors can be
% stacked into a single set with the equalities first, then the inequalities.
% ::
%
%       Q = [Qe; Qi]
%       B = [Be; Bi]
%       d = [de; di]
%
% Inputs:
%   Q (double) : :math:`n_q \times 1`, cell array with sparse :math:`n \times n`
%       symmetric matrices holding the quadratic parameters of the quadratic
%       constraints
%   B (double) : :math:`n_q \times n`, matrix whose rows represent the linear
%                parameters of the quadratic constraints
%   l (double) : :math:`n_q \times 1`, quadratic constraint lower bound vector
%   u (double) : :math:`n_q \times 1`, quadratic constraint upper bound vector
%
% Outputs:
%   ieq (integer) : vector of indices of equality constraints
%   igt (integer) : vector of indices of lower bounded inequality constraints
%   ilt (integer) : vector of indices of upper bounded inequality constraints
%   Qe (double)   : cell array of quadratic parameters for equality constraints
%   Ae (double)   : matrix of linear parameters for equality constraints
%   de (double)   : equality constraint RHS
%   Qi (double)   : cell array of quadratic parameters for inequality constraints
%   Ai (double)   : matrix of linear terms for inequality constraints
%   di (double)   : inequality constraint RHS
%   Q (double)    : stacked cell of quadratic parameters
%   B (double)    : stacked constraint matrix of linear parameters
%   d (double)    : stacked RHS
%
% See also convert_constraint_multipliers.

%   MP-Opt-Model
%   Copyright (c) 2019-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model..
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

eq = abs(u-l) <= eps;
ieq = find( eq );               %% equality
igt = find( l > -1e10 & ~eq );  %% inequality w/finite lower bound
ilt = find( u <  1e10 & ~eq );  %% inequality w/finite upper bound

if nargout < 9      %% stacked case
    Qe = [Q([ieq; ilt]); cellfun(@(x)(-1*x), Q(igt), 'UniformOutput', false)];
    Be = [B([ieq; ilt],:); -B(igt,:)];
    de = [u([ieq; ilt]);  -l(igt)];
else                %% equality/inequality pairs
    Qe = Q(ieq);
    Be = B(ieq,:);
    de = u(ieq);
    Qi = [Q(ilt); cellfun(@(x)(-1*x), Q(igt), 'UniformOutput', false)];
    Bi = [B(ilt,:); -B(igt,:)];
    di = [u(ilt);  -l(igt)];
end