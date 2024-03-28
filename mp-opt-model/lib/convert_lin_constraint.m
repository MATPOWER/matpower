function [ieq, igt, ilt, Ae, be, Ai, bi] = convert_lin_constraint(A, l, u)
% convert_lin_constraint - Convert from bounded to equality/inequality pair.
% ::
%
%   [ieq, igt, ilt, Ae, be, Ai, bi] = convert_lin_constraint(A, l, u)
%   [ieq, igt, ilt, A, b] = convert_lin_constraint(A, l, u)
%
% Convert a set of constraints of the form::
%
%   l <= A * x <= u
%
% to::
%
%   Ae * x = be
%   Ai * x <= bi
%
% where::
%
%   Ae = A(ieq, :);
%   be = u(ieq, 1);
%   Ai  = [ A(ilt, :); -A(igt, :) ];
%   bi  = [ u(ilt, 1); -l(igt, 1) ];
%
% Alternatively, the returned matrices and RHS vectors can be stacked into
% a single set with the equalities first, then the inequalities.
% ::
%
%   A = [Ae; Ai]
%   b = [be; bi]
%
% Inputs:
%   A (double) : m x n, (possibly sparse) constraint matrix
%   l (double): m x 1, constraint lower bound vector
%   u (double): m x 1, constraint upper bound vector
%
% Outputs:
%   ieq (integer) : vector of indices of equality constraints
%   igt (integer) : vector of indices of lower bounded inequality constraints
%   ilt (integer) : vector of indices of upper bounded inequality constraints
%   Ae (double) : (possibly sparse) equality constraint matrix
%   be (double) : equality constraint RHS
%   Ai (double) : (possibly sparse) inequality constraint matrix
%   bi (double) : inequality constraint RHS
%   A (double) : (possibly sparse) stacked constraint matrix
%   b (double) : stacked RHS
%
% See also convert_lin_constraint_multipliers.

%   MP-Opt-Model
%   Copyright (c) 2010-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

eq = abs(u-l) <= eps;
ieq = find( eq );               %% equality
igt = find( l > -1e10 & ~eq );  %% inequality w/finite lower bound
ilt = find( u <  1e10 & ~eq );  %% inequality w/finite upper bound

if nargout < 7
    Ae  = [ A([ieq; ilt], :); -A(igt, :) ];     %% A
    be  = [ u([ieq; ilt], 1); -l(igt, 1) ];     %% b
else
    Ae = A(ieq, :);
    be = u(ieq, 1);
    Ai  = [ A(ilt, :); -A(igt, :) ];
    bi  = [ u(ilt, 1); -l(igt, 1) ];
end
