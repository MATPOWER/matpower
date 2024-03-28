function [p, dp] = pne_pfcn_natural(x, xp, step, zp)
% pne_pfcn_natural - Natural parameterization function for PNE.
% ::
%
%   [P, DP] = PNE_PFCN_NATURAL(X, XP, STEP)
%   [P, DP] = PNE_PFCN_NATURAL(X, XP, STEP, ZP)
%
%   This function defines a "natural" parameterization for a PNE, where the
%   step size is the difference between current and previous values of the
%   parameter lambda.
%
%   Inputs:
%       X    : solution vector x (last element is parameter lambda)
%       XP   : previous solution vector
%       STEP : continuation parameter step size
%       ZP   : normalized tangent vector at XP (not used by this
%              parameterization)
%
%   Outputs:
%       P : value of parameterization function
%       DP : Jacobian of paramerization function (transpose of gradient)
%
% See also pne_pfcn_arc_len, pne_pfcn_pseudo_arc_len.

%   MP-Opt-Model
%   Copyright (c) 2013-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% function
N = length(x);
p = abs(x(end) - xp(end)) - step;

%% derivative
if nargout > 1
    dp = zeros(1, N);
    if x(end) >= xp(end)
        dp(end) = 1.0;
    else
        dp(end) = -1.0;
    end
end
