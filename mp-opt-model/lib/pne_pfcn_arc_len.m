function [p, dp] = pne_pfcn_arc_len(x, xp, step, zp)
% pne_pfcn_arc_len - Arc length parameterization function for PNE.
% ::
%
%   [P, DP] = PNE_PFCN_ARC_LEN(X, XP, STEP)
%   [P, DP] = PNE_PFCN_ARC_LEN(X, XP, STEP, ZP)
%
%   Inputs:
%       X    : solution vector x (last element is parameter lambda)
%       XP   : previous solution vector
%       STEP : continuation parameter step size
%       ZP   : normalized tangent vector at XP (not used by this
%              parameterization)
%
%   This function defines an arc length parameterization for a PNE, where the
%   step size is the 2-norm distance of the current solution from the previous
%   solution.
%
%   Outputs:
%       P : value of parameterization function
%       DP : Jacobian of paramerization function (transpose of gradient)
%
% See also pne_pfcn_natural, pne_pfcn_pseudo_arc_len.

%   MP-Opt-Model
%   Copyright (c) 2013-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% function
p = sum((x - xp).^2) - step^2;

%% derivative
if nargout > 1
    dp = 2 * [x - xp]';
    if x(end) == xp(end)    %% first step
        dp(end) = 1.0;      %% avoid singular Jacobian from dp = 0
    end
end
