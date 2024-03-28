function [p, dp] = pne_pfcn_pseudo_arc_len(x, xp, step, zp)
% pne_pfcn_pseudo_arc_len - Pseudo arc length parameterization function for PNE.
% ::
%
%   [P, DP] = PNE_PFCN_PSEUDO_ARC_LEN(X, XP, STEP, ZP)
%
%   Inputs:
%       X    : solution vector x (last element is parameter lambda)
%       XP   : previous solution vector
%       STEP : continuation parameter step size
%       ZP   : normalized tangent vector at XP
%
%   This function defines a pseudo arc length parameterization for a PNE,
%   where the current point on the solution curve is constrained to lie in
%   the hyperplane running through the predicted solution orthogonal to the
%   tangent line from the previous corrected solution.
%
%   Outputs:
%       P : value of parameterization function
%       DP : Jacobian of paramerization function (transpose of gradient)
%
% See also pne_pfcn_natural, pne_pfcn_arc_len.

%   MP-Opt-Model
%   Copyright (c) 2013-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

dp = zp';                   %% derivative
p = dp * (x - xp) - step;   %% function
