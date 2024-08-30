function d2G = eval_nln_constraint_hess(om, x, lam, iseq)
% eval_nln_constraint_hess - Builds and returns Hessian of nonlinear constraints.
%
% .. note::
%    .. deprecated:: 4.3 Please use mp.sm_nln_constraint.eval_hess instead, as
%       in ``om.nle.eval_hess(...)`` or ``om.nli.eval_hess(...)``.
%
% ::
%
%   D2G = OM.EVAL_NLN_CONSTRAINT_HESS(X, LAM, ISEQ)
%   Builds the Hessian of the full set of nonlinear equality or inequality
%   constraints (ISEQ equal to 1 or 0, respectively) for given values of
%   the optimization vector X and dual variables LAM, based on constraints
%   added by ADD_NLN_CONSTRAINT.
%
%       g(X) = 0
%       h(X) <= 0
%
%   Example:
%       d2G = om.eval_nln_constraint_hess(x, lam, 1)
%       d2H = om.eval_nln_constraint_hess(x, lam, 0)
%
% See also opt_model, add_nln_constraint, eval_nln_constraint.

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% get constraint type
if iseq         %% equality constraints
    om_nlx = om.nle;
else            %% inequality constraints
    om_nlx = om.nli;
end

d2G = om_nlx.eval_hess(om.var, x, lam);
