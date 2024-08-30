function varargout = eval_nln_constraint(om, x, iseq, varargin)
% eval_nln_constraint - Builds and returns full set of nonlinear constraints.
%
% .. note::
%    .. deprecated:: 4.3 Please use mp.sm_nln_constraint.eval instead, as
%       in ``om.nle.eval(...)`` or ``om.nli.eval(...)``.
%
% ::
%
%   G = OM.EVAL_NLN_CONSTRAINT(X, ISEQ)
%   G = OM.EVAL_NLN_CONSTRAINT(X, ISEQ, NAME)
%   G = OM.EVAL_NLN_CONSTRAINT(X, ISEQ, NAME, IDX_LIST)
%   [G, DG] = OM.EVAL_NLN_CONSTRAINT(...)
%   Builds the nonlinear equality or inequality constraints (ISEQ equal to
%   1 or 0, respectively) and optionally their gradients for the full set
%   of constraints or an individual named subset for a given value of the
%   vector X, based on constraints added by ADD_NLN_CONSTRAINT.
%
%       g(X) = 0
%       h(X) <= 0
%
%   Example:
%       [g, dg] = om.eval_nln_constraint(x, 1)
%       [h, dh] = om.eval_nln_constraint(x, 0)
%
% See also opt_model, add_nln_constraint, eval_nln_constraint_hess.

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

[varargout{1:nargout}] = om_nlx.eval(om.var, x, varargin{:});
