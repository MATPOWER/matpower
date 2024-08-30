function varargout = eval_quad_cost(om, varargin)
% eval_quad_cost - Evaluates individual or full set of quadratic costs.
%
% .. note::
%    .. deprecated:: 4.3 Please use mp.sm_quad_cost.eval instead, as
%       in ``om.qdc.eval(...)``.
%
% ::
%
%   F = OM.EVAL_QUAD_COST(X ...)
%   [F, DF] = OM.EVAL_QUAD_COST(X ...)
%   [F, DF, D2F] = OM.EVAL_QUAD_COST(X ...)
%   [F, DF, D2F] = OM.EVAL_QUAD_COST(X, NAME)
%   [F, DF, D2F] = OM.EVAL_QUAD_COST(X, NAME, IDX_LIST)
%   Evaluates the cost function and its derivatives for an individual named
%   set or the full set of quadratic costs for a given value of the
%   optimization vector X, based on costs added by ADD_QUAD_COST.
%
%   Example:
%       [f, df, d2f] = om.eval_quad_cost(x)
%       [f, df, d2f] = om.eval_quad_cost(x, name)
%       [f, df, d2f] = om.eval_quad_cost(x, name, idx_list)
%
% See also opt_model, add_quad_cost, params_quad_cost.

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

[varargout{1:nargout}] = om.qdc.eval(om.var, varargin{:});
