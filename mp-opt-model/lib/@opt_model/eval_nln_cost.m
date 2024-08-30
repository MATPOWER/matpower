function varargout = eval_nln_cost(om, varargin)
% eval_nln_cost - Evaluates individual or full set of general nonlinear costs.
%
% .. note::
%    .. deprecated:: 4.3 Please use mp.sm_nln_cost.eval instead, as
%       in ``om.nlc.eval(...)``.
%
% ::
%
%   [F, DF, D2F] = OM.EVAL_NLN_COST(X)
%   [F, DF, D2F] = OM.EVAL_NLN_COST(X, NAME)
%   [F, DF, D2F] = OM.EVAL_NLN_COST(X, NAME, IDX_LIST)
%   Evaluates the cost function and its derivatives for an individual named
%   set or the full set of general nonlinear costs for a given value of the
%   optimization vector X, based on costs added by ADD_NLN_COST.
%
%   Example:
%       [f, df, d2f] = om.eval_nln_cost(x)
%       [f, df, d2f] = om.eval_nln_cost(x, name)
%       [f, df, d2f] = om.eval_nln_cost(x, name, idx_list)
%
% See also opt_model, add_nln_cost, params_nln_cost.

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

[varargout{1:nargout}] = om.nlc.eval(om.var, varargin{:});
