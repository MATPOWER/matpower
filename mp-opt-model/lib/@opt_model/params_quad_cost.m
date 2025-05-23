function varargout = params_quad_cost(om, varargin)
% params_quad_cost - Returns the cost parameters for quadratic costs.
%
% .. note::
%    .. deprecated:: 5.0 Please use mp.sm_quad_cost.params instead, as
%       in ``om.qdc.params(...)``.
%
% ::
%
%   [Q, C] = OM.PARAMS_QUAD_COST()
%   [Q, C] = OM.PARAMS_QUAD_COST(NAME)
%   [Q, C] = OM.PARAMS_QUAD_COST(NAME, IDX_LIST)
%   [Q, C, K] = OM.PARAMS_QUAD_COST(...)
%   [Q, C, K, VS] = OM.PARAMS_QUAD_COST(...)
%
%   With no input parameters, it assembles and returns the parameters
%   for the aggregate quadratic cost from all quadratic cost sets added
%   using ADD_QUAD_COST. The values of these parameters are cached
%   for subsequent calls. The parameters are Q, C, and optionally K,
%   where the quadratic cost is of the form
%       F(X) = 1/2 * X'*Q*X + C'*X + K
%
%   If a NAME is provided then it simply returns the parameters for the
%   corresponding named set. Likewise for indexed named sets specified
%   by NAME and IDX_LIST. In this case, Q and K may be vectors, corresponding
%   to a cost function of the form
%       F(X) = 1/2 * Q .* X.^2 + C .* X + K
%
%   An optional 4th output argument VS indicates the variable sets used by
%   this cost set. The size of Q and C will be consistent with VS.
%
% See also opt_model, add_quad_cost.

%   MP-Opt-Model
%   Copyright (c) 2017-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

[varargout{1:nargout}] = om.qdc.params(om.var, varargin{:});
