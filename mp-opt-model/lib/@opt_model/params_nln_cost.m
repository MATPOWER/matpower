function varargout = params_nln_cost(om, varargin)
% params_nln_cost - Returns cost parameters for general nonlinear costs.
%
% .. note::
%    .. deprecated:: 4.3 Please use mp.sm_nln_cost.params instead, as
%       in ``om.nlc.params(...)``.
%
% ::
%
%   [N, FCN] = OM.PARAMS_NLN_COST(NAME)
%   [N, FCN] = OM.PARAMS_NLN_COST(NAME, IDX_LIST)
%   [N, FCN, VS] = OM.PARAMS_NLN_COST(...)
%
%   Returns the parameters N and FCN provided when the corresponding
%   named general nonlinear cost set was added to the model. Likewise
%   for indexed named sets specified by NAME and IDX_LIST.
%
%   An optional 3rd output argument VS indicates the variable sets used by
%   this cost set.
%
% See also opt_model, add_nln_cost, eval_nln_cost.

%   MP-Opt-Model
%   Copyright (c) 2017-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

[varargout{1:nargout}] = om.nlc.params(om.var, varargin{:});
