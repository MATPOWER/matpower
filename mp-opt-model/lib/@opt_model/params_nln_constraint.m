function varargout = params_nln_constraint(om, iseq, varargin)
% params_nln_constraint - Returns parameters for general nonlinear constraints.
%
% .. note::
%    .. deprecated:: 5.0 Please use mp.sm_nln_constraint.params instead, as
%       in ``om.nle.params(...)`` or ``om.nli.params(...)``.
%
% ::
%
%   N = OM.PARAMS_NLN_CONSTRAINT(ISEQ, NAME)
%   N = OM.PARAMS_NLN_CONSTRAINT(ISEQ, NAME, IDX_LIST)
%   [N, FCN] = OM.PARAMS_NLN_CONSTRAINT(...)
%   [N, FCN, HESS] = OM.PARAMS_NLN_CONSTRAINT(...)
%   [N, FCN, HESS, VS] = OM.PARAMS_NLN_CONSTRAINT(...)
%   [N, FCN, HESS, VS, INCLUDE] = OM.PARAMS_NLN_CONSTRAINT(...)
%
%   Returns the parameters N, and optionally FCN, and HESS provided when
%   the corresponding named nonlinear constraint set was added to the
%   model. Likewise for indexed named sets specified by NAME and IDX_LIST.
%   The ISEQ input should be set to 1 for equality constrainst and to 0
%   for inequality constraints.
%
%   An optional 4th output argument VS indicates the variable sets used by
%   this constraint set.

%   And, for constraint sets whose functions compute the constraints for
%   another set, an optional 5th output argument returns a struct with a
%   cell array of set names in the 'name' field and an array of
%   corresponding dimensions in the 'N' field.
%
% See also opt_model, add_nln_constraint, eval_nln_constraint.

%   MP-Opt-Model
%   Copyright (c) 2017-2024, Power Systems Engineering Research Center (PSERC)
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

[varargout{1:nargout}] = om_nlx.params(om.var, varargin{:});
