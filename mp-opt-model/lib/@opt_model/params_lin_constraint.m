function varargout = params_lin_constraint(om, varargin)
% params_lin_constraint - Builds and returns linear constraint parameters.
%
% .. note::
%    .. deprecated:: 5.0 Please use mp.sm_lin_constraint.params instead, as
%       in ``om.lin.params(...)``.
%
% ::
%
%   [A, L, U] = OM.PARAMS_LIN_CONSTRAINT()
%   [A, L, U] = OM.PARAMS_LIN_CONSTRAINT(NAME)
%   [A, L, U] = OM.PARAMS_LIN_CONSTRAINT(NAME, IDX_LIST)
%   [A, L, U, VS] = OM.PARAMS_LIN_CONSTRAINT(...)
%   [A, L, U, VS, I1, IN] = OM.PARAMS_LIN_CONSTRAINT(...)
%   [A, L, U, VS, I1, IN, TR] = OM.PARAMS_LIN_CONSTRAINT(NAME ...)
%
%   With no input parameters, it assembles and returns the parameters
%   for the aggregate linear constraints from all linear constraint sets
%   added using ADD_LIN_CONSTRAINT. The values of these parameters are
%   cached for subsequent calls. The parameters are A, L and U where the
%   linear constraint is of the form
%       L <= A * x <= U
%
%   If a NAME is provided then it simply returns the parameters for the
%   corresponding named set. Likewise for indexed named sets specified
%   by NAME and IDX_LIST. An optional 4th output argument VS indicates the
%   variable sets used by this constraint set. The size of A will be
%   consistent with VS. Optional 5th and 6th output arguments I1 and IN
%   indicate the starting and ending row indices of the corresponding
%   constraint set in the full aggregate constraint matrix. Finally, TR
%   will be true if it was the transpose of the A matrix that was
%   supplied/stored for this set (NAME or NAME/IDX_LIST must be supplied).
%
%   Examples:
%       [A, l, u] = om.params_lin_constraint();
%       [A, l, u, vs, i1, iN] = om.params_lin_constraint('Pmis');
%
% See also opt_model, add_lin_constraint.

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

[varargout{1:nargout}] = om.lin.params(om.var, varargin{:});
