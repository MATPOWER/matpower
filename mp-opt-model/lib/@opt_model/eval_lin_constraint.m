function varargout = eval_lin_constraint(om, varargin)
% eval_lin_constraint - Builds and returns full set of linear constraints.
%
% .. note::
%    .. deprecated:: 4.3 Please use mp.sm_lin_constraint.eval instead, as
%       in ``om.lin.eval(...)``.
%
% ::
%
%   AX_U = OM.EVAL_LIN_CONSTRAINT(X)
%   AX_U = OM.EVAL_LIN_CONSTRAINT(X, NAME)
%   AX_U = OM.EVAL_LIN_CONSTRAINT(X, NAME, IDX_LIST)
%   [AX_U, L_AX] = OM.EVAL_LIN_CONSTRAINT(...)
%   [AX_U, L_AX, A] = OM.EVAL_LIN_CONSTRAINT(...)
%   Builds the linear constraints for the full set of constraints or an
%   individual named subset for a given value of the vector X, based on
%   constraints added by ADD_LIN_CONSTRAINT.
%
%       l <= A * x <= u
%
%   Returns A*X - U, and optionally L - A*X and A.
%
%   Example:
%       [Ax_u, l_Ax, A] = om.eval_lin_constraint(x)
%
% See also opt_model, add_lin_constraint, params_lin_constraint.

%   MP-Opt-Model
%   Copyright (c) 2020-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

[varargout{1:nargout}] = om.lin.eval(om.var, varargin{:});
