function varargout = params_var(om, varargin)
% params_var - Returns initial value, lower bound and upper bound for opt variables.
%
% .. note::
%    .. deprecated:: 4.3 Please use mp.sm_variable.params instead, as
%       in ``om.var.params(...)``.
%
% ::
%
%   [V0, VL, VU] = OM.PARAMS_VAR()
%   [V0, VL, VU] = OM.PARAMS_VAR(NAME)
%   [V0, VL, VU] = OM.PARAMS_VAR(NAME, IDX_LIST)
%   [V0, VL, VU, VT] = PARAMS_VAR(...)
%
%   Returns the initial value V0, lower bound VL and upper bound VU for
%   the full set of optimization variables, or for a specific named or named
%   and indexed variable set. Values for the full set are cached for
%   subsequent calls. Optionally also returns a corresponding char
%   vector VT of variable types, where 'C', 'I' and 'B' represent continuous
%   integer and binary variables, respectively.
%
%   Examples:
%       [x0, xmin, xmax] = om.params_var();
%       [Pg0, Pmin, Pmax] = om.params_var('Pg');
%       [zij0, zijmin, zijmax, ztype] = om.params_var('z', {i, j});
%   
% See also opt_model, add_var.

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

[varargout{1:nargout}] = om.var.params(varargin{:});
