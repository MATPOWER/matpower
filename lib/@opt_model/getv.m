function varargout = getv(om, varargin)
%GETV  Returns initial value, lower bound and upper bound for opt variables.
%
%   -----  DEPRECATED - Please use PARAMS_VAR instead  -----
%
%   [V0, VL, VU] = OM.GETV()
%   [V0, VL, VU] = OM.GETV(NAME)
%   [V0, VL, VU] = OM.GETV(NAME, IDX)
%   [V0, VL, VU, VT] = GETV(...)
%   Returns the initial value V0, lower bound VL and upper bound VU for
%   the full optimization variable vector, or for a specific named or named
%   and indexed variable set. Optionally also returns a corresponding char
%   vector VT of variable types, where 'C', 'I' and 'B' represent continuous
%   integer and binary variables, respectively.
%
%   Examples:
%       [x, xmin, xmax] = om.getv();
%       [Pg, Pmin, Pmax] = om.getv('Pg');
%       [zij0, zijmin, zijmax, ztype] = om.getv('z', {i, j});
%   
%   See also OPT_MODEL, PARAMS_VAR.

%   MATPOWER
%   Copyright (c) 2008-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

[varargout{1:nargout}] = om.params_var(varargin{:});
