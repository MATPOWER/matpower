function varargout = add_vars(om, varargin)
%ADD_VARS  Adds a set of variables to the model.
%
%   -----  DEPRECATED - use ADD_VAR instead  -----
%
%   OM.ADD_VARS(NAME, N, V0, VL, VU, VT)
%   OM.ADD_VARS(NAME, N, V0, VL, VU)
%   OM.ADD_VARS(NAME, N, V0, VL)
%   OM.ADD_VARS(NAME, N, V0)
%   OM.ADD_VARS(NAME, N)
%   OM.ADD_VARS(NAME, DIM_LIST) (deprecated, use INIT_INDEXED_NAME instead)
%   OM.ADD_VARS(NAME, IDX_LIST, N, V0, VL, VU, VT)
%   OM.ADD_VARS(NAME, IDX_LIST, N, V0, VL, VU)
%   OM.ADD_VARS(NAME, IDX_LIST, N, V0, VL)
%   OM.ADD_VARS(NAME, IDX_LIST, N, V0)
%   OM.ADD_VARS(NAME, IDX_LIST, N)
%   
%   Adds a set of variables to the model, where N is the number of
%   variables in the set, V0 is the initial value of those variables,
%   VL and VU are the lower and upper bounds on the variables and VT
%   is the variable type. The accepted values for elements of VT are:
%       'C' - continuous
%       'I' - integer
%       'B' - binary
%   V0, VL and VU are N x 1 column vectors, VT is a scalar or a 1 x N row
%   vector. The defaults for the last four arguments, which are all optional,
%   are for all values to be initialized to zero (V0 = 0), unbounded
%   (VL = -Inf, VU = Inf), and continuous (VT = 'C').
%
%   Examples:
%       om.add_vars('V', nb, V0, Vmin, Vmax, 'C');
%
%       om.init_indexed_name('x', {2, 3});
%       for i = 1:2
%         for j = 1:3
%           om.add_vars('x', {i, j}, nx(i,j), ...);
%         end
%       end
%
%   See also OPT_MODEL, PARAMS_VAR.

%   MATPOWER
%   Copyright (c) 2008-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

[varargout{1:nargout}] = om.add_var(varargin{:});
