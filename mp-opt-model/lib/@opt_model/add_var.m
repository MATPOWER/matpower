function om = add_var(om, name, idx, varargin)
% add_var - Adds a set of variables to the model.
% ::
%
%   OM.ADD_VAR(NAME, N, V0, VL, VU, VT)
%   OM.ADD_VAR(NAME, N, V0, VL, VU)
%   OM.ADD_VAR(NAME, N, V0, VL)
%   OM.ADD_VAR(NAME, N, V0)
%   OM.ADD_VAR(NAME, N)
%   OM.ADD_VAR(NAME, DIM_LIST) (deprecated, use INIT_INDEXED_NAME instead)
%   OM.ADD_VAR(NAME, IDX_LIST, N, V0, VL, VU, VT)
%   OM.ADD_VAR(NAME, IDX_LIST, N, V0, VL, VU)
%   OM.ADD_VAR(NAME, IDX_LIST, N, V0, VL)
%   OM.ADD_VAR(NAME, IDX_LIST, N, V0)
%   OM.ADD_VAR(NAME, IDX_LIST, N)
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
%       om.add_var('V', nb, V0, Vmin, Vmax, 'C');
%
%       om.init_indexed_name('x', {2, 3});
%       for i = 1:2
%         for j = 1:3
%           om.add_var('x', {i, j}, nx(i,j), ...);
%         end
%       end
%
% See also opt_model, params_var.

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% Due to a bug related to inheritance in constructors in
%% Octave 5.2 and earlier (https://savannah.gnu.org/bugs/?52614),
%% INIT_SET_TYPES() cannot be called directly in the
%% MP_IDX_MANAGER constructor, as desired.
%%
%% WORKAROUND:  Initialize MP_IDX_MANAGER fields here, if needed,
%%              after object construction, but before object use.
if isempty(om.var)          %% only if not already initialized
    om.init_set_types();
end

%% set up default args
if iscell(idx) && isempty(varargin) %% just setting dimensions for indexed set
    om.init_indexed_name('var', name, idx);
else
    if iscell(idx)
        N = varargin{1};
        args = varargin(2:end);
    else
        N = idx;
        idx = {};
        args = varargin;
    end
    nargs = length(args);
    
    v0 = []; vl = []; vu = []; vt = [];
    if nargs >= 1
        v0 = args{1};
        if N > 1 && length(v0) == 1         %% expand from scalar as needed
            v0 = v0 * ones(N, 1);
        end
        if nargs >= 2
            vl = args{2};
            if N > 1 && length(vl) == 1     %% expand from scalar as needed
                vl = vl * ones(N, 1);
            end
            if nargs >= 3
                vu = args{3};
                if N > 1 && length(vu) == 1 %% expand from scalar as needed
                    vu = vu * ones(N, 1);
                end
                if nargs >= 4
                    vt = args{4};
                end
            end
        end
    end
    if isempty(v0)
        v0 = zeros(N, 1);   %% init to zero by default
    end
    if isempty(vl)
        vl = -Inf(N, 1);    %% unbounded below by default
    end
    if isempty(vu)
        vu = Inf(N, 1);     %% unbounded above by default
    end
    if isempty(vt) && N > 0
        vt = 'C';           %% all continuous by default
    end

    %% add the named variable set
    om.add_named_set('var', name, idx, N, v0, vl, vu, vt);
end
