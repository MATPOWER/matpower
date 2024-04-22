function obj = init_indexed_name(obj, set_type, name, dim_list)
% init_indexed_name - Initializes the dimensions for an indexed named set.
% ::
%
%   OBJ.INIT_INDEXED_NAME(SET_TYPE, NAME, DIM_LIST)
%
%   Initializes the dimensions for an indexed named set.
%
%   For example, variables and linear constraints could be referenced in
%   terms of named sets, where the type of named set being referenced is
%   given by SET_TYPE, with the following valid options:
%       SET_TYPE = 'var'   => variable set
%       SET_TYPE = 'lin'   => linear constraint set
%
%   Indexed Named Sets
%
%   In this case a variable or constraint set can be identified by a single
%   NAME, such as 'Pmismatch', or by a name that is indexed by one or more
%   indices, such as 'Pmismatch(3,4)'. For an indexed named set, before adding
%   the indexed variable or constraint sets themselves, the dimensions of
%   the indexed set must be set by calling INIT_INDEXED_NAME, where DIM_LIST
%   is a cell array of the dimensions.
%
%   Examples:
%       %% linear constraints with indexed named set 'R(i,j)'
%       obj.init_indexed_name('lin', 'R', {2, 3});
%       for i = 1:2
%         for j = 1:3
%           obj.add_lin_constraint('R', {i, j}, A{i,j}, ...);
%         end
%       end
%
%   See OPT_MODEL, and its methods INIT_INDEXED_NAME, ADD_VAR,
%           ADD_LIN_CONSTRAINT, ADD_NLN_CONSTRAINT, ADD_QUAD_COST and
%           ADD_NLN_COST.

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% check for valid type for named set
st_label = obj.valid_named_set_type(set_type);
if st_label
    ff = set_type;
else
    ff = fieldnames(obj.set_types);
    stypes = sprintf('\n  ''%s''', ff{:});
    error('mp_idx_manager.init_indexed_name: ''%s'' is not a valid SET_TYPE, must be one of the following:%s', set_type, stypes);
end

%% prevent duplicate name in set of specified type
if isfield(obj.(ff).idx.N, name)
    error('mp_idx_manager.init_indexed_name: %s set named ''%s'' already exists', ...
        st_label, name);
end

%% use column vector if single dimension
if length(dim_list) == 1
    dim_list = {dim_list{:}, 1};
end

%% add general info about this named set
zero_vector = zeros(dim_list{:});
obj.(ff).idx.i1.(name)  = zero_vector;  %% starting index
obj.(ff).idx.iN.(name)  = zero_vector;  %% ending index
obj.(ff).idx.N.(name)   = zero_vector;  %% number of vars/constraints/costs
