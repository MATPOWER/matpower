function obj = add_named_set(obj, set_type, name, idx, N, varargin)
% add_named_set - Adds a named set of a particular type to the object.
% ::
%
%   -----  PRIVATE METHOD  -----
%
%   This method is intended to be a private method, used internally by
%   the public methods for each set type, e.g. ADD_VAR, ADD_LIN_CONSTRAINT,
%   etc. This method handles the indexing part. The set-type-specific methods
%   that call it need to handle any data that goes with each set added.
%
%   E.g.
%
%   Variable Set
%       OBJ.ADD_NAMED_SET('var', NAME, IDX_LIST, N, V0, VL, VU, VT);
%
%   Linear Constraint Set
%       OBJ.ADD_NAMED_SET('lin', NAME, IDX_LIST, N, A, L, U, VARSETS);
%
%   Nonlinear Inequality Constraint Set
%       OBJ.ADD_NAMED_SET('nle', NAME, IDX_LIST, N, FCN, HESS, COMPUTED_BY, VARSETS);
%
%   Nonlinear Inequality Constraint Set
%       OBJ.ADD_NAMED_SET('nli', NAME, IDX_LIST, N, FCN, HESS, COMPUTED_BY, VARSETS);
%
%   Quadratic Cost Set
%       OBJ.ADD_NAMED_SET('qdc', NAME, IDX_LIST, N, CP, VARSETS);
%
%   General Nonlinear Cost Set
%       OBJ.ADD_NAMED_SET('nlc', NAME, IDX_LIST, N, FCN, VARSETS);
%
%   See OPT_MODEL and its methods ADD_VAR, ADD_LIN_CONSTRAINT,
%       ADD_NLN_CONSTRAINT, ADD_QUAD_COST and ADD_NLN_COST.
%
% See also opt_model.

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
    obj_ff = obj.(ff);
    obj.(ff) = [];
else
    ff = fieldnames(obj.set_types);
    stypes = sprintf('\n  ''%s''', ff{:});
    error('mp_idx_manager.add_named_set: ''%s'' is not a valid SET_TYPE, must be one of the following:%s', set_type, stypes);
end

%% add general indexing info about this named set
if isempty(idx)     %% simple named set
    %% prevent duplicate name in set of specified type
    if isfield(obj_ff.idx.N, name)
        error('mp_idx_manager.add_named_set: %s set named ''%s'' already exists', st_label, name);
    end

    %% add indexing info about this set
    obj_ff.idx.i1.(name)  = obj_ff.N + 1;   %% starting index
    obj_ff.idx.iN.(name)  = obj_ff.N + N;   %% ending index
    obj_ff.idx.N.(name)   = N;              %% number in set
    obj_ff.N  = obj_ff.idx.iN.(name);       %% number of elements of this type
    obj_ff.NS = obj_ff.NS + 1;              %% number of sets of this type
    obj_ff.order(obj_ff.NS).name = name;    %% add name to ordered list of sets
    obj_ff.order(obj_ff.NS).idx  = {};
else                %% indexed named set
    %% calls to substruct() are relatively expensive, so we pre-build the
    %% struct for addressing numeric array fields
    %% sn = substruct('.', name, '()', idx);
    sn = struct('type', {'.', '()'}, 'subs', {name, idx});  %% num array field

    %% prevent duplicate name in set of specified type
    if subsref(obj_ff.idx.i1, sn) ~= 0
        str = '%d'; for m = 2:length(idx), str = [str ',%d']; end
        nname = sprintf(['%s(' str, ')'], name, idx{:});
        error('mp_idx_manager.add_named_set: %s set named ''%s'' already exists', st_label, nname);
    end

    %% add indexing info about this set
    obj_ff.idx.i1  = subsasgn(obj_ff.idx.i1, sn, obj_ff.N + 1); %% starting index
    obj_ff.idx.iN  = subsasgn(obj_ff.idx.iN, sn, obj_ff.N + N); %% ending index
    obj_ff.idx.N   = subsasgn(obj_ff.idx.N,  sn, N);            %% number in set
    obj_ff.N  = subsref(obj_ff.idx.iN, sn); %% number of elements of this type
    obj_ff.NS = obj_ff.NS + 1;              %% number of sets of this type
    obj_ff.order(obj_ff.NS).name = name;    %% add name to ordered list of sets
    obj_ff.order(obj_ff.NS).idx  = idx;     %% add indices to ordered list of sets
end

obj.(ff) = obj_ff;
