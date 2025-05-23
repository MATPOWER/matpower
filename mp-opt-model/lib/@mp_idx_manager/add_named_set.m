function obj = add_named_set(obj, set_type, name, idx, N, varargin)
% add_named_set - Adds a named set of a particular type to the object.
%
% .. note::
%    .. deprecated:: 5.0 Please use mp.set_manager.add instead.
%
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
%   Nonlinear Equality Constraint Set
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
    obj.(set_type).add(name, idx, N, varargin{:});
else
    ff = fieldnames(obj.set_types);
    stypes = sprintf('\n  ''%s''', ff{:});
    error('mp_idx_manager.add_named_set: ''%s'' is not a valid SET_TYPE, must be one of the following:%s', set_type, stypes);
end
