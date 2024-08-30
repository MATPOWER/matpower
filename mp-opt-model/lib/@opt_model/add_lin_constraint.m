function om = add_lin_constraint(om, name, idx, varargin)
% add_lin_constraint - Adds a set of linear constraints to the model.
%
% .. note::
%    .. deprecated:: 4.3 Please use mp.sm_lin_constraint.add instead, as
%       in ``om.lin.add(...)``.
%
% ::
%
%   OM.ADD_LIN_CONSTRAINT(NAME, A, L, U);
%   OM.ADD_LIN_CONSTRAINT(NAME, A, L, U, VARSETS);
%   OM.ADD_LIN_CONSTRAINT(NAME, A, L, U, VARSETS, TR);
%
%   OM.ADD_LIN_CONSTRAINT(NAME, IDX_LIST, A, L, U);
%   OM.ADD_LIN_CONSTRAINT(NAME, IDX_LIST, A, L, U, VARSETS);
%   OM.ADD_LIN_CONSTRAINT(NAME, IDX_LIST, A, L, U, VARSETS, TR);
%
%   Linear constraints are of the form L <= A * x <= U, where x is a
%   vector made of the vars specified in VARSETS (in the order given).
%   This allows the A matrix to be defined in terms of only the relevant
%   variables without the need to manually create a lot of zero columns.
%   If VARSETS is empty, x is taken to be the full vector of all
%   optimization variables. If L or U are empty, they are assumed to be
%   appropriately sized vectors of -Inf and Inf, respectively. If TR is
%   present and true, it means that A' is supplied/stored rather than A.
%   In some contexts this can be used to save significant memory.
%
%   Indexed Named Sets
%       A constraint set can be identified by a single NAME, as described
%       above, such as 'Pmismatch', or by a name that is indexed by one
%       or more indices, such as 'Pmismatch(3,4)'. For an indexed named
%       set, before adding the constraint sets themselves, the dimensions
%       of the indexed set must be set by calling INIT_INDEXED_NAME.
%
%       The constraints are then added using the following, where
%       all arguments are as described above, except IDX_LIST is a cell
%       array of the indices for the particular constraint set being added.
%
%       OM.ADD_LIN_CONSTRAINT(NAME, IDX_LIST, A, L, U);
%       OM.ADD_LIN_CONSTRAINT(NAME, IDX_LIST, A, L, U, VARSETS);
%
%   Examples:
%       %% linear constraint
%       om.add_lin_constraint('vl', Avl, lvl, uvl, {'Pg', 'Qg'});
%
%       %% linear constraints with indexed named set 'R(i,j)'
%       om.init_indexed_name('lin', 'R', {2, 3});
%       for i = 1:2
%         for j = 1:3
%           om.add_lin_constraint('R', {i, j}, A{i,j}, ...);
%         end
%       end
%
% See also opt_model, params_lin_constraint.

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

om.lin.add(om.var, name, idx, varargin{:});
