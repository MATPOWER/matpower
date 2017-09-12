function om = add_lin_constraint(om, name, idx, A, l, u, varsets)
%ADD_LIN_CONSTRAINT  Adds a set of linear constraints to the model.
%
%   OM.ADD_LIN_CONSTRAINT(NAME, A, L, U);
%   OM.ADD_LIN_CONSTRAINT(NAME, A, L, U, VARSETS);
%
%   OM.ADD_LIN_CONSTRAINT(NAME, IDX_LIST, A, L, U);
%   OM.ADD_LIN_CONSTRAINT(NAME, IDX_LIST, A, L, U, VARSETS);
%
%   Linear constraints are of the form L <= A * x <= U, where x is a
%   vector made of the vars specified in VARSETS (in the order given).
%   This allows the A matrix to be defined only in terms of the relevant
%   variables without the need to manually create a lot of zero columns.
%   If VARSETS is empty, x is taken to be the full vector of all
%   optimization variables. If L or U are empty, they are assumed to be
%   appropriately sized vectors of -Inf and Inf, respectively.
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
%   See also OPT_MODEL, PARAMS_LIN_CONSTRAINT.

%   MATPOWER
%   Copyright (c) 2008-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% initialize input arguments
if iscell(idx)          %% indexed named set
    if nargin < 7
        varsets = {};
    end
else                    %% simple named set
    if nargin < 6
        varsets = {};
    else
        varsets = u;
    end
    u = l;
    l = A;
    A = idx;
    idx = {};
end

[N, M] = size(A);
if isempty(l)                   %% default l is -Inf
    l = -Inf(N, 1);
end
if isempty(u)                   %% default u is Inf
    u = Inf(N, 1);
end

%% check sizes
if size(l, 1) ~= N || size(u, 1) ~= N
    error('@opt_model/add_lin_constraint: sizes of A, l and u must match');
end

%% convert varsets from cell to struct array if necessary
varsets = om.varsets_cell2struct(varsets);
nv = om.varsets_len(varsets);   %% number of variables

%% check consistency of varsets with size of A
if M ~= nv
    error('@opt_model/add_lin_constraint: number of columns of A does not match\nnumber of variables, A is %d x %d, nv = %d\n', N, M, nv);
end

%% add the named linear constraint set
om.add_named_set('lin', name, idx, N, A, l, u, varsets);
