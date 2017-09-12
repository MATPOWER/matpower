function om = add_quad_cost(om, name, idx, Q, c, k, varsets)
%ADD_QUAD_COST  Adds a set of user costs to the model.
%   OM.ADD_QUAD_COST(NAME, Q, C);
%   OM.ADD_QUAD_COST(NAME, Q, C, K);
%   OM.ADD_QUAD_COST(NAME, Q, C, K, VARSETS);
%   OM.ADD_QUAD_COST(NAME, IDX_LIST, Q, C);
%   OM.ADD_QUAD_COST(NAME, IDX_LIST, Q, C, K);
%   OM.ADD_QUAD_COST(NAME, IDX_LIST, Q, C, K, VARSETS);
%
%   Adds a named block of quadratic costs to the model. Costs are of the
%   form
%       F(X) = 1/2 * X'*Q*X + C'*X + K
%   where Q is an NX x NX matrix (possibly sparse), C is an NX x 1 vector,
%   K is a scalar and NX is the number of elements in X. Here X is the vector
%   formed by combining the specified VARSETS (the full optimization vector
%   by default). Alternatively, if Q is an NX x 1 vector, then K and F(X)
%   are also NX x 1 vectors
%       F(X) = 1/2 * Q .* X.^2 + C .* X + K
%
%   If Q is empty, then F(X) is also an NX x 1 vector, unless K is scalar
%   and non-zero.
%
%   Indexed Named Sets
%       A cost set can be identified by a single NAME, as described
%       above, such as 'PgCost', or by a name that is indexed by one
%       or more indices, such as 'PgCost(3,4)'. For an indexed named
%       set, before adding the cost sets themselves, the dimensions
%       of the indexed set must be set by calling INIT_INDEXED_NAME.
%
%       The constraints are then added using the following, where
%       all arguments are as described above, except IDX_LIST is a cell
%       array of the indices for the particular cost set being added.
%
%       OM.ADD_QUAD_COST(NAME, IDX_LIST, Q, C, K);
%       OM.ADD_QUAD_COST(NAME, IDX_LIST, Q, C, K, VARSETS);
%
%   Examples:
%       om.add_quad_cost('quad_cost1', Q1, c1, 0);
%       om.add_quad_cost('lin_cost2',  [], c2, k2, {'Vm', 'Pg', 'z'});
%
%       om.init_indexed_name('c', {2, 3});
%       for i = 1:2
%         for j = 1:3
%           om.add_quad_cost('c', {i, j}, Q{i,j}, ...);
%         end
%       end
%
%   See also OPT_MODEL, PARAMS_QUAD_COST, EVAL_QUAD_COST.

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
        varsets = k;
    end
    if nargin < 5
        k = 0;
    else
        k = c;
    end
    c = Q;
    Q = idx;
    idx = {};
end

%% convert varsets from cell to struct array if necessary
varsets = om.varsets_cell2struct(varsets);
nv = om.varsets_len(varsets);   %% number of variables

%% check sizes
[MQ, NQ] = size(Q);
[Mc, Nc] = size(c);
if MQ
    if NQ ~= MQ && NQ ~= 1
        error('@opt_model/add_quad_cost: Q (%d x %d) must be square or a column vector (or empty)', MQ, NQ);
    end
end
if Mc && Nc ~= 1
    error('@opt_model/add_quad_cost: c (%d x %d) must be a column vector (or empty)', Mc, Nc);
end
if MQ
    if Mc && Mc ~= MQ
        error('@opt_model/add_quad_cost: dimensions of Q (%d x %d) and c (%d x %d) are not compatible', MQ, NQ, Mc, Nc);
    end
    nx = MQ;
else
    if ~Mc
        error('@opt_model/add_quad_cost: Q and c cannot both be empty');
    end
    nx = Mc;
end
if nx ~= nv
    error('@opt_model/add_quad_cost: dimensions of Q (%d x %d) and c (%d x %d) do not match\nnumber of variables (%d)\n', MQ, NQ, Mc, Nc, nv);
end

%% size of named cost set
if NQ == 1 || (isempty(Q) && (length(k) > 1 || k == 0))
    %% Q is a column vector (cost is element-wise, i.e. a vector)
    %% OR Q is empty and k is either a vector or zero
    N = nx;
else            %% Q is a square matrix (cost is a scalar)
    N = 1;
end

%% add the named quadratic cost set
om.add_named_set('qdc', name, idx, N, Q, c, k, varsets);
