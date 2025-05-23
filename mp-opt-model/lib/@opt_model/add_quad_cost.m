function om = add_quad_cost(om, name, idx, varargin)
% add_quad_cost - Adds a set of user costs to the model.
%
% .. note::
%    .. deprecated:: 5.0 Please use mp.sm_quad_cost.add instead, as
%       in ``om.qdc.add(...)``.
%
% ::
%
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
%   by default). Alternatively, if Q is an NX x 1 vector or empty, then F(X)
%   is also NX x 1, and K can be either NX x 1 or scalar.
%       F(X) = 1/2 * Q .* X.^2 + C .* X + K
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
% See also opt_model, params_quad_cost, eval_quad_cost.

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

om.qdc.add(om.var, name, idx, varargin{:});
