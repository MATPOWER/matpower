function om = add_nln_cost(om, name, idx, varargin)
% add_nln_cost - Adds a set of general nonlinear costs to the model.
%
% .. note::
%    .. deprecated:: 5.0 Please use mp.sm_nln_cost.add instead, as
%       in ``om.nlc.add(...)``.
%
% ::
%
%   OM.ADD_NLN_COST(NAME, N, FCN);
%   OM.ADD_NLN_COST(NAME, N, FCN, VARSETS);
%   OM.ADD_NLN_COST(NAME, IDX_LIST, N, FCN);
%   OM.ADD_NLN_COST(NAME, IDX_LIST, N, FCN, VARSETS);
%
%   Adds a named block of general nonlinear costs to the model. FCN is
%   a handle to function that evaluates the cost, its gradient and Hessian
%   as described below.
%
%   The N parameter specifies the dimension for vector valued cost
%   functions, which are not yet implemented. Currently N must equal 1
%   or it will throw an error.
%
%   For a cost function f(x), FCN should point to a function with the
%   following interface:
%       F = FCN(X)
%       [F, DF] = FCN(X)
%       [F, DF, D2F] = FCN(X)
%   where F is a scalar with the value of the function, DF is the NX x 1
%   gradient, and D2F is the NX x NX Hessian and NX is the number of
%   elements in X.
%
%   The first input argument X can take two forms. If the constraint set
%   is added with VARSETS empty or missing, then X will be the full
%   optimization vector. Otherwise it will be a cell array of vectors
%   corresponding to the variable sets specified in VARSETS.
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
%       OM.ADD_NLN_COST(NAME, IDX_LIST, FCN);
%       OM.ADD_NLN_COST(NAME, IDX_LIST, FCN, VARSETS);
%
%   Examples:
%
%       fcn1 = @(x)my_cost_function1(x, other_args)
%       fcn2 = @(x)my_cost_function2(x, other_args)
%       om.add_nln_cost('mycost1', 1, fcn1);
%       om.add_nln_cost('mycost2', 1, fcn2, {'Vm', 'Pg', 'Qg', 'z'});
%
%       om.init_indexed_name('c', {2, 3});
%       for i = 1:2
%         for j = 1:3
%           om.add_nln_cost('c', {i, j}, 1, fcn(i,j), ...);
%         end
%       end
%
% See also opt_model, eval_nln_cost.

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

om.nlc.add(om.var, name, idx, varargin{:});
