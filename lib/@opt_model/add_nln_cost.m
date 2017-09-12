function om = add_nln_cost(om, name, idx, N, fcn, varsets)
%ADD_NLN_COST  Adds a set of general nonlinear costs to the model.
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
%   where F is a scalar with the value of the function, DG is the 1 x NX
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
%   See also OPT_MODEL, EVAL_NLN_COST.

%   MATPOWER
%   Copyright (c) 2008-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% initialize input arguments
if iscell(idx)          %% indexed named set
    if nargin < 6
        varsets = {};
    end
else                    %% simple named set
    if nargin < 5
        varsets = {};
    else
        varsets = fcn;
    end
    fcn = N;
    N = idx;
    idx = {};
end

if N ~= 1
    error('@opt_model/add_nln_cost: not yet implemented for vector valued functions (i.e. N currently must equal 1)');
end

%% convert varsets from cell to struct array if necessary
varsets = om.varsets_cell2struct(varsets);

%% add the named general nonlinear cost set
om.add_named_set('nlc', name, idx, N, fcn, varsets);
