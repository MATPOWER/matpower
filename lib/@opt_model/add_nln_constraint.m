function om = add_nln_constraint(om, name, idx, N, iseq, fcn, hess, varsets)
%ADD_NLN_CONSTRAINT  Adds a set of nonlinear constraints to the model.
%
%   OM.ADD_NLN_CONSTRAINT(NAME, N, ISEQ, FCN, HESS);
%   OM.ADD_NLN_CONSTRAINT(NAME, N, ISEQ, FCN, HESS, VARSETS);
%
%   OM.ADD_NLN_CONSTRAINT(NAME, IDX_LIST, N, ISEQ, FCN, HESS);
%   OM.ADD_NLN_CONSTRAINT(NAME, IDX_LIST, N, ISEQ, FCN, HESS, VARSETS);
%
%   N specifies the number of constraints in the set, ISEQ is a one
%   character string specifying whether it is an equality or inequality
%   constraint set ('=' or '<' respectively), FCN is the handle of a
%   function that evaluates the constraint and its gradients, and HESS is
%   the handle of a function that evaluates the Hessian of the constraints.
%
%   For a constraint G(x) = 0, FCN should point to a function with the
%   following interface:
%       G = FCN(X)
%       [G, DG] = FCN(X)
%   where G is an N x 1 vector and DG is the N x NX gradient, where
%       DG(i, j) = dG(i)/dX(j)
%   and NX is the number of elements in X.
%
%   HESS should point to a function that returns an NX x NX matrix of
%   derivatives of DG * LAMBDA, with the following interface:
%       D2G = HESS(X, LAMBDA)
%
%   For both functions, the first input argument X can take two forms. If
%   the constraint set is added with VARSETS empty or missing, then X will
%   be the full optimization vector. Otherwise it will be a cell array of
%   vectors corresponding to the variable sets specified in VARSETS.
%
%   For simple (not indexed) named sets, NAME can be a cell array of
%   constraint set names, in which case N is a vector, specifying the number
%   of constraints in each corresponding set. FCN and HESS are still single
%   function handles for functions that compute the values for the entire
%   collection of constraint sets together.
%
%   Likewise, if FCN or HESS are empty, it also indicates a placeholder in
%   the indexing for a constraint set whose implmentation is included in
%   another constraint set. This functionality is only intended to be used
%   internally to handle constraint/gradient and Hessian functions that
%   compute the values for more than one constraint set simultaneously.
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
%       OM.ADD_NLN_CONSTRAINT(NAME, IDX_LIST, N, ISEQ, FCN, HESS);
%       OM.ADD_NLN_CONSTRAINT(NAME, IDX_LIST, N, ISEQ, FCN, HESS, VARSETS);
%
%   Examples:
%       %% nonlinear equality constraint with constraint/gradient and Hessian
%       %% evaluation functions provided
%       om.add_nln_constraint('Qmis', nb, 1, fcn, hess);
%
%       %% nonlinear inequality constraints with indexed named set 'S(i,j)'
%       om.init_indexed_name('nli', 'S', {2, 3});
%       for i = 1:2
%         for j = 1:3
%           om.add_nln_constraint('S', {i, j}, N{i,j}, ...);
%         end
%       end
%
%   See also OPT_MODEL, EVAL_NLN_CONSTRAINT.

%   MATPOWER
%   Copyright (c) 2008-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% initialize input arguments
if iscell(idx)          %% indexed named set
    if nargin < 8
        varsets = {};
    end
else                    %% simple named set
    if nargin < 7
        varsets = {};
    else
        varsets = hess;
    end
    if nargin > 4
        hess = fcn;
        fcn = iseq;
    end
    iseq = N;
    N = idx;
    idx = {};
end
if iseq
    ff = 'nle';     %% nonlinear equality
else
    ff = 'nli';     %% nonlinear inequality
end

%% convert varsets from cell to struct array if necessary
varsets = om.varsets_cell2struct(varsets);

%% add the named nonlinear constraint set
if iscell(name)
    if length(name) ~= length(N)
        error('@opt_model/add_nln_constraint: dimensions of NAME and N must match');
    end
    om.add_named_set(ff, name{1}, idx, N(1), fcn, hess, '', varsets);
    for k = 2:length(name)
        om.add_named_set(ff, name{k}, idx, N(k), [], [], name{1}, varsets);
    end
else
    om.add_named_set(ff, name, idx, N, fcn, hess, '', varsets);
end
