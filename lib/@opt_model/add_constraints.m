function om = add_constraints(om, name, varargin)
%ADD_CONSTRAINTS  Adds a set of constraints to the model.
%
%   -----  DEPRECATED - Please use one of the following instead:      -----
%   -----  ADD_LIN_CONSTRAINT, ADD_NLN_CONSTRAINT, INIT_INDEXED_NAME  -----
%
%   OM.ADD_CONSTRAINTS(NAME, A, L, U);
%   OM.ADD_CONSTRAINTS(NAME, A, L, U, VARSETS);
%
%   OM.ADD_CONSTRAINTS(NAME, N, ISEQ, FCN, HESS);
%   OM.ADD_CONSTRAINTS(NAME, N, ISEQ, FCN, HESS, VARSETS);
%
%   OM.ADD_CONSTRAINTS(NAME, DIM_LIST);           %% linear
%   OM.ADD_CONSTRAINTS(NAME, DIM_LIST, 'lin');    %% linear
%   OM.ADD_CONSTRAINTS(NAME, DIM_LIST, 'nle');    %% nonlinear equality
%   OM.ADD_CONSTRAINTS(NAME, DIM_LIST, 'nli');    %% nonlinear inequality
%
%   OM.ADD_CONSTRAINTS(NAME, IDX_LIST, A, L, U);
%   OM.ADD_CONSTRAINTS(NAME, IDX_LIST, A, L, U, VARSETS);
%
%   OM.ADD_CONSTRAINTS(NAME, IDX_LIST, N, ISEQ, FCN, HESS);
%   OM.ADD_CONSTRAINTS(NAME, IDX_LIST, N, ISEQ, FCN, HESS, VARSETS);
%
%   Linear Constraints
%
%       OM.ADD_CONSTRAINTS(NAME, A, L, U);
%       OM.ADD_CONSTRAINTS(NAME, A, L, U, VARSETS);
%
%       Linear constraints are of the form L <= A * x <= U, where
%       x is a vector made of the vars specified in VARSETS (in
%       the order given). This allows the A matrix to be defined only
%       in terms of the relevant variables without the need to manually
%       create a lot of zero columns. If VARSETS is empty, x is taken
%       to be the full vector of all optimization variables. If L or 
%       U are empty, they are assumed to be appropriately sized vectors
%       of -Inf and Inf, respectively.
%
%   Nonlinear Constraints
%
%       OM.ADD_CONSTRAINTS(NAME, N, ISEQ, FCN, HESS);
%       OM.ADD_CONSTRAINTS(NAME, N, ISEQ, FCN, HESS, VARSETS);
%
%       For nonlinear constraints, N specifies the number of constraints
%       in the set, ISEQ is a one character string specifying whether it
%       is an equality or inequality constraint set ('=' or '<' respectively),
%       FCN is the handle of a function that evaluates the constraint
%       and its gradients, and HESS is the handle of a function that
%       evaluates the Hessian of the constraints.
%
%       For a constraint G(x) = 0, FCN should point to a function with the
%       following interface:
%           G = FCN(X)
%           [G, DG] = FCN(X)
%       where G is an N x 1 vector and DG is the N x NX gradient, where
%           DG(i, j) = dG(i)/dX(j)
%       and NX is the number of elements in X.
%
%       HESS should point to a function that returns an NX x NX matrix
%       of derivatives of DG * LAMBDA, with the following interface:
%           D2G = HESS(X, LAMBDA)
%
%       For both functions, the first input argument X can take two forms.
%       If the constraint set is added VARSETS empty or missing, then X
%       will be the full optimization vector. Otherwise it will be a cell
%       array of vectors corresponding to the variable sets specified in
%       VARSETS.
%
%       For nonlinear constraints, NAME can be a cell array of constraint
%       set names, in which case N is a vector, specifying the number of
%       constraints in each set. FCN and HESS are still single function
%       handles for functions that compute the values for the entire
%       collection of constraint sets together.
%
%       Similarly, if the third argument is a string containing 'nle' or
%       'nli', it indicates a placeholder in the indexing for a constraint
%       set whose implmentation is included in another constraint set. This
%       functionality is only intended to be used internally to handle
%       constraint/gradient and Hessian functions that handle more than one
%       constraint set simultaneously.
%
%   Indexed Named Sets
%       A constraint set can be identified by a single NAME, as described
%       above, such as 'Pmismatch', or by a name that is indexed by one
%       or more indices, such as 'Pmismatch(3,4)'. For an indexed named
%       set, before adding the constraint sets themselves, the dimensions
%       of the indexed set must be set by calling one of the following,
%       where DIM_LIST is a cell array of the dimensions, and the 4th
%       argument indicates the type of constraint set (linear by default).
%
%       OM.ADD_CONSTRAINTS(NAME, DIM_LIST);           %% linear
%       OM.ADD_CONSTRAINTS(NAME, DIM_LIST, 'lin');    %% linear
%       OM.ADD_CONSTRAINTS(NAME, DIM_LIST, 'nle');    %% nonlin equality
%       OM.ADD_CONSTRAINTS(NAME, DIM_LIST, 'nli');    %% nonlin inequality
%
%       The constraints are then added using one of the following, where
%       all arguments are as described above, except IDX_LIST is a cell
%       array of the indices for the particular constraint set being added.
%
%       OM.ADD_CONSTRAINTS(NAME, IDX_LIST, A, L, U);
%       OM.ADD_CONSTRAINTS(NAME, IDX_LIST, A, L, U, VARSETS);
%       OM.ADD_CONSTRAINTS(NAME, IDX_LIST, N, ISEQ, FCN, HESS);
%       OM.ADD_CONSTRAINTS(NAME, IDX_LIST, N, ISEQ, FCN, HESS, VARSETS);
%
%   Examples:
%       %% linear constraint
%       om.add_constraints('vl', Avl, lvl, uvl, {'Pg', 'Qg'});
%
%       %% nonlinear equality constraint with constraint/gradient and Hessian
%       %% evaluation functions provided
%       om.add_constraints('Qmis', nb, 1, fcn, hess);
%
%       %% linear constraints with indexed named set 'R(i,j)'
%       om.add_constraints('R', {2, 3}, 'lin');
%       for i = 1:2
%         for j = 1:3
%           om.add_constraints('R', {i, j}, A{i,j}, ...);
%         end
%       end
%
%   See also OPT_MODEL, LINEAR_CONSTRAINTS, EVAL_NLN_CONSTRAINT.

%   MATPOWER
%   Copyright (c) 2008-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% determine whether it's linear or nonlinear, simple or indexed
%% and if indexed, just setting dimensions or actually adding constraints
%%      linear    : nonlin = 0, ff = 'lin'
%%      nonlinear : nonlin = 1
%%          equality :          ff = 'nle'
%%          inequality :        ff = 'nli'
%%      simple    : idx empty, args not empty
%%      indexed   : idx not empty
%%          set dims only : args empty
%%          add idx'd cons: args not empty
nonlin = 0;
ff = 'lin';
if iscell(varargin{1})          %% indexed named set
    idx = varargin{1};          %% dimensions or indices
    if length(varargin) < 3     %% set dimensions for indexed set
        args = {};
        if length(varargin) > 1
            ff = lower(varargin{2});
            if strcmp(ff, 'nle')
                nonlin = 1;
            elseif strcmp(ff, 'nli')
                nonlin = 1;
            elseif ~strcmp(ff, 'lin')
                error('@opt_model/add_constraints: ''%s'' is not a valid type (''lin'', ''nle'', or ''nli'') for an indexed constraint set', ff);
            end
        end
    else                        %% indexed set
        args = varargin(2:end);
        if isa(args{3}, 'function_handle')
            if args{2}
                ff = 'nle';
            else
                ff = 'nli';
            end
            nonlin = 1;
        end
    end
else                            %% simple named set
    idx = {};
    args = varargin;

    if length(args) < 3 || isa(args{3}, 'function_handle')
        if args{2}
            ff = 'nle';
        else
            ff = 'nli';
        end
        nonlin = 1;
    end
end
nargs = length(args);

if ~isempty(idx) && nargs == 0      %% just set dimensions for indexed set
    om.init_indexed_name(ff, name, idx);
else
    if nonlin           %% nonlinear
        if nargs == 2
            [N, iseq, fcn, hess, varsets] = deal(args{:}, [], [], {});
        elseif nargs < 5
            [N, iseq, fcn, hess, varsets] = deal(args{1:4}, {});
        else
            [N, iseq, fcn, hess, varsets] = deal(args{1:5});
        end
        om.add_nln_constraint(name, idx, N, iseq, fcn, hess, varsets);
    else                %% linear
        if nargs < 4
            [A, l, u, varsets] = deal(args{1:3}, {});
        else
            [A, l, u, varsets] = deal(args{1:4});
        end
        om.add_lin_constraint(name, idx, A, l, u, varsets);
    end
end
