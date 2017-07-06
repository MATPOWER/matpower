function om = add_constraints(om, name, varargin)
%ADD_CONSTRAINTS  Adds a set of constraints to the model.
%
%   OM.ADD_CONSTRAINTS(NAME, A, L, U);
%   OM.ADD_CONSTRAINTS(NAME, A, L, U, VARSETS);
%
%   OM.ADD_CONSTRAINTS(NAME, N, ISEQ, G_DG, D2G);
%   OM.ADD_CONSTRAINTS(NAME, N, ISEQ, G_DG, D2G, VARSETS);
%   OM.ADD_CONSTRAINTS(NAME, N, 'nonlinear');  (legacy)
%
%   OM.ADD_CONSTRAINTS(NAME, DIM_LIST);           %% linear
%   OM.ADD_CONSTRAINTS(NAME, DIM_LIST, 'lin');    %% linear
%   OM.ADD_CONSTRAINTS(NAME, DIM_LIST, 'nle');    %% nonlinear equality
%   OM.ADD_CONSTRAINTS(NAME, DIM_LIST, 'nli');    %% nonlinear inequality
%
%   OM.ADD_CONSTRAINTS(NAME, IDX_LIST, A, L, U);
%   OM.ADD_CONSTRAINTS(NAME, IDX_LIST, A, L, U, VARSETS);
%
%   OM.ADD_CONSTRAINTS(NAME, IDX_LIST, N, ISEQ, G_DG, D2G);
%   OM.ADD_CONSTRAINTS(NAME, IDX_LIST, N, ISEQ, G_DG, D2G, VARSETS);
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
%       OM.ADD_CONSTRAINTS(NAME, N, ISEQ, G_DG, D2G);
%       OM.ADD_CONSTRAINTS(NAME, N, ISEQ, G_DG, D2G, VARSETS);
%       OM.ADD_CONSTRAINTS(NAME, N, 'nonlinear');  (legacy)
%
%       For nonlinear constraints, N specifies the number of constraints
%       in the set, ISEQ is a one character string specifying whether it
%       is an equality or inequality constraint set ('=' or '<' respectively),
%       G_DG is the handle of a function that evaluates the constraint
%       and its gradients, and D2G is the handle of a function that
%       evaluates the Hessian of the constraints.
%
%       For a constraint G(x) = 0, G_DG should point to a function with the
%       following interface:
%           G = G_DG(X)
%           [G, DG] = G_DG(X)
%       where G is an N x 1 vector and DG is the N x NX gradient, where
%           DG(i, j) = dG(i)/dX(j)
%       and NX is the number of elements in X.
%
%       D2G should point to a function that returns an NX x NX matrix
%       of derivatives of DG * LAMBDA, with the following interface:
%           D2G = D2G(X, LAMBDA)
%
%       For both functions, the first input argument X can take two forms.
%       If the constraint set is added VARSETS empty or missing, then X
%       will be the full optimization vector. Otherwise it will be a cell
%       array of vectors corresponding to the variable sets specified in
%       VARSETS.
%
%       For backward compatibility, ADD_CONSTRAINTS can also simply include
%       the number of constraints N followed by the string 'nonlinear' to
%       indicate a placeholder in the indexing for constraints that will
%       be implemented outside of the OPT_MODEL.
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
%       OM.ADD_CONSTRAINTS(NAME, IDX_LIST, N, ISEQ, G_DG, D2G);
%       OM.ADD_CONSTRAINTS(NAME, IDX_LIST, N, ISEQ, G_DG, D2G, VARSETS);
%
%   Examples:
%       %% linear constraint
%       om.add_constraints('vl', Avl, lvl, uvl, {'Pg', 'Qg'});
%
%       %% legacy nonlinear constraint (indexing placeholder)
%       om.add_constraints('Pmis', nb, 'nonlinear');
%
%       %% nonlinear equality constraint with constraint/gradient and Hessian
%       %% evaluation functions provided
%       om.add_constraints('Qmis', nb, 1, g_dg, d2G);
%
%       %% linear constraints with indexed named set 'R(i,j)'
%       om.add_constraints('R', {2, 3}, 'lin');
%       for i = 1:2
%         for j = 1:3
%           om.add_constraints('R', {i, j}, A{i,j}, ...);
%         end
%       end
%
%   See also OPT_MODEL, LINEAR_CONSTRAINTS, NONLIN_CONSTRAINTS.

%   MATPOWER
%   Copyright (c) 2008-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% determine whether it's linear or nonlinear, simple or indexed
%% and if indexed, just setting dimensions or actually adding constraints
%%      linear    : nonlin = 0, ff = 'lin', label = 'linear'
%%      nonlinear : nonlin = 1
%%          equality :          ff = 'nle', label = 'nonlinear equality'
%%          inequality :        ff = 'nli', label = 'nonlinear inequality'
%%      simple    : idx empty, args not empty
%%      indexed   : idx not empty
%%          set dims only : args empty
%%          add idx'd cons: args not empty
nonlin = 0;
label = 'linear';
ff = 'lin';
if iscell(varargin{1})          %% indexed named set
    idx = varargin{1};          %% dimensions or indices
    if length(varargin) < 3     %% set dimensions for indexed set
        args = {};
        if length(varargin) > 1
            ff = lower(varargin{2});
            if strcmp(ff, 'nle')
                label = 'nonlinear equality';
                nonlin = 1;
            elseif strcmp(ff, 'nli')
                label = 'nonlinear inequality';
                nonlin = 1;
            elseif ~strcmp(ff, 'lin')
                error('@opt_model/add_constraints: ''%s'' is not a valid type (''lin'', ''nle'', or ''nli'') for an indexed constraint set', ff);
            end
        end
        %% prevent duplicate named constraint sets
        if isfield(om.(ff).idx.N, name)
            error('@opt_model/add_constraints: %s constraint set named ''%s'' already exists', label, name);
        end
    else                        %% indexed set
        args = varargin(2:end);

        % (calls to substruct() are relatively expensive ...
        % s1 = substruct('.', name, '()', idx);
        % s2 = substruct('.', name, '{}', idx);
        % ... so replace them with these more efficient lines)
        s1 = struct('type', {'.', '()'}, 'subs', {name, idx});
        s2 = s1;
        s2(2).type = '{}';

        %% prevent duplicate named constraint sets
        if isa(args{3}, 'function_handle')
            if args{2}
                label = 'nonlinear equality';
                ff = 'nle';
            else
                label = 'nonlinear inequality';
                ff = 'nli';
            end
            nonlin = 1;
        end
        if subsref(om.(ff).idx.i1, s1) ~= 0
            str = '%d'; for m = 2:length(idx), str = [str ',%d']; end
            nname = sprintf(['%s(' str, ')'], name, idx{:});
            error('@opt_model/add_constraints: %s constraint set named ''%s'' already exists', label, nname);
        end
    end
else                            %% simple named set
    idx = {};
    args = varargin;

    if length(args) < 3 || isa(args{3}, 'function_handle')
        if length(args) < 3         %% legacy
            label = 'nonlinear';
            ff = 'nln';
        elseif args{2}
            label = 'nonlinear equality';
            ff = 'nle';
        else
            label = 'nonlinear inequality';
            ff = 'nli';
        end
        nonlin = 1;
    end

    %% prevent duplicate named constraint sets
    if isfield(om.(ff).idx.N, name)
        error('@opt_model/add_constraints: %s constraint set named ''%s'' already exists', label, name);
    end
end
nargs = length(args);

if ~isempty(idx) && nargs == 0      %% just set dimensions for indexed set
    %% use column vector if single dimension
    if length(idx) == 1
        idx = {idx{:}, 1};
    end

    %% add info about this constraint set
    om.(ff).idx.i1.(name)    = zeros(idx{:});   %% starting index
    om.(ff).idx.iN.(name)    = zeros(idx{:});   %% ending index
    om.(ff).idx.N.(name)     = zeros(idx{:});   %% number of constraints
    om.(ff).data.vs.(name)   = cell(idx{:});
    if nonlin
        %% add info about this nonlinear constraint set
        om.(ff).data.g_dg.(name) = cell(idx{:});
        om.(ff).data.d2G.(name)  = cell(idx{:});
    else
        %% add info about this linear constraint set
        om.(ff).data.A.(name)    = cell(idx{:});
        om.(ff).data.l.(name)    = cell(idx{:});
        om.(ff).data.u.(name)    = cell(idx{:});
    end
else
    if nonlin           %% nonlinear
        if nargs == 2
            [N, iseq, g_dg, d2G, varsets] = deal(args{:}, [], [], {});
        elseif nargs < 5
            [N, iseq, g_dg, d2G, varsets] = deal(args{1:4}, {});
        else
            [N, iseq, g_dg, d2G, varsets] = deal(args{1:5});
        end
    else                %% linear
        if nargs < 4
            [A, l, u, varsets] = deal(args{1:3}, {});
        else
            [A, l, u, varsets] = deal(args{1:4});
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
            error('@opt_model/add_constraints: sizes of A, l and u must match');
        end
    end

    if ~isempty(varsets) && iscell(varsets)
        empty_cells = cell(1, length(varsets));
        [empty_cells{:}] = deal({});    %% empty cell arrays
        varsets = struct('name', varsets, 'idx', empty_cells);
    end
    if isempty(varsets)
        nv = om.var.N;
    else
        nv = 0;
        s = struct('type', {'.', '()'}, 'subs', {'', 1});
        for k = 1:length(varsets)
            % (calls to substruct() are relatively expensive ...
            % s = substruct('.', varsets(k).name, '()', varsets(k).idx);
            % ... so replace it with these more efficient lines)
            s(1).subs = varsets(k).name;
            s(2).subs = varsets(k).idx;
            nv = nv + subsref(om.var.idx.N, s);
        end
    end
    if ~nonlin && M ~= nv
        error('@opt_model/add_constraints: number of columns of A does not match\nnumber of variables, A is %d x %d, nv = %d\n', N, M, nv);
    end

    if isempty(idx)     %% simple
        %% add info about this constraint set
        om.(ff).idx.i1.(name)  = om.(ff).N + 1; %% starting index
        om.(ff).idx.iN.(name)  = om.(ff).N + N; %% ending index
        om.(ff).idx.N.(name)   = N;             %% number of constraints
        om.(ff).data.vs.(name) = varsets;

        %% update number of constraints and constraint sets
        om.(ff).N  = om.(ff).idx.iN.(name);
        om.(ff).NS = om.(ff).NS + 1;

        %% add to ordered list of constraint sets
        om.(ff).order(om.(ff).NS).name = name;
        om.(ff).order(om.(ff).NS).idx  = {};

        if nonlin
            if ~isempty(g_dg)
                om.(ff).data.g_dg.(name) = g_dg;
                om.(ff).data.d2G.(name)  = d2G;
            end
        else
            om.(ff).data.A.(name)  = A;
            om.(ff).data.l.(name)  = l;
            om.(ff).data.u.(name)  = u;
        end
    else                %% indexed
        %% add info about this constraint set
        om.(ff).idx.i1  = subsasgn(om.(ff).idx.i1, s1, om.(ff).N + 1);  %% starting index
        om.(ff).idx.iN  = subsasgn(om.(ff).idx.iN, s1, om.(ff).N + N);  %% ending index
        om.(ff).idx.N   = subsasgn(om.(ff).idx.N,  s1, N);              %% number of constraints
        om.(ff).data.vs = subsasgn(om.(ff).data.vs, s2, varsets);

        %% update number of constraints and constraint sets
        om.(ff).N  = subsref(om.(ff).idx.iN, s1);
        om.(ff).NS = om.(ff).NS + 1;

        %% add to ordered list of constraint sets
        om.(ff).order(om.(ff).NS).name = name;
        om.(ff).order(om.(ff).NS).idx  = idx;

        if nonlin
            %% add info about this nonlinear constraint set
            om.(ff).data.g_dg   = subsasgn(om.(ff).data.g_dg, s2, g_dg);
            om.(ff).data.d2G    = subsasgn(om.(ff).data.d2G, s2, d2G);
        else
            om.(ff).data.A  = subsasgn(om.(ff).data.A, s2, A);
            om.(ff).data.l  = subsasgn(om.(ff).data.l, s2, l);
            om.(ff).data.u  = subsasgn(om.(ff).data.u, s2, u);
        end
    end
end
