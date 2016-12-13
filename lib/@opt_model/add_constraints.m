function om = add_constraints(om, name, idx, varargin)
%ADD_CONSTRAINTS  Adds a set of constraints to the model.
%   OM = ADD_CONSTRAINTS(OM, NAME, A, L, U);
%   OM = ADD_CONSTRAINTS(OM, NAME, A, L, U, VARSETS);
%   OM = ADD_CONSTRAINTS(OM, NAME, DIM_LIST);
%   OM = ADD_CONSTRAINTS(OM, NAME, IDX_LIST, A, L, U);
%   OM = ADD_CONSTRAINTS(OM, NAME, IDX_LIST, A, L, U, VARSETS);
%   OM = ADD_CONSTRAINTS(OM, NAME, N, 'NON-LINEAR');
%
%   Linear constraints are of the form L <= A * x <= U, where
%   x is a vector made of of the vars specified in VARSETS (in
%   the order given). This allows the A matrix to be defined only
%   in terms of the relevant variables without the need to manually
%   create a lot of zero columns. If VARSETS is empty, x is taken
%   to be the full vector of all optimization variables. If L or 
%   U are empty, they are assumed to be appropriately sized vectors
%   of -Inf and Inf, respectively.
%
%   For nonlinear constraints, the 3rd argument, N, is the number
%   of constraints in the set. Currently, this is used internally
%   by MATPOWER, but there is no way for the user to specify
%   additional nonlinear constraints.
%
%   Examples:
%       om = add_constraints(om, 'vl', Avl, lvl, uvl, {'Pg', 'Qg'});
%       om = add_constraints(om, 'Pmis', nb, 'nonlinear');
%
%       om = add_constraints(om, 'R', {2, 3});
%       for i = 1:2
%         for j = 1:3
%           om = add_constraints(om, 'R', {i, j}, A{i,j}, ...);
%         end
%       end
%
%   See also OPT_MODEL, LINEAR_CONSTRAINTS.

%   MATPOWER
%   Copyright (c) 2008-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

nonlin = 0;
if iscell(idx)
    if ~isempty(varargin)       %% linear: indexed named set
        % (calls to substruct() are relatively expensive ...
        % s1 = substruct('.', name, '()', idx);
        % s2 = substruct('.', name, '{}', idx);
        % ... so replace them with these more efficient lines)
        s1 = struct('type', {'.', '()'}, 'subs', {name, idx});
        s2 = s1;
        s2(2).type = '{}';
        
        %% prevent duplicate named constraint sets
        if subsref(om.lin.idx.i1, s1) ~= 0
            str = '%d'; for m = 2:length(idx), str = [str ',%d']; end
            nname = sprintf(['%s(' str, ')'], name, idx{:});
            error('@opt_model/add_constraints: linear constraint set named ''%s'' already exists', nname);
        end
        
        A = varargin{1};
        args = varargin(2:end);
    else                        %% linear: just setting dimensions for indexed set
        %% prevent duplicate named constraint sets
        if isfield(om.lin.idx.N, name)
            error('@opt_model/add_constraints: linear constraint set named ''%s'' already exists', name);
        end
        
        A = sparse(0,0);
        args = {};
    end
else
    if length(varargin) == 1    %% non-linear
        %% prevent duplicate named constraint sets
        if isfield(om.nln.idx.N, name)
            error('@opt_model/add_constraints: nonlinear constraint set named ''%s'' already exists', name);
        end
        
        nonlin = 1;
        N = idx;
        args = {};
    else                        %% linear: simple named set
        %% prevent duplicate named constraint sets
        if isfield(om.lin.idx.N, name)
            error('@opt_model/add_constraints: linear constraint set named ''%s'' already exists', name);
        end
        
        A = idx;
        args = varargin;
    end
    idx = {};
end
nargs = length(args);

if nonlin           %% nonlinear
    %% add info about this nonlinear constraint set
    om.nln.idx.i1.(name) = om.nln.N + 1;    %% starting index
    om.nln.idx.iN.(name) = om.nln.N + N;    %% ending index
    om.nln.idx.N.(name)  = N;               %% number of constraints
    
    %% update number of nonlinear constraints and constraint sets
    om.nln.N  = om.nln.idx.iN.(name);
    om.nln.NS = om.nln.NS + 1;
    
    %% add to ordered list of nonlinear constraint sets
    om.nln.order(om.nln.NS).name = name;
    om.nln.order(om.nln.NS).idx  = {};
elseif nargs == 0   %% linear: just setting dimensions for indexed set
    %% use column vector if single dimension
    if length(idx) == 1
        idx = {idx{:}, 1};
    end
    
    %% add info about this linear constraint set
    om.lin.idx.i1.(name)  = zeros(idx{:});  %% starting index
    om.lin.idx.iN.(name)  = zeros(idx{:});  %% ending index
    om.lin.idx.N.(name)   = zeros(idx{:});  %% number of constraints
    om.lin.data.A.(name)  = cell(idx{:});
    om.lin.data.l.(name)  = cell(idx{:});
    om.lin.data.u.(name)  = cell(idx{:});
    om.lin.data.vs.(name) = cell(idx{:});
else                %% linear
    if nargs >= 3
        [l, u, varsets] = deal(args{1:3});
    else
        varsets = {};
        if nargs >= 2
            [l, u] = deal(args{1:2});
        else
            u = [];
            if nargs >= 1
                l = args{1};
            else
                l = [];
            end
        end
    end
    
    [N, M] = size(A);
    if isempty(l)                   %% default l is -Inf
        l = -Inf(N, 1);
    end
    if isempty(u)                   %% default u is Inf
        u = Inf(N, 1);
    end
    if ~isempty(varsets) && iscell(varsets)
        empty_cells = cell(1, length(varsets));
        [empty_cells{:}] = deal({});    %% empty cell arrays
        varsets = struct('name', varsets, 'idx', empty_cells);
    end
    
    %% check sizes
    if size(l, 1) ~= N || size(u, 1) ~= N
        error('@opt_model/add_constraints: sizes of A, l and u must match');
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
    if M ~= nv
        error('@opt_model/add_constraints: number of columns of A does not match\nnumber of variables, A is %d x %d, nv = %d\n', N, M, nv);
    end
    if isempty(idx)     %% linear: simple named set
        %% add info about this linear constraint set
        om.lin.idx.i1.(name)  = om.lin.N + 1;   %% starting index
        om.lin.idx.iN.(name)  = om.lin.N + N;   %% ending index
        om.lin.idx.N.(name)   = N;              %% number of constraints
        om.lin.data.A.(name)  = A;
        om.lin.data.l.(name)  = l;
        om.lin.data.u.(name)  = u;
        om.lin.data.vs.(name) = varsets;
        
        %% update number of linear constraints and constraint sets
        om.lin.N  = om.lin.idx.iN.(name);
        om.lin.NS = om.lin.NS + 1;
        
        %% add to ordered list of linear constraint sets
        om.lin.order(om.lin.NS).name = name;
        om.lin.order(om.lin.NS).idx  = {};
    else                %% linear: indexed named set
        %% add info about this linear constraint set
        om.lin.idx.i1  = subsasgn(om.lin.idx.i1, s1, om.lin.N + 1); %% starting index
        om.lin.idx.iN  = subsasgn(om.lin.idx.iN, s1, om.lin.N + N); %% ending index
        om.lin.idx.N   = subsasgn(om.lin.idx.N,  s1, N);            %% number of constraints
        om.lin.data.A  = subsasgn(om.lin.data.A, s2, A);
        om.lin.data.l  = subsasgn(om.lin.data.l, s2, l);
        om.lin.data.u  = subsasgn(om.lin.data.u, s2, u);
        om.lin.data.vs = subsasgn(om.lin.data.vs, s2, varsets);
        
        %% update number of linear constraints and constraint sets
        om.lin.N  = subsref(om.lin.idx.iN, s1);
        om.lin.NS = om.lin.NS + 1;
        
        %% add to ordered list of linear constraint sets
        om.lin.order(om.lin.NS).name = name;
        om.lin.order(om.lin.NS).idx  = idx;
    end
end