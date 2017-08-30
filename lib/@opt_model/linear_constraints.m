function [A, l, u] = linear_constraints(om)
%LINEAR_CONSTRAINTS  Builds and returns the full set of linear constraints.
%   [A, L, U] = OM.LINEAR_CONSTRAINTS()
%   Builds the full set of linear constraints based on those added by
%   ADD_LIN_CONSTRAINTS.
%
%       L <= A * x <= U
%
%   Example:
%       [A, l, u] = om.linear_constraints();
%
%   See also OPT_MODEL, ADD_LIN_CONSTRAINTS.

%   MATPOWER
%   Copyright (c) 2008-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.


%% initialize A, l and u
nnzA = 0;

%% calls to substruct() are relatively expensive, so we pre-build the
%% structs for addressing cell and numeric array fields, updating only
%% the subscripts before use
sc = struct('type', {'.', '{}'}, 'subs', {'', 1});  %% cell array field
sn = sc; sn(2).type = '()';                         %% num array field

for k = 1:om.lin.NS
    name = om.lin.order(k).name;
    idx  = om.lin.order(k).idx;
    if isempty(idx)
        nnzA = nnzA + nnz(om.lin.data.A.(name));
    else
        % (calls to substruct() are relatively expensive ...
        % sc = substruct('.', name, '{}', idx);
        % ... so replace it with these more efficient lines)
        sc(1).subs = name;
        sc(2).subs = idx;
        nnzA = nnzA + nnz(subsref(om.lin.data.A, sc));
    end
end
At = sparse([], [], [], om.var.N, om.lin.N, nnzA);  %% use A transpose for speed
u = Inf(om.lin.N, 1);
l = -u;

%% fill in each piece
for k = 1:om.lin.NS
    name = om.lin.order(k).name;
    idx  = om.lin.order(k).idx;
    if isempty(idx)
        N = om.lin.idx.N.(name);
    else
        % (calls to substruct() are relatively expensive ...
        % sn = substruct('.', name, '()', idx);
        % sc = substruct('.', name, '{}', idx);
        % ... so replace them with these more efficient lines)
        sn(1).subs = name;
        sn(2).subs = idx;
        sc(1).subs = name;
        sc(2).subs = idx;
        N = subsref(om.lin.idx.N, sn);
    end
    if N                                %% non-zero number of rows to add
        if isempty(idx)
            Ak = om.lin.data.A.(name);          %% A for kth linear constrain set
            i1 = om.lin.idx.i1.(name);          %% starting row index
            iN = om.lin.idx.iN.(name);          %% ending row index
            vs = om.lin.data.vs.(name);         %% var sets
        else
            Ak = subsref(om.lin.data.A, sc);    %% A for kth linear constrain set
            i1 = subsref(om.lin.idx.i1, sn);    %% starting row index
            iN = subsref(om.lin.idx.iN, sn);    %% ending row index
            vs = subsref(om.lin.data.vs, sc);   %% var sets
        end
        if isempty(vs)          %% full rows
            if size(Ak,2) == om.var.N
                At(:, i1:iN) = Ak';     %% assign as columns in transpose for speed
            else                %% must have added vars since adding
                                %% this constraint set
                At(1:size(Ak,2), i1:iN) = Ak';  %% assign as columns in transpose for speed
            end
        else                    %% selected columns
            jj = om.varsets_idx(vs);    %% column indices for var set
            Ai = sparse(N, om.var.N);
            Ai(:, jj) = Ak;
            At(:, i1:iN) = Ai';     %% assign as columns in transpose for speed
        end

        if isempty(idx)
            l(i1:iN) = om.lin.data.l.(name);
            u(i1:iN) = om.lin.data.u.(name);
        else
            l(i1:iN) = subsref(om.lin.data.l, sc);
            u(i1:iN) = subsref(om.lin.data.u, sc);
        end
    end
end
A = At';
