function om = build_cost_params(om, force)
%BUILD_COST_PARAMS  Builds and saves the full generalized cost parameters.
%   OM.BUILD_COST_PARAMS()
%   OM.BUILD_COST_PARAMS('force')
%
%   Builds the full set of cost parameters from the individual named
%   sub-sets added via ADD_COSTS. Skips the building process if it has
%   already been done, unless a second input argument is present.
%
%   These cost parameters can be retrieved by calling GET_COST_PARAMS
%   and the user-defined costs evaluated by calling COMPUTE_COST.
%
%   See also OPT_MODEL, ADD_COSTS, GET_COST_PARAMS, COMPUTE_COST.

%   MATPOWER
%   Copyright (c) 2008-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin > 1 || isempty(om.cost.params)
    %% initialize parameters
    nw = om.cost.N;
    nnzN = 0;
    nnzH = 0;

    %% calls to substruct() are relatively expensive, so we pre-build the
    %% structs for addressing cell and numeric array fields, updating only
    %% the subscripts before use
    sc = struct('type', {'.', '{}'}, 'subs', {'', 1});  %% cell array field
    sn = sc; sn(2).type = '()';                         %% num array field

    for k = 1:om.cost.NS
        name = om.cost.order(k).name;
        idx  = om.cost.order(k).idx;
        if isempty(idx)
            nnzN = nnzN + nnz(om.cost.data.N.(name));
            if isfield(om.cost.data.H, name)
                nnzH = nnzH + nnz(om.cost.data.H.(name));
            end
        else
            % (calls to substruct() are relatively expensive ...
            % sc = substruct('.', name, '{}', idx);
            % ... so replace it with these more efficient lines)
            sc(1).subs = name;
            sc(2).subs = idx;
            nnzN = nnzN + nnz(subsref(om.cost.data.N, sc));
            if isfield(om.cost.data.H, name)
                nnzH = nnzH + nnz(subsref(om.cost.data.H, sc));
            end
        end
    end
    NNt = sparse([], [], [], om.var.N, nw, nnzN);   %% use NN transpose for speed
    Cw = zeros(nw, 1);
    H = sparse([], [], [], nw, nw, nnzH);   %% default => no quadratic term
    dd = ones(nw, 1);                       %% default => linear
    rh = Cw;                                %% default => no shift
    kk = Cw;                                %% default => no dead zone
    mm = dd;                                %% default => no scaling
    
    %% fill in each piece
    for k = 1:om.cost.NS
        name = om.cost.order(k).name;
        idx  = om.cost.order(k).idx;
        if isempty(idx)
            N = om.cost.idx.N.(name);       %% number of rows to add
        else
            % (calls to substruct() are relatively expensive ...
            % sn = substruct('.', name, '()', idx);
            % sc = substruct('.', name, '{}', idx);
            % ... so replace them with these more efficient lines)
            sn(1).subs = name;
            sn(2).subs = idx;
            sc(1).subs = name;
            sc(2).subs = idx;
            N = subsref(om.cost.idx.N, sn); %% number of rows to add
        end
        if N                                %% non-zero number of rows to add
            if isempty(idx)
                Nk = om.cost.data.N.(name);         %% N for kth cost set
                i1 = om.cost.idx.i1.(name);         %% starting row index
                iN = om.cost.idx.iN.(name);         %% ending row index
                vs = om.cost.data.vs.(name);        %% var sets
            else
                Nk = subsref(om.cost.data.N, sc);   %% N for kth cost set
                i1 = subsref(om.cost.idx.i1, sn);   %% starting row index
                iN = subsref(om.cost.idx.iN, sn);   %% ending row index
                vs = subsref(om.cost.data.vs, sc);  %% var sets
            end
            if isempty(vs)          %% full rows
                if size(Nk,2) == om.var.N
                    NNt(:, i1:iN) = Nk';     %% assign as columns in transpose for speed
                else                %% must have added vars since adding
                                    %% this cost set
                    NNt(1:size(Nk,2), i1:iN) = Nk';  %% assign as columns in transpose for speed
                end
            else                    %% selected columns
                jj = om.varsets_idx(vs);    %% column indices for var set
                Ni = sparse(N, om.var.N);
                Ni(:, jj) = Nk;
                NNt(:, i1:iN) = Ni';    %% assign as columns in transpose for speed
            end

            if isempty(idx)
                Cw(i1:iN) = om.cost.data.Cw.(name);
                if isfield(om.cost.data.H, name)
                    H(i1:iN, i1:iN) = om.cost.data.H.(name);
                end
                if isfield(om.cost.data.dd, name)
                    dd(i1:iN) = om.cost.data.dd.(name);
                end
                if isfield(om.cost.data.rh, name)
                    rh(i1:iN) = om.cost.data.rh.(name);
                end
                if isfield(om.cost.data.kk, name)
                    kk(i1:iN) = om.cost.data.kk.(name);
                end
                if isfield(om.cost.data.mm, name)
                    mm(i1:iN) = om.cost.data.mm.(name);
                end
            else
                Cw(i1:iN) = subsref(om.cost.data.Cw, sc);
                if isfield(om.cost.data.H, name) && ~isempty(subsref(om.cost.data.H, sc))
                    H(i1:iN, i1:iN) = subsref(om.cost.data.H, sc);
                end
                if isfield(om.cost.data.dd, name) && ~isempty(subsref(om.cost.data.dd, sc))
                    dd(i1:iN) = subsref(om.cost.data.dd, sc);
                end
                if isfield(om.cost.data.rh, name) && ~isempty(subsref(om.cost.data.rh, sc))
                    rh(i1:iN) = subsref(om.cost.data.rh, sc);
                end
                if isfield(om.cost.data.kk, name) && ~isempty(subsref(om.cost.data.kk, sc))
                    kk(i1:iN) = subsref(om.cost.data.kk, sc);
                end
                if isfield(om.cost.data.mm, name) && ~isempty(subsref(om.cost.data.mm, sc))
                    mm(i1:iN) = subsref(om.cost.data.mm, sc);
                end
            end
        end
    end

    %% save in object   
    om.cost.params = struct( ...
        'N', NNt', 'Cw', Cw, 'H', H, 'dd', dd, 'rh', rh, 'kk', kk, 'mm', mm );
end
