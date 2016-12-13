function om = build_cost_params(om, force)
%BUILD_COST_PARAMS  Builds and saves the full generalized cost parameters.
%   OM = BUILD_COST_PARAMS(OM)
%   OM = BUILD_COST_PARAMS(OM, 'force')
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

if nargin > 1 || ~isfield(om.cost.params, 'N')
    %% initialize parameters
    nw = om.cost.N;
    nnzN = 0;
    nnzH = 0;
    s = struct('type', {'.', '{}'}, 'subs', {'', 1});
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
            % s = substruct('.', name, '{}', idx);
            % ... so replace it with these more efficient lines)
            s(1).subs = name;
            s(2).subs = idx;
            nnzN = nnzN + nnz(subsref(om.cost.data.N, s));
            if isfield(om.cost.data.H, name)
                nnzH = nnzH + nnz(subsref(om.cost.data.H, s));
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
    s2 = s;
    s(2).type = '()';
    s1 = s;
    for k = 1:om.cost.NS
        name = om.cost.order(k).name;
        idx  = om.cost.order(k).idx;
        if isempty(idx)
            N = om.cost.idx.N.(name);       %% number of rows to add
        else
            % (calls to substruct() are relatively expensive ...
            % s1 = substruct('.', name, '()', idx);
            % s2 = substruct('.', name, '{}', idx);
            % ... so replace them with these more efficient lines)
            s1(1).subs = name;
            s1(2).subs = idx;
            s2(1).subs = name;
            s2(2).subs = idx;
            N = subsref(om.cost.idx.N, s1); %% number of rows to add
        end
        if N                                %% non-zero number of rows to add
            if isempty(idx)
                Nk = om.cost.data.N.(name);         %% N for kth cost set
                i1 = om.cost.idx.i1.(name);         %% starting row index
                iN = om.cost.idx.iN.(name);         %% ending row index
                vsl = om.cost.data.vs.(name);       %% var set list
            else
                Nk = subsref(om.cost.data.N, s2);   %% N for kth cost set
                i1 = subsref(om.cost.idx.i1, s1);   %% starting row index
                iN = subsref(om.cost.idx.iN, s1);   %% ending row index
                vsl = subsref(om.cost.data.vs, s2); %% var set list
            end
            if isempty(vsl)         %% full rows
                if size(Nk,2) == om.var.N
                    NNt(:, i1:iN) = Nk';     %% assign as columns in transpose for speed
                else                %% must have added vars since adding
                                    %% this cost set
                    NNt(1:size(Nk,2), i1:iN) = Nk';  %% assign as columns in transpose for speed
                end
            else                    %% selected columns
                kN = 0;                             %% initialize last col of Nk used
                Ni = sparse(N, om.var.N);
                for v = 1:length(vsl)
                    % (calls to substruct() are relatively expensive ...
                    % s = substruct('.', vsl(v).name, '()', vsl(v).idx);
                    % ... so replace it with these more efficient lines)
                    s(1).subs = vsl(v).name;
                    s(2).subs = vsl(v).idx;
                    j1 = subsref(om.var.idx.i1, s); %% starting column in N
                    jN = subsref(om.var.idx.iN, s); %% ending column in N
                    k1 = kN + 1;                    %% starting column in Nk
                    kN = kN + subsref(om.var.idx.N, s);%% ending column in Nk
                    Ni(:, j1:jN) = Nk(:, k1:kN);
                end
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
                Cw(i1:iN) = subsref(om.cost.data.Cw, s2);
                if isfield(om.cost.data.H, name) && ~isempty(subsref(om.cost.data.H, s2))
                    H(i1:iN, i1:iN) = subsref(om.cost.data.H, s2);
                end
                if isfield(om.cost.data.dd, name) && ~isempty(subsref(om.cost.data.dd, s2))
                    dd(i1:iN) = subsref(om.cost.data.dd, s2);
                end
                if isfield(om.cost.data.rh, name) && ~isempty(subsref(om.cost.data.rh, s2))
                    rh(i1:iN) = subsref(om.cost.data.rh, s2);
                end
                if isfield(om.cost.data.kk, name) && ~isempty(subsref(om.cost.data.kk, s2))
                    kk(i1:iN) = subsref(om.cost.data.kk, s2);
                end
                if isfield(om.cost.data.mm, name) && ~isempty(subsref(om.cost.data.mm, s2))
                    mm(i1:iN) = subsref(om.cost.data.mm, s2);
                end
            end
        end
    end

    %% save in object   
    om.cost.params = struct( ...
        'N', NNt', 'Cw', Cw, 'H', H, 'dd', dd, 'rh', rh, 'kk', kk, 'mm', mm );
end
