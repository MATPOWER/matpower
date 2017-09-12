function [cp, vs, i1, iN] = params_legacy_cost(om, name, idx)
%PARAMS_LEGACY_COST  Returns cost parameters for legacy user-defined costs.
%   CP = OM.PARAMS_LEGACY_COST()
%   CP = OM.PARAMS_LEGACY_COST(NAME)
%   CP = OM.PARAMS_LEGACY_COST(NAME, IDX)
%   [CP, VS] = OM.PARAMS_LEGACY_COST(...)
%   [CP, VS, I1, IN] = OM.PARAMS_LEGACY_COST(...)
%
%   With no input parameters, it assembles and returns the parameters
%   for the aggregate legacy user-defined cost from all legacy cost sets
%   added using ADD_LEGACY_COST. The values of these parameters are cached
%   for subsequent calls. The parameters are contained in the struct CP,
%   described below.
%
%   If a NAME is provided then it simply returns parameter struct CP for the
%   corresponding named set. Likewise for indexed named sets specified
%   by NAME and IDX.
%
%   An optional 2nd output argument VS indicates the variable sets used by
%   this cost set. The size of CP.N will be consistent with VS.
%
%   If NAME is provided, optional 3rd and 4th output arguments I1 and IN
%   indicate the starting and ending row indices of the corresponding
%   cost set in the full aggregate cost matrix CP.N.
%
%   Let X refer to the vector formed by combining the corresponding varsets VS,
%   and F_U(X, CP) be the cost at X corresponding to the cost parameters
%   contained in CP, where CP is a struct with the following fields:
%       N      - nw x nx sparse matrix (optional, identity matrix by default)
%       Cw     - nw x 1 vector
%       H      - nw x nw sparse matrix (optional, all zeros by default)
%       dd, mm - nw x 1 vectors (optional, all ones by default)
%       rh, kk - nw x 1 vectors (optional, all zeros by default)
%
%   These parameters are used as follows to compute F_U(X, CP)
%
%       R  = N*X - rh
%
%               /  kk(i),  R(i) < -kk(i)
%       K(i) = <   0,     -kk(i) <= R(i) <= kk(i)
%               \ -kk(i),  R(i) > kk(i)
%
%       RR = R + K
%
%       U(i) =  /  0, -kk(i) <= R(i) <= kk(i)
%               \  1, otherwise
%
%       DDL(i) = /  1, dd(i) = 1
%                \  0, otherwise
%
%       DDQ(i) = /  1, dd(i) = 2
%                \  0, otherwise
%
%       Dl = diag(mm) * diag(U) * diag(DDL)
%       Dq = diag(mm) * diag(U) * diag(DDQ)
%
%       w = (Dl + Dq * diag(RR)) * RR
%
%       F_U(X, CP) = 1/2 * w'*H*w + Cw'*w
%
%   See also OPT_MODEL, ADD_LEGACY_COST, EVAL_LEGACY_COST.

%   MATPOWER
%   Copyright (c) 2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin > 1       %% individual set
    if nargin < 3
        idx = {};
    end
    if isempty(idx)                 %% name, no index provided
        if prod(size(om.cost.idx.i1.(name))) == 1   %% simple named set
            N  = om.cost.data.N.( name);
            Cw = om.cost.data.Cw.(name);
            [nw, nx] = size(N);
            e1 = ones(nw, 1);
            e0 = zeros(nw, 1);
            cp = struct(...
                'N',    N, ...
                'Cw',   Cw, ...
                'H',    sparse(nw, nw), ... %% default => no quadratic term
                'dd',   e1, ...             %% default => linear
                'rh',   e0, ...             %% default => no shift
                'kk',   e0, ...             %% default => no dead zone
                'mm',   e1 );               %% default => no scaling
            if isfield(om.cost.data.H, name) && ~isempty(om.cost.data.H.(name))
                cp.H = om.cost.data.H.(name);
            end
            if isfield(om.cost.data.dd, name) && ~isempty(om.cost.data.dd.(name))
                cp.dd = om.cost.data.dd.(name);
            end
            if isfield(om.cost.data.rh, name) && ~isempty(om.cost.data.rh.(name))
                cp.rh = om.cost.data.rh.(name);
            end
            if isfield(om.cost.data.kk, name) && ~isempty(om.cost.data.kk.(name))
                cp.kk = om.cost.data.kk.(name);
            end
            if isfield(om.cost.data.mm, name) && ~isempty(om.cost.data.mm.(name))
                cp.mm = om.cost.data.mm.(name);
            end
            if nargout > 1
                vs = om.cost.data.vs.(name);
                if nargout > 3
                    i1 = om.cost.idx.i1.(name);     %% starting row index
                    iN = om.cost.idx.iN.(name);     %% ending row index
                end
            end
        else                                        %% indexing required
            error('@opt_model/params_legacy_cost: legacy cost set ''%s'' requires an IDX arg', name);
        end
    else                            %% indexed named set
        % (calls to substruct() are relatively expensive ...
        % s = substruct('.', name, '{}', idx);
        % ... so replace it with these more efficient lines)
        sc = struct('type', {'.', '{}'}, 'subs', {name, idx});
        N  = subsref(om.cost.data.N,  sc);
        Cw = subsref(om.cost.data.Cw, sc);
        [nw, nx] = size(N);
        e1 = ones(nw, 1);
        e0 = zeros(nw, 1);
        cp = struct(...
            'N',    N, ...
            'Cw',   Cw, ...
            'H',    sparse(nw, nw), ... %% default => no quadratic term
            'dd',   e1, ...             %% default => linear
            'rh',   e0, ...             %% default => no shift
            'kk',   e0, ...             %% default => no dead zone
            'mm',   e1 );               %% default => no scaling
        if isfield(om.cost.data.H, name) && ~isempty(subsref(om.cost.data.H, sc))
            cp.H = subsref(om.cost.data.H,  sc);
        end
        if isfield(om.cost.data.dd, name) && ~isempty(subsref(om.cost.data.dd, sc))
            cp.dd = subsref(om.cost.data.dd, sc);
        end
        if isfield(om.cost.data.rh, name) && ~isempty(subsref(om.cost.data.rh, sc))
            cp.rh = subsref(om.cost.data.rh, sc);
        end
        if isfield(om.cost.data.kk, name) && ~isempty(subsref(om.cost.data.kk, sc))
            cp.kk = subsref(om.cost.data.kk, sc);
        end
        if isfield(om.cost.data.mm, name) && ~isempty(subsref(om.cost.data.mm, sc))
            cp.mm = subsref(om.cost.data.mm, sc);
        end
        if nargout > 1
            vs = subsref(om.cost.data.vs, sc);
            if nargout > 3
                sn = sc; sn(2).type = '()';         %% num array field
                i1 = subsref(om.cost.idx.i1, sn);   %% starting row index
                iN = subsref(om.cost.idx.iN, sn);   %% ending row index
            end
        end
    end
else                %% aggregate
    cache = om.cost.params;
    if isempty(cache)       %% build the aggregate
        nx = om.var.N;          %% number of variables
        nw = om.cost.N;         %% number of cost rows
        Nt = sparse(nx, nw);    %% use N transpose for speed
        Cw = zeros(nw, 1);
        H = sparse(nw, nw);     %% default => no quadratic term
        dd = ones(nw, 1);       %% default => linear
        rh = Cw;                %% default => no shift
        kk = Cw;                %% default => no dead zone
        mm = dd;                %% default => no scaling

        %% fill in each piece
        for k = 1:om.cost.NS
            name = om.cost.order(k).name;
            idx  = om.cost.order(k).idx;
            [cp, vs, i1, iN] = om.params_legacy_cost(name, idx);
            [mk, nk] = size(cp.N);      %% size of Nk
            if mk
                Nkt_full = sparse(nx, nw);
                if isempty(vs)
                    if nk == nx     %% full size
                        Nkt_full(:, i1:iN) = cp.N';
                    else            %% vars added since adding this cost set
                        Nk_all_cols = sparse(mk, nx);
                        Nk_all_cols(:, 1:nk) = cp.N;
                        Nkt_full(:, i1:iN) = Nk_all_cols';
                    end
                else
                    jj = om.varsets_idx(vs);    %% indices for var set
                    Nk_all_cols = sparse(mk, nx);
                    Nk_all_cols(:, jj) = cp.N;
                    Nkt_full(:, i1:iN) = Nk_all_cols';
                end
                Nt = Nt + Nkt_full;
                Cw(i1:iN) = cp.Cw;
                H_all_cols = sparse(mk, nw);
                H_all_cols(:, i1:iN) = cp.H';
                H(:, i1:iN) = H_all_cols';
                dd(i1:iN) = cp.dd;
                rh(i1:iN) = cp.rh;
                kk(i1:iN) = cp.kk;
                mm(i1:iN) = cp.mm;
            end
        end

        %% cache aggregated parameters
        om.cost.params = struct( ...
            'N', Nt', 'Cw', Cw, 'H', H, 'dd', dd, 'rh', rh, 'kk', kk, 'mm', mm );
        cp = om.cost.params;
    else                    %% return cached values
        cp = cache;
    end
    if nargout > 1
        vs = {};
    end
end
