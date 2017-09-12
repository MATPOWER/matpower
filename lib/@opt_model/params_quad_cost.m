function [Q, c, k, vs] = params_quad_cost(om, name, idx)
%PARAMS_QUAD_COST  Returns the cost parameters for quadratic costs.
%   [Q, C] = OM.PARAMS_QUAD_COST()
%   [Q, C] = OM.PARAMS_QUAD_COST(NAME)
%   [Q, C] = OM.PARAMS_QUAD_COST(NAME, IDX)
%   [Q, C, K] = OM.PARAMS_QUAD_COST(...)
%   [Q, C, K, VS] = OM.PARAMS_QUAD_COST(...)
%
%   With no input parameters, it assembles and returns the parameters
%   for the aggregate quadratic cost from all quadratic cost sets added
%   using ADD_QUAD_COST. The values of these parameters are cached
%   for subsequent calls. The parameters are Q, C, and optionally K,
%   where the quadratic cost is of the form
%       F(X) = 1/2 * X'*Q*X + C'*X + K
%
%   If a NAME is provided then it simply returns the parameters for the
%   corresponding named set. Likewise for indexed named sets specified
%   by NAME and IDX. In this case, Q and K may be vectors, corresponding
%   to a cost function of the form
%       F(X) = 1/2 * Q .* X.^2 + C .* X + K
%
%   An optional 4th output argument VS indicates the variable sets used by
%   this cost set. The size of Q and C will be consistent with VS.
%
%   See also OPT_MODEL, ADD_QUAD_COST.

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
        if prod(size(om.qdc.idx.i1.(name))) == 1    %% simple named set
            Q = om.qdc.data.Q.(name);
            c = om.qdc.data.c.(name);
            k = om.qdc.data.k.(name);
            if nargout > 3
                vs = om.qdc.data.vs.(name);
            end
        else                                        %% indexing required
            error('@opt_model/params_quad_cost: quadratic cost set ''%s'' requires an IDX arg', name);
        end
    else                            %% indexed named set
        % (calls to substruct() are relatively expensive ...
        % s = substruct('.', name, '{}', idx);
        % ... so replace it with these more efficient lines)
        sc = struct('type', {'.', '{}'}, 'subs', {name, idx});
        Q = subsref(om.qdc.data.Q, sc);
        c = subsref(om.qdc.data.c, sc);
        k = subsref(om.qdc.data.k, sc);
        if nargout > 3
            vs = subsref(om.qdc.data.vs, sc);
        end
    end
else                %% aggregate
    cache = om.qdc.params;
    if isempty(cache)       %% build the aggregate
        nx = om.var.N;          %% number of variables
        Qt = sparse(nx, nx);    %% transpose of quadratic coefficients
        c = zeros(nx, 1);       %% linear coefficients
        k = 0;                  %% constant term
        for i = 1:om.qdc.NS
            name = om.qdc.order(i).name;
            idx  = om.qdc.order(i).idx;
            [Qk, ck, kk, vs] = om.params_quad_cost(name, idx);
            haveQ = ~isempty(Qk);
            havec = ~isempty(ck);
            nk = max(size(Qk, 1), size(ck, 1));     %% size of Qk and/or ck
            if isempty(vs)
                if nk == nx     %% full size
                    if size(Qk, 2) == 1     %% Qk is a column vector
                        Qkt_full = spdiags(Qk, 0, nx, nx);
                    elseif haveQ            %% Qk is a matrix
                        Qkt_full = Qk';
                    end
                    if havec
                        ck_full = ck;
                    end
                else            %% vars added since adding this cost set
                    if size(Qk, 2) == 1     %% Qk is a column vector
                        Qkt_full = sparse(1:nk, 1:nk, Qk, nx, nx);
                    elseif haveQ            %% Qk is a matrix
                        Qk_all_cols = sparse(nk, nx);
                        Qk_all_cols(:, 1:nk) = Qk;
                        Qkt_full(:, 1:nk) = Qk_all_cols';
                    end
                    if havec
                        ck_full = zeros(nx, 1);
                        ck_full(1:nk) = ck;
                    end
                end
            else
                jj = om.varsets_idx(vs);    %% indices for var set
                if size(Qk, 2) == 1     %% Qk is a column vector
                    Qkt_full = sparse(jj, jj, Qk, nx, nx);
                elseif haveQ            %% Qk is a matrix
                    Qk_all_cols = sparse(nk, nx);
                    Qk_all_cols(:, jj) = Qk;
                    Qkt_full = sparse(nx, nx);
                    Qkt_full(:, jj) = Qk_all_cols';
                end
                if havec
                    ck_full = zeros(nx, 1);
                    ck_full(jj) = ck;
                end
            end
            if haveQ
                Qt = Qt + Qkt_full;
            end
            if havec
                c = c + ck_full;
            end
            k = k + sum(kk);
        end
        Q = Qt';

        %% cache aggregated parameters
        om.qdc.params = struct('Q', Q, 'c', c, 'k', k);
    else                    %% return cached values
        Q = cache.Q;
        c = cache.c;
        k = cache.k;
    end
    if nargout > 3
        vs = {};
    end
end
