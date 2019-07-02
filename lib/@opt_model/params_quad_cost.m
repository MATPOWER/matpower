function [Q, c, K, vs] = params_quad_cost(om, name, idx)
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
%   See https://matpower.org for more info.

if nargin > 1       %% individual set
    if nargin < 3
        idx = {};
    end
    if isempty(idx)                 %% name, no index provided
        if numel(om.qdc.idx.i1.(name)) == 1     %% simple named set
            Q = om.qdc.data.Q.(name);
            c = om.qdc.data.c.(name);
            K = om.qdc.data.k.(name);
            if nargout > 3
                vs = om.qdc.data.vs.(name);
            end
        else                                    %% indexing required
            error('@opt_model/params_quad_cost: quadratic cost set ''%s'' requires an IDX arg', name);
        end
    else                            %% indexed named set
        % (calls to substruct() are relatively expensive ...
        % s = substruct('.', name, '{}', idx);
        % ... so replace it with these more efficient lines)
        sc = struct('type', {'.', '{}'}, 'subs', {name, idx});
        Q = subsref(om.qdc.data.Q, sc);
        c = subsref(om.qdc.data.c, sc);
        K = subsref(om.qdc.data.k, sc);
        if nargout > 3
            vs = subsref(om.qdc.data.vs, sc);
        end
    end
else                %% aggregate
    cache = om.qdc.params;
    if isempty(cache)       %% build the aggregate
        nx = om.var.N;          %% number of variables
        Q_ijv = cell(om.qdc.NS, 3); %% indices/values to construct Q
        c_ijv = cell(om.qdc.NS, 3); %% indices/values to construct c
        K = 0;                  %% constant term
        for k = 1:om.qdc.NS
            name = om.qdc.order(k).name;
            idx  = om.qdc.order(k).idx;
            [Qk, ck, kk, vs] = om.params_quad_cost(name, idx);
            haveQ = ~isempty(Qk);
            havec = ~isempty(ck);
            if haveQ
                [i, j, v] = find(Qk);
            end
            if havec
                [ic, jc, vc] = find(ck);
            end
            if isempty(vs)
                if size(Qk, 2) == 1     %% Qk is a column vector
                    Q_ijv(k, :) = {i, i, v};
                elseif haveQ            %% Qk is a matrix
                    Q_ijv(k, :) = {i, j, v};
                end
                if havec
                    c_ijv(k, :) = {ic, jc, vc};
                end
            else
                jj = om.varsets_idx(vs)';    %% indices for var set
                if size(Qk, 2) == 1     %% Qk is a column vector
                    Q_ijv(k, :) = {jj(i), jj(i), v};
                elseif haveQ            %% Qk is a matrix
                    Q_ijv(k, :) = {jj(i), jj(j), v};
                end
                if havec
                    c_ijv(k, :) = {jj(ic), jc, vc};
                end
            end
            K = K + sum(kk);
        end
        Q = sparse( vertcat(Q_ijv{:,1}), ...
                    vertcat(Q_ijv{:,2}), ...
                    vertcat(Q_ijv{:,3}), nx, nx);
        c = accumarray(vertcat(c_ijv{:,1}), vertcat(c_ijv{:,3}), [nx 1]);

        %% cache aggregated parameters
        om.qdc.params = struct('Q', Q, 'c', c, 'k', K);
    else                    %% return cached values
        Q = cache.Q;
        c = cache.c;
        K = cache.k;
    end
    if nargout > 3
        vs = {};
    end
end
