function [A, l, u, vs, i1, iN] = params_lin_constraint(om, name, idx)
%PARAMS_LIN_CONSTRAINT  Builds and returns linear constraint parameters.
%   [A, L, U] = OM.PARAMS_LIN_CONSTRAINT()
%   [A, L, U] = OM.PARAMS_LIN_CONSTRAINT(NAME)
%   [A, L, U] = OM.PARAMS_LIN_CONSTRAINT(NAME, IDX_LIST)
%   [A, L, U, VS] = OM.PARAMS_LIN_CONSTRAINT(...)
%   [A, L, U, VS, I1, IN] = OM.PARAMS_LIN_CONSTRAINT(...)
%
%   With no input parameters, it assembles and returns the parameters
%   for the aggregate linear constraints from all linear constraint sets
%   added using ADD_LIN_CONSTRAINT. The values of these parameters are
%   cached for subsequent calls. The parameters are A, L and U where the
%   linear constraint is of the form
%       L <= A * x <= U
%
%   If a NAME is provided then it simply returns the parameters for the
%   corresponding named set. Likewise for indexed named sets specified
%   by NAME and IDX_LIST. An optional 4th output argument VS indicates the
%   variable sets used by this constraint set. The size of A will be
%   consistent with VS. Optional 5th and 6th output arguments I1 and IN
%   indicate the starting and ending row indices of the corresponding
%   constraint set in the full aggregate constraint matrix.
%
%   Examples:
%       [A, l, u] = om.params_lin_constraint();
%       [A, l, u, vs, i1, iN] = om.params_lin_constraint('Pmis');
%
%   See also OPT_MODEL, ADD_LIN_CONSTRAINT.

%   MP-Opt-Model
%   Copyright (c) 2008-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if nargin > 1       %% individual set
    if nargin < 3
        idx = {};
    end
    if isempty(idx)
        if numel(om.lin.idx.i1.(name)) == 1     %% simple named set
            A = om.lin.data.A.(name);
            l = om.lin.data.l.(name);
            u = om.lin.data.u.(name);
            if nargout > 3
                vs = om.lin.data.vs.(name);
                if nargout > 5
                    i1 = om.lin.idx.i1.(name);      %% starting row index
                    iN = om.lin.idx.iN.(name);      %% ending row index
                end
            end
        else                                    %% indexing required
            error('@opt_model/params_lin_constraint: linear constraint set ''%s'' requires an IDX_LIST arg', name);
        end
    else
        % (calls to substruct() are relatively expensive ...
        % s = substruct('.', name, '{}', idx);
        % ... so replace it with these more efficient lines)
        sc = struct('type', {'.', '{}'}, 'subs', {name, idx});
        A = subsref(om.lin.data.A, sc);
        l = subsref(om.lin.data.l, sc);
        u = subsref(om.lin.data.u, sc);
        if nargout > 3
            vs = subsref(om.lin.data.vs, sc);
            if nargout > 5
                sn = sc; sn(2).type = '()';         %% num array field
                i1 = subsref(om.lin.idx.i1, sn);    %% starting row index
                iN = subsref(om.lin.idx.iN, sn);    %% ending row index
            end
        end
    end
else                %% aggregate
    cache = om.lin.params;
    if isempty(cache)       %% build the aggregate
        nx = om.var.N;          %% number of variables
        nlin = om.lin.N;        %% number of linear constraints
        if om.lin.NS < 25 || om.lin.NS < 100 && nx < 300
            %% METHOD 1: Add sparse matrices (original method)
            At = sparse(nx, nlin);  %% transpose of constraint matrix
            u = Inf(nlin, 1);       %% upper bound
            l = -u;                 %% lower bound

            %% fill in each piece
            for k = 1:om.lin.NS
                name = om.lin.order(k).name;
                idx  = om.lin.order(k).idx;
                [Ak, lk, uk, vs, i1, iN] = om.params_lin_constraint(name, idx);
                [mk, nk] = size(Ak);        %% size of Ak
                if mk
                    Akt_full = sparse(nx, nlin);
                    if isempty(vs)
                        if nk == nx     %% full size
                            Akt_full(:, i1:iN) = Ak';
                        else            %% vars added since adding this constraint set
                            Ak_all_cols = sparse(mk, nx);
                            Ak_all_cols(:, 1:nk) = Ak;
                            Akt_full(:, i1:iN) = Ak_all_cols';
                        end
                    else
                        jj = om.varsets_idx(vs);    %% indices for var set
                        Ak_all_cols = sparse(mk, nx);
                        Ak_all_cols(:, jj) = Ak;
                        Akt_full(:, i1:iN) = Ak_all_cols';
                    end
                    At = At + Akt_full;
                    l(i1:iN) = lk;
                    u(i1:iN) = uk;
                end
            end
            A = At';
        else
            %% METHOD 2: construct using single call to sparse()
            A_ijv = cell(om.lin.NS, 3); %% indices/values to construct A
            u = Inf(nlin, 1);       %% upper bound
            l = -u;                 %% lower bound

            %% fill in each piece
            for k = 1:om.lin.NS
                name = om.lin.order(k).name;
                idx  = om.lin.order(k).idx;
                [Ak, lk, uk, vs, i1, iN] = om.params_lin_constraint(name, idx);
                [mk, nk] = size(Ak);        %% size of Ak
                if mk
                    % find nonzero sub indices and values
                    [i, j, v] = find(Ak);
                    if mk == 1  %% force col vectors for single row Ak
                        i = i'; j = j'; v = v';
                    end

                    if isempty(vs)
                        A_ijv(k,:) = {i+(i1-1), j, v};
                    else
                        jj = om.varsets_idx(vs)';    %% indices for var set
                        A_ijv(k,:) = {i+(i1-1), jj(j), v};
                    end
                    l(i1:iN) = lk;
                    u(i1:iN) = uk;
                end
            end
            A = sparse( vertcat(A_ijv{:,1}), ...
                        vertcat(A_ijv{:,2}), ...
                        vertcat(A_ijv{:,3}), nlin, nx);
        end

        %% cache aggregated parameters
        om.lin.params = struct('A', A, 'l', l, 'u', u);
    else                    %% return cached values
        A = cache.A;
        l = cache.l;
        u = cache.u;
    end
    if nargout > 3
        vs = {};
    end
end
