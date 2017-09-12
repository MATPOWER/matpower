function [A, l, u, vs, i1, iN] = params_lin_constraint(om, name, idx)
%PARAMS_LIN_CONSTRAINT  Builds and returns linear constraint parameters.
%   [A, L, U] = OM.PARAMS_LIN_CONSTRAINT()
%   [A, L, U] = OM.PARAMS_LIN_CONSTRAINT(NAME)
%   [A, L, U] = OM.PARAMS_LIN_CONSTRAINT(NAME, IDX)
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
%   by NAME and IDX. 
%
%   An optional 4th output argument VS indicates the variable sets used by
%   this cost set. The size of A will be consistent with VS.
%
%   If NAME is provided, optional 5th and 6th output arguments I1 and IN
%   indicate the starting and ending row indices of the corresponding
%   constraint set in the full aggregate constraint matrix.
%
%   Examples:
%       [A, l, u] = om.params_lin_constraint();
%       [A, l, u, vs, i1, i2] = om.params_lin_constraint('Pmis');
%
%   See also OPT_MODEL, ADD_LIN_CONSTRAINT.

%   MATPOWER
%   Copyright (c) 2008-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin > 1       %% individual set
    if nargin < 3
        idx = {};
    end
    if isempty(idx)
        if prod(size(om.lin.idx.i1.(name))) == 1
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
        else
            error('@opt_model/params_lin_constraint: linear constraint set ''%s'' requires an IDX arg', name);
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
                    else            %% vars added since adding this cost set
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
