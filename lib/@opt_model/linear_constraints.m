function [A, l, u] = linear_constraints(om)
%LINEAR_CONSTRAINTS  Builds and returns the full set of linear constraints.
%   [A, L, U] = LINEAR_CONSTRAINTS(OM)
%   Builds the full set of linear constraints based on those added by
%   ADD_CONSTRAINTS.
%
%       L <= A * x <= U
%
%   Example:
%       [A, l, u] = linear_constraints(om);
%
%   See also OPT_MODEL, ADD_CONSTRAINTS.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2008-2012 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.


%% initialize A, l and u
nnzA = 0;
for k = 1:om.lin.NS
    name = om.lin.order(k).name;
    idx  = om.lin.order(k).idx;
    if isempty(idx)
        nnzA = nnzA + nnz(om.lin.data.A.(name));
    else
        s = substruct('.', name, '{}', idx);
        nnzA = nnzA + nnz(subsref(om.lin.data.A, s));
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
        s1 = substruct('.', name, '()', idx);
        s2 = substruct('.', name, '{}', idx);
        N = subsref(om.lin.idx.N, s1);
    end
    if N                                %% non-zero number of rows to add
        if isempty(idx)
            Ak = om.lin.data.A.(name);          %% A for kth linear constrain set
            i1 = om.lin.idx.i1.(name);          %% starting row index
            iN = om.lin.idx.iN.(name);          %% ending row index
            vsl = om.lin.data.vs.(name);        %% var set list
        else
            Ak = subsref(om.lin.data.A, s2);    %% A for kth linear constrain set
            i1 = subsref(om.lin.idx.i1, s1);    %% starting row index
            iN = subsref(om.lin.idx.iN, s1);    %% ending row index
            vsl = subsref(om.lin.data.vs, s2);  %% var set list
        end
        if isempty(vsl)         %% full rows
            if size(Ak,2) == om.var.N
                At(:, i1:iN) = Ak';     %% assign as columns in transpose for speed
            else                %% must have added vars since adding
                                %% this constraint set
                At(1:size(Ak,2), i1:iN) = Ak';  %% assign as columns in transpose for speed
            end
        else                    %% selected columns
            kN = 0;                             %% initialize last col of Ak used
            Ai = sparse(N, om.var.N);
            for v = 1:length(vsl)
                s = substruct('.', vsl(v).name, '()', vsl(v).idx);
                j1 = subsref(om.var.idx.i1, s); %% starting column in A
                jN = subsref(om.var.idx.iN, s); %% ending column in A
                k1 = kN + 1;                    %% starting column in Ak
                kN = kN + subsref(om.var.idx.N, s);%% ending column in Ak
                Ai(:, j1:jN) = Ak(:, k1:kN);
            end
            At(:, i1:iN) = Ai';     %% assign as columns in transpose for speed
        end

        if isempty(idx)
            l(i1:iN) = om.lin.data.l.(name);
            u(i1:iN) = om.lin.data.u.(name);
        else
            l(i1:iN) = subsref(om.lin.data.l, s2);
            u(i1:iN) = subsref(om.lin.data.u, s2);
        end
    end
end
A = At';
