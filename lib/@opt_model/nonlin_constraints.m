function [g, dg] = nonlin_constraints(om, x, iseq)
%NONLIN_CONSTRAINTS  Builds and returns full set of nonlinear constraints.
%   [G, DG] = OM.NONLIN_CONSTRAINTS(X, ISEQ)
%   Builds a full set of nonlinear equality or inequality constraints and
%   their gradients for a given value of the optimization vector X,
%   based on constraints added by ADD_CONSTRAINTS.
%
%       g(X) = 0
%       h(X) <= 0
%
%   Example:
%       [g, dg] = om.nonlin_constraints(x, 1)
%       [h, dh] = om.nonlin_constraints(x, 0)
%
%   See also OPT_MODEL, ADD_CONSTRAINTS.

%   MATPOWER
%   Copyright (c) 2008-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% get constraint type
if iseq         %% equality constraints
    ff = 'nle';
else            %% inequality constraints
    ff = 'nli';
end

%% initialize g, dg
g = NaN(om.(ff).N, 1);
dgt = sparse(om.var.N, om.(ff).N);  %% use transpose of dg for speed

%% fill in each piece
s = struct('type', {'.', '()'}, 'subs', {'', 1});
s1 = s;
s2 = s;
s2(2).type = '{}';
for k = 1:om.(ff).NS
    name = om.(ff).order(k).name;
    idx  = om.(ff).order(k).idx;
    if isempty(idx)
        N = om.(ff).idx.N.(name);
    else
        % (calls to substruct() are relatively expensive ...
        % s1 = substruct('.', name, '()', idx);
        % s2 = substruct('.', name, '{}', idx);
        % ... so replace them with these more efficient lines)
        s1(1).subs = name;
        s1(2).subs = idx;
        s2(1).subs = name;
        s2(2).subs = idx;
        N = subsref(om.(ff).idx.N, s1);
    end
    if N                                %% non-zero number of rows
        if isempty(idx)
            fcn = om.(ff).data.fcn.(name);      %% fcn for kth constraint set
            i1 = om.(ff).idx.i1.(name);         %% starting row index
            iN = om.(ff).idx.iN.(name);         %% ending row index
            vsl = om.(ff).data.vs.(name);       %% var set list
        else
            fcn = subsref(om.(ff).data.fcn, s2);%% fcn for kth constraint set
            i1 = subsref(om.(ff).idx.i1, s1);   %% starting row index
            iN = subsref(om.(ff).idx.iN, s1);   %% ending row index
            vsl = subsref(om.(ff).data.vs, s2); %% var set list
        end
        if isempty(vsl)         %% all rows of x
            xx = x;
        else                    %% selected rows of x
            xx = cell(size(vsl));
            for v = 1:length(vsl)
                % (calls to substruct() are relatively expensive ...
                % s = substruct('.', vsl(v).name, '()', vsl(v).idx);
                % ... so replace it with these more efficient lines)
                s(1).subs = vsl(v).name;
                s(2).subs = vsl(v).idx;
                j1 = subsref(om.var.idx.i1, s); %% starting row in full x
                jN = subsref(om.var.idx.iN, s); %% ending row in full x
                xx{v} = x(j1:jN);
            end
        end
        [gk, dgk] = fcn(xx);    %% evaluate kth constraint and gradient
        g(i1:iN) = gk;          %% assign kth constraint
        
        if isempty(vsl)         %% all rows of x
            if size(dgk, 2) == om.var.N
                dgt(:, i1:iN) = dgk';     %% assign as columns in transpose for speed
            else                %% must have added vars since adding
                                %% this constraint set
                dgt(1:size(dgk, 2), i1:iN) = dgk';  %% assign as columns in transpose for speed
            end
        else                    %% selected rows of x
            kN = 0;                             %% initialize last col of dgk used
            dgi = sparse(N, om.var.N);
            for v = 1:length(vsl)
                % (calls to substruct() are relatively expensive ...
                % s = substruct('.', vsl(v).name, '()', vsl(v).idx);
                % ... so replace it with these more efficient lines)
                s(1).subs = vsl(v).name;
                s(2).subs = vsl(v).idx;
                j1 = subsref(om.var.idx.i1, s); %% starting row in full x
                jN = subsref(om.var.idx.iN, s); %% ending row in full x
                k1 = kN + 1;                    %% starting column in dgk
                kN = kN + subsref(om.var.idx.N, s);%% ending column in dgk
                dgi(:, j1:jN) = dgk(:, k1:kN);
            end
            dgt(:, i1:iN) = dgi';     %% assign as columns in transpose for speed
        end
    end
end
dg = dgt';
