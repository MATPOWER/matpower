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
    om_nlx = om.nle;
else            %% inequality constraints
    om_nlx = om.nli;
end

%% initialize g, dg
g = NaN(om_nlx.N, 1);
dgt = sparse(om.var.N, om_nlx.N);  %% use transpose of dg for speed

%% fill in each piece
s = struct('type', {'.', '()'}, 'subs', {'', 1});
s1 = s;
s2 = s;
s2(2).type = '{}';
for k = 1:om_nlx.NS
    name = om_nlx.order(k).name;
    idx  = om_nlx.order(k).idx;
    if isempty(idx)
        if ~isfield(om_nlx.data.fcn, name)
            continue;   %% skip, there is no function handle stored here,
                        %% the function value for this named set was included
                        %% in the value computed by a previous named set
        end
        N = om_nlx.idx.N.(name);    %% number of constraint functions
                                    %% evaluated for this named set
        if isfield(om_nlx.data.include, name)
            N = N + sum(om_nlx.data.include.(name).N);
        end
    else
        % (calls to substruct() are relatively expensive ...
        % s1 = substruct('.', name, '()', idx);
        % s2 = substruct('.', name, '{}', idx);
        % ... so replace them with these more efficient lines)
        s1(1).subs = name;
        s1(2).subs = idx;
        s2(1).subs = name;
        s2(2).subs = idx;
        N = subsref(om_nlx.idx.N, s1);
    end
    if N                                %% non-zero number of rows
        if isempty(idx)
            fcn = om_nlx.data.fcn.(name);       %% fcn for kth constraint set
            i1 = om_nlx.idx.i1.(name);          %% starting row index
            iN = i1 + N - 1;                    %% ending row index
            vsl = om_nlx.data.vs.(name);        %% var set list
        else
            fcn = subsref(om_nlx.data.fcn, s2); %% fcn for kth constraint set
            i1 = subsref(om_nlx.idx.i1, s1);    %% starting row index
            iN = subsref(om_nlx.idx.iN, s1);    %% ending row index
            vsl = subsref(om_nlx.data.vs, s2);  %% var set list
        end
        if isempty(vsl)         %% all rows of x
            xx = x;
        else                    %% selected rows of x
            xx = cell(size(vsl));
            for v = 1:length(vsl)
                vidx = vsl(v).idx;
                if isempty(vidx)
                    j1 = om.var.idx.i1.(vsl(v).name);
                    jN = om.var.idx.iN.(vsl(v).name);
                else
                    % (calls to substruct() are relatively expensive ...
                    % s = substruct('.', vsl(v).name, '()', vsl(v).idx);
                    % ... so replace it with these more efficient lines)
                    s(1).subs = vsl(v).name;
                    s(2).subs = vsl(v).idx;
                    j1 = subsref(om.var.idx.i1, s); %% starting row in full x
                    jN = subsref(om.var.idx.iN, s); %% ending row in full x
                end
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
                vidx = vsl(v).idx;
                if isempty(vidx)
                    j1 = om.var.idx.i1.(vsl(v).name);
                    jN = om.var.idx.iN.(vsl(v).name);
                    k1 = kN + 1;                    %% starting column in dgk
                    kN = kN + om.var.idx.N.(vsl(v).name);%% ending column in dgk
                else
                    % (calls to substruct() are relatively expensive ...
                    % s = substruct('.', vsl(v).name, '()', vsl(v).idx);
                    % ... so replace it with these more efficient lines)
                    s(1).subs = vsl(v).name;
                    s(2).subs = vsl(v).idx;
                    j1 = subsref(om.var.idx.i1, s); %% starting row in full x
                    jN = subsref(om.var.idx.iN, s); %% ending row in full x
                    k1 = kN + 1;                    %% starting column in dgk
                    kN = kN + subsref(om.var.idx.N, s);%% ending column in dgk
                end
                dgi(:, j1:jN) = dgk(:, k1:kN);
            end
            dgt(:, i1:iN) = dgi';     %% assign as columns in transpose for speed
        end
    end
end
dg = dgt';
