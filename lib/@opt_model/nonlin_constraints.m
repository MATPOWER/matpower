function [g, dg] = nonlin_constraints(om, x, iseq)
%NONLIN_CONSTRAINTS  Builds and returns full set of nonlinear constraints.
%   [G, DG] = OM.NONLIN_CONSTRAINTS(X, ISEQ)
%   Builds a full set of nonlinear equality or inequality constraints and
%   their gradients for a given value of the optimization vector X,
%   based on constraints added by ADD_NLN_CONSTRAINTS.
%
%       g(X) = 0
%       h(X) <= 0
%
%   Example:
%       [g, dg] = om.nonlin_constraints(x, 1)
%       [h, dh] = om.nonlin_constraints(x, 0)
%
%   See also OPT_MODEL, ADD_NLN_CONSTRAINTS, NONLIN_CONSTRAINT_HESS.

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
dg = sparse(0, om.var.N);   %% build gradient by stacking

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
            vs = om_nlx.data.vs.(name);         %% var sets
        else
            fcn = subsref(om_nlx.data.fcn, s2); %% fcn for kth constraint set
            i1 = subsref(om_nlx.idx.i1, s1);    %% starting row index
            iN = subsref(om_nlx.idx.iN, s1);    %% ending row index
            vs = subsref(om_nlx.data.vs, s2);   %% var sets
        end
        xx = om.varsets_x(x, vs);
        [gk, dgk] = fcn(xx);    %% evaluate kth constraint and gradient
        g(i1:iN) = gk;          %% assign kth constraint
        
        if isempty(vs)          %% all rows of x
            if size(dgk, 2) == om.var.N
                dg = [dg; dgk];
            else                %% must have added vars since adding
                                %% this constraint set
                dg(i1:iN, 1:size(dgk, 2)) = dgk;
            end
        else                    %% selected rows of x
            jj = om.varsets_idx(vs);    %% column indices for var set
            dgi = sparse(N, om.var.N);
            dgi(:, jj) = dgk;
            dg = [dg; dgi];
        end
    end
end
