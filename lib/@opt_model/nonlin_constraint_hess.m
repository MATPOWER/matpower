function d2G = nonlin_constraint_hess(om, x, lam, iseq)
%NONLIN_CONSTRAINT_HESS  Builds and returns Hessian of nonlinear constraints.
%   D2G = OM.NONLIN_CONSTRAINT_HESS(X, LAM, ISEQ)
%   Builds the Hessian of the full set of nonlinear equality or inequality
%   constraints for given values of the optimization vector X and dual
%   variables LAM, based on constraints added by ADD_CONSTRAINTS.
%
%       g(X) = 0
%       h(X) <= 0
%
%   Example:
%       d2G = om.nonlin_constraint_hess(x, lam, 1)
%       d2H = om.nonlin_constraint_hess(x, lam, 0)
%
%   See also OPT_MODEL, ADD_CONSTRAINTS, NONLIN_CONSTRAINTS.

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
d2G = sparse(om.var.N, om.var.N);

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
            d2G_fcn = om.(ff).data.hess.(name); %% Hessian fcn for kth constraint set
            i1 = om.(ff).idx.i1.(name);         %% starting row index
            iN = om.(ff).idx.iN.(name);         %% ending row index
            vsl = om.(ff).data.vs.(name);       %% var set list
        else
            d2G_fcn = subsref(om.(ff).data.hess, s2);  %% Hessian fcn for kth constraint set
            i1 = subsref(om.(ff).idx.i1, s1);    %% starting row index
            iN = subsref(om.(ff).idx.iN, s1);    %% ending row index
            vsl = subsref(om.(ff).data.vs, s2);  %% var set list
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
        d2Gk = d2G_fcn(xx, lam(i1:iN));     %% evaluate kth Hessian
        
        if isempty(vsl)         %% all rows of x
            Nk = size(d2Gk, 2);
            if Nk == om.var.N
                d2Gk_full = d2Gk;
            else                %% must have added vars since adding
                                %% this constraint set
                d2Gk_full = sparse(om.var.N, om.var.N);
                d2Gk_full(1:Nk, 1:Nk) = d2Gk;
            end
        else                    %% selected rows of x
            kN = 0;                             %% initialize last col of d2Gk used
            dgi = sparse(N, om.var.N);
            ii = [];
            for v = 1:length(vsl)
                % (calls to substruct() are relatively expensive ...
                % s = substruct('.', vsl(v).name, '()', vsl(v).idx);
                % ... so replace it with these more efficient lines)
                s(1).subs = vsl(v).name;
                s(2).subs = vsl(v).idx;
                j1 = subsref(om.var.idx.i1, s); %% starting row in full x
                jN = subsref(om.var.idx.iN, s); %% ending row in full x
                ii = [ii; (j1:jN)'];
            end
            d2Gk_full = sparse(om.var.N, om.var.N);
            d2Gk_full(ii, ii) = d2Gk;
        end
        d2G = d2G + d2Gk_full;
    end
end
