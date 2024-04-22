function [f, df, d2f] = eval_nln_cost(om, x, name, idx)
% eval_nln_cost - Evaluates individual or full set of general nonlinear costs.
% ::
%
%   [F, DF, D2F] = OM.EVAL_NLN_COST(X)
%   [F, DF, D2F] = OM.EVAL_NLN_COST(X, NAME)
%   [F, DF, D2F] = OM.EVAL_NLN_COST(X, NAME, IDX_LIST)
%   Evaluates the cost function and its derivatives for an individual named
%   set or the full set of general nonlinear costs for a given value of the
%   optimization vector X, based on costs added by ADD_NLN_COST.
%
%   Example:
%       [f, df, d2f] = om.eval_nln_cost(x)
%       [f, df, d2f] = om.eval_nln_cost(x, name)
%       [f, df, d2f] = om.eval_nln_cost(x, name, idx_list)
%
% See also opt_model, add_nln_cost, params_nln_cost.

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% initialize
if nargin < 4
    idx = {};
end
nlc = om.nlc;

if nargin < 3                       %% full set
    f = 0;
    nx = om.var.N;          %% number of variables
    if nargout > 1
        df = zeros(nx, 1);  %% column vector
        if nargout > 2
            d2f = sparse(nx, nx);   %% symmetric, so transpose (used for speed)
        end                         %% is equal
    end

    for k = 1:nlc.NS
        name = nlc.order(k).name;
        idx  = nlc.order(k).idx;
        [N, fcn, vs] = om.params_nln_cost(name, idx);
        if N ~= 1
            error('opt_model.eval_nln_cost: not yet implemented for vector valued functions');
        end
        xx = om.varsets_x(x, vs);
        if nargout == 3
            [fk, dfk, d2fk] = fcn(xx);  %% evaluate kth cost, gradient, Hessian
        elseif nargout == 2
            [fk, dfk] = fcn(xx);        %% evaluate kth cost and gradient
        else
            fk = fcn(xx);               %% evaluate kth cost
        end

        f = f + fk;
        
        if nargout > 1              %% assemble gradient
            nk = length(dfk);
            if isempty(vs)          %% all rows of x
                if nk == nx
                    df = df + dfk;
                    if nargout > 2  %% assemble Hessian
                        d2f = d2f + d2fk;
                    end
                else                %% must have added vars since adding
                                    %% this cost set
                    df(1:nk) = df(1:nk) + dfk;
                    if nargout > 2      %% assemble Hessian
                        d2fk_all_cols = sparse(nk, nx);
                        d2fk_all_cols(:, 1:nk) = d2fk;
                        d2fk_full = sparse(nx, nx);
                        d2f(:, 1:nk) = d2f(:, 1:nk) + d2fk_all_cols';
                    end
                end
            else                    %% selected rows of x
                jj = om.varsets_idx(vs);    %% indices for var set
                df(jj) = df(jj) + dfk;
                if nargout > 2      %% assemble Hessian
                    d2fk_all_cols = sparse(nk, nx);
                    d2fk_all_cols(:, jj) = d2fk;
                    d2fk_full = sparse(nx, nx);
                    d2f(:, jj) = d2f(:, jj) + d2fk_all_cols';
                end
            end
        end
    end
else                                %% individual named set
    dims = size(nlc.idx.i1.(name));
    if ~isempty(idx) || prod(dims) == 1 %% indexed, or simple named set
        [N, fcn, vs] = om.params_nln_cost(name, idx);
        if N ~= 1
            error('opt_model.eval_nln_cost: not yet implemented for vector valued functions');
        end
        xx = om.varsets_x(x, vs);
        if nargout == 3
            [f, df, d2f] = fcn(xx);     %% evaluate kth cost, gradient, Hessian
        elseif nargout == 2
            [f, df] = fcn(xx);          %% evaluate kth cost and gradient
        else
            f = fcn(xx);                %% evaluate kth cost
        end
    elseif nargout == 1             %% sum over all indices for name
        done = 0;
        f = 0;          %% initialize cumulative cost
        idx = num2cell(ones(size(dims))); %% initialize idx
        while ~done     %% call eval_nln_cost() recursively
            f = f + sum(om.eval_nln_cost(x, name, idx));
        
            %% increment idx
            D = length(dims);
            idx{D} = idx{D} + 1;    %% increment last dimension
            for d = D:-1:2          %% increment next dimension, if necessary
                if idx{d} > dims(d)
                    idx{d} = 1;
                    idx{d-1} = idx{d-1} + 1;
                end
            end
            if idx{1} > dims(1)     %% check if done
                done = 1;
            end
        end
    else
        error('opt_model.eval_nln_cost: general nonlinear cost set ''%s'' requires an IDX_LIST arg when requesting DF output', name)
    end
end
