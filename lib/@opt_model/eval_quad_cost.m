function [f, df, d2f] = eval_quad_cost(om, x, name, idx)
%EVAL_QUAD_COST  Evaluates individual or full set of quadratic costs.
%   F = OM.EVAL_QUAD_COST(X ...)
%   [F, DF] = OM.EVAL_QUAD_COST(X ...)
%   [F, DF, D2F] = OM.EVAL_QUAD_COST(X ...)
%   [F, DF, D2F] = OM.EVAL_QUAD_COST(X, NAME)
%   [F, DF, D2F] = OM.EVAL_QUAD_COST(X, NAME, IDX)
%   Evaluates an individual named set or the full set of quadratic
%   costs and their derivatives for a given value of the optimization vector
%   X, based on costs added by ADD_QUADRATIC_COSTS.
%
%   Example:
%       [f, df, d2f] = om.eval_quad_cost(x)
%       [f, df, d2f] = om.eval_quad_cost(x, name)
%       [f, df, d2f] = om.eval_quad_cost(x, name, idx)
%
%   See also OPT_MODEL, ADD_QUADRATIC_COSTS.

%   MATPOWER
%   Copyright (c) 2008-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if om.qdc.N
    %% collect cost parameters
    if nargin > 2       %% individual set
        if nargin < 4 || isempty(idx)   %% simple named set
            [Q, c, k, vs] = om.params_quad_cost(name);
            N = om.getN('qdc', name);
        else                            %% indexed named set
            [Q, c, k, vs] = om.params_quad_cost(name, idx);
            N = om.getN('qdc', name, idx);
        end
    else                %% full set
        [Q, c, k, vs] = om.params_quad_cost();
        N = 1;
    end

    %% assemble appropriately-sized x vector
    xx = om.varsets_x(x, vs, 'vector');
    
    %% compute/assemble f
    if N == 1               %% f is scalar (Q is matrix, k is scalar)
        f = k;                  %% start with k term
        if ~isempty(c)
            f = f + c'*xx;      %% add c term
        end
        if ~isempty(Q)          %% add Q term
            f = f + (xx'*Q*xx)/2;
        end
    else                    %% f is vector (Q is vector, k is vector or 0)
        if isempty(c)           %% Q, k terms only
            f = (Q .* xx.^2)/2 + k;
        else
            if isempty(Q)       %% c, k terms only
                f = c .* xx + k;
            else                %% Q, c, k terms
                f = (Q .* xx.^2)/2 + c .* xx + k;
            end
        end
    end

    if nargout > 1
        %% compute/assemble df
        if ~isempty(c)
            df = c;             %% start with c term
        else
            df = 0;             %% start with nothing
        end
        if ~isempty(Q)
            if N == 1           %% f is scalar (Q is matrix, k is scalar)
                df = df + Q*xx;     %% add Q term
            else                %% f is vector (Q is vector, k is vector or 0)
                df = df + Q.*xx;    %% add Q term
            end
        end

        %% assemble d2f
        if nargout > 2
            if isempty(Q)
                nx = length(xx);
                if N == 1   %% f is scalar (Q is matrix, k is scalar)
                    d2f = sparse(nx, nx);
                else        %% f is vector (Q is vector, k is vector or 0)
                    d2f = sparse(nx, 1);
                end
            else
                d2f = Q;
            end
        end
    end
else
    f = [];
    df = [];
    d2f = [];
end
