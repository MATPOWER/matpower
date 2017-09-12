function [f, df, d2f] = eval_quad_cost(om, x, name, idx)
%EVAL_QUAD_COST  Evaluates individual or full set of quadratic costs.
%   F = OM.EVAL_QUAD_COST(X ...)
%   [F, DF] = OM.EVAL_QUAD_COST(X ...)
%   [F, DF, D2F] = OM.EVAL_QUAD_COST(X ...)
%   [F, DF, D2F] = OM.EVAL_QUAD_COST(X, NAME)
%   [F, DF, D2F] = OM.EVAL_QUAD_COST(X, NAME, IDX)
%   Evaluates an individual named set or the full set of quadratic
%   costs and their derivatives for a given value of the optimization vector
%   X, based on costs added by ADD_QUAD_COST.
%
%   Example:
%       [f, df, d2f] = om.eval_quad_cost(x)
%       [f, df, d2f] = om.eval_quad_cost(x, name)
%       [f, df, d2f] = om.eval_quad_cost(x, name, idx)
%
%   See also OPT_MODEL, ADD_QUAD_COST, PARAMS_QUAD_COST.

%   MATPOWER
%   Copyright (c) 2008-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if om.qdc.N
    done = 0;

    %% collect cost parameters
    if nargin < 3                       %% full set
        [Q, c, k, vs] = om.params_quad_cost();
        N = 1;
    elseif nargin < 4 || isempty(idx)   %% name, no idx provided
        dims = size(om.qdc.idx.i1.(name));
        if prod(dims) == 1              %% simple named set
            [Q, c, k, vs] = om.params_quad_cost(name);
            N = om.getN('qdc', name);
        elseif nargout == 1             %% indexing required, recurse
            f = 0;          %% initialize cumulative cost
            idx = num2cell(ones(size(dims))); %% initialize idx
            while ~done     %% call eval_quad_cost() recursively
                f = f + sum(om.eval_quad_cost(x, name, idx));
            
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
            error('@opt_model/eval_quad_cost: quadratic cost set ''%s'' requires an IDX arg when requesting DF output', name)
        end
    else                                %% indexed named set
        [Q, c, k, vs] = om.params_quad_cost(name, idx);
        N = om.getN('qdc', name, idx);
    end
    
    if ~done
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
        end     %% nargout > 1
    end         %% ~done
else
    f = 0;
    if nargout > 1
%         nx = length(x);
%         df = zeros(nx, 1);
        df = [];
        if nargout > 2
%             d2f = sparse(nx, nx);
            d2f = [];
        end
    end
end
