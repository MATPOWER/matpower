function [f, df, d2f] = eval_legacy_cost(om, x, name, idx)
%EVAL_LEGACY_COST  Evaluates individual or full set of legacy user costs.
%   F = OM.EVAL_LEGACY_COST(X ...)
%   [F, DF] = OM.EVAL_LEGACY_COST(X ...)
%   [F, DF, D2F] = OM.EVAL_LEGACY_COST(X ...)
%   [F, DF, D2F] = OM.EVAL_LEGACY_COST(X, NAME)
%   [F, DF, D2F] = OM.EVAL_LEGACY_COST(X, NAME, IDX)
%   Evaluates an individual named set or the full set of legacy user
%   costs and their derivatives for a given value of the optimization vector
%   X, based on costs added by ADD_LEGACY_COST.
%
%   Example:
%       [f, df, d2f] = om.eval_legacy_cost(x)
%       [f, df, d2f] = om.eval_legacy_cost(x, name)
%       [f, df, d2f] = om.eval_legacy_cost(x, name, idx)
%
%   See also OPT_MODEL, ADD_LEGACY_COST, PARAMS_LEGACY_COST.

%   MATPOWER
%   Copyright (c) 2008-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if om.cost.N
    done = 0;

    %% collect cost parameters
    if nargin < 3                       %% full set
        [cp, vs] = om.params_legacy_cost();
    elseif nargin < 4 || isempty(idx)   %% name, no idx provided
        dims = size(om.cost.idx.i1.(name));
        if prod(dims) == 1              %% simple named set
            [cp, vs] = om.params_legacy_cost(name);
        elseif nargout == 1             %% indexing required, recurse
            f = 0;          %% initialize cumulative cost
            idx = num2cell(ones(size(dims))); %% initialize idx
            while ~done     %% call eval_legacy_cost() recursively
                f = f + om.eval_legacy_cost(x, name, idx);
            
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
            error('@opt_model/eval_legacy_cost: legacy cost set ''%s'' requires an IDX arg when requesting DF output', name)
        end
    else                                %% indexed named set
        [cp, vs] = om.params_legacy_cost(name, idx);
    end

    if ~done
        %% assemble appropriately-sized x vector
        xx = om.varsets_x(x, vs, 'vector');

        %% compute function & derivatives
        if nargout == 1
            f = opf_legacy_user_cost_fcn(xx, cp);
        elseif nargout == 2
            [f, df] = opf_legacy_user_cost_fcn(xx, cp);
        else    %% nargout == 3
            [f, df, d2f] = opf_legacy_user_cost_fcn(xx, cp);
        end
    end
else
    f = 0;
    if nargout > 1
        df = [];
        if nargout > 2
            d2f = [];
        end
    end
end
