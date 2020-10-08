function [Ax_u, l_Ax, A] = eval_lin_constraint(om, x, name, idx)
%EVAL_LIN_CONSTRAINT  Builds and returns full set of linear constraints.
%   AX_U = OM.EVAL_LIN_CONSTRAINT(X)
%   AX_U = OM.EVAL_LIN_CONSTRAINT(X, NAME)
%   AX_U = OM.EVAL_LIN_CONSTRAINT(X, NAME, IDX_LIST)
%   [AX_U, L_AX] = OM.EVAL_LIN_CONSTRAINT(...)
%   [AX_U, L_AX, A] = OM.EVAL_LIN_CONSTRAINT(...)
%   Builds the linear constraints for the full set of constraints or an
%   individual named subset for a given value of the vector X, based on
%   constraints added by ADD_LIN_CONSTRAINT.
%
%       l <= A * x <= u
%
%   Returns A*X - U, and optionally L - A*X and A.
%
%   Example:
%       [Ax_u, l_Ax, A] = om.eval_lin_constraint(x)
%
%   See also OPT_MODEL, ADD_LIN_CONSTRAINT, PARAMS_LIN_CONSTRAINT.

%   MP-Opt-Model
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if om.lin.N
    %% collect cost parameters
    if nargin < 3                       %% full set
        [A, l, u, vs] = om.params_lin_constraint();
        N = 1;
    elseif nargin < 4 || isempty(idx)   %% name, no idx provided
        dims = size(om.lin.idx.i1.(name));
        if prod(dims) == 1              %% simple named set
            [A, l, u, vs] = om.params_lin_constraint(name);
            N = om.getN('lin', name);
        else
            error('@opt_model/eval_quad_cost: quadratic cost set ''%s'' requires an IDX_LIST arg when requesting DF output', name)
        end
    else                                %% indexed named set
        [A, l, u, vs] = om.params_lin_constraint(name, idx);
        N = om.getN('lin', name, idx);
    end
    
    %% assemble appropriately-sized x vector
    xx = om.varsets_x(x, vs, 'vector');

    %% compute constraints
    Ax = A * xx;
    Ax_u = Ax - u;
    if nargout > 1
        l_Ax = l - Ax;
    end
else
    Ax_u = [];
    if nargout > 1
        l_Ax = [];
    end
end
