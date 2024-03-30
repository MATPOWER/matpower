function om = init_indexed_name(om, set_type, name, dim_list)
% init_indexed_name - Initializes the dimensions for an indexed named set.
% ::
%
%   OM.INIT_INDEXED_NAME(SET_TYPE, NAME, DIM_LIST)
%
%   Initializes the dimensions for an indexed named variable, constraint
%   or cost set.
%
%   Variables, constraints and costs are referenced in OPT_MODEL in terms
%   of named sets. The specific type of named set being referenced is
%   given by SET_TYPE, with the following valid options:
%       SET_TYPE = 'var'   => variable set
%       SET_TYPE = 'lin'   => linear constraint set
%       SET_TYPE = 'nle'   => nonlinear equality constraint set
%       SET_TYPE = 'nli'   => nonlinear inequality constraint set
%       SET_TYPE = 'qdc'   => quadratic cost set
%       SET_TYPE = 'nlc'   => nonlinear cost set
%       SET_TYPE = 'cost'  => legacy cost set
%
%   Indexed Named Sets
%
%   A variable, constraint or cost set can be identified by a single NAME,
%   such as 'Pmismatch', or by a name that is indexed by one or more indices,
%   such as 'Pmismatch(3,4)'. For an indexed named set, before adding the
%   indexed variable, constraint or cost sets themselves, the dimensions of
%   the indexed set must be set by calling INIT_INDEXED_NAME, where
%   DIM_LIST is a cell array of the dimensions.
%
%   Examples:
%       %% linear constraints with indexed named set 'R(i,j)'
%       om.init_indexed_name('lin', 'R', {2, 3});
%       for i = 1:2
%         for j = 1:3
%           om.add_lin_constraint('R', {i, j}, A{i,j}, ...);
%         end
%       end
%
% See also opt_model, add_var, add_lin_constraint, add_nln_constraint,
% add_quad_cost, add_nln_cost, add_legacy_cost.

%   MATPOWER
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% call parent method (also checks for valid type for named set)
om = init_indexed_name@opt_model(om, set_type, name, dim_list);

%% add type-specific info about this named set
switch set_type
    case 'cost'         %% cost set
        %% use column vector if single dimension
        if length(dim_list) == 1
            dim_list = {dim_list{:}, 1};
        end
        empty_cell  = cell(dim_list{:});
        om.cost.data.N.(name)   = empty_cell;
        om.cost.data.Cw.(name)  = empty_cell;
        om.cost.data.vs.(name)  = empty_cell;
end
