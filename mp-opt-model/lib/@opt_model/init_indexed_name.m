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
% add_quad_cost, add_nln_cost.

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% Due to a bug related to inheritance in constructors in
%% Octave 5.2 and earlier (https://savannah.gnu.org/bugs/?52614),
%% INIT_SET_TYPES() cannot be called directly in the
%% MP_IDX_MANAGER constructor, as desired.
%%
%% WORKAROUND:  Initialize MP_IDX_MANAGER fields here, if needed,
%%              after object construction, but before object use.
if isempty(om.var)          %% only if not already initialized
    om.init_set_types();
end

%% use column vector if single dimension
if length(dim_list) == 1
    dim_list = {dim_list{:}, 1};
end

%% call parent method (also checks for valid type for named set)
om = init_indexed_name@mp_idx_manager(om, set_type, name, dim_list);

%% add type-specific info about this named set
empty_cell  = cell(dim_list{:});
switch set_type
    case 'var'          %% variable set
        om.var.data.v0.(name)   = empty_cell;   %% initial value
        om.var.data.vl.(name)   = empty_cell;   %% lower bound
        om.var.data.vu.(name)   = empty_cell;   %% upper bound
        om.var.data.vt.(name)   = empty_cell;   %% variable type
    case 'lin'          %% linear constraint set
        om.lin.data.A.(name)    = empty_cell;
        om.lin.data.l.(name)    = empty_cell;
        om.lin.data.u.(name)    = empty_cell;
        om.lin.data.vs.(name)   = empty_cell;
    case 'nle'          %% nonlinear equality constraint set
        om.nle.data.fcn.(name)  = empty_cell;
        om.nle.data.hess.(name) = empty_cell;
        om.nle.data.vs.(name)   = empty_cell;
    case 'nli'          %% nonlinear inequality constraint set
        om.nli.data.fcn.(name)  = empty_cell;
        om.nli.data.hess.(name) = empty_cell;
        om.nli.data.vs.(name)   = empty_cell;
    case 'nlc'          %% nonlinear cost set
        om.nlc.data.fcn.(name)  = empty_cell;
        om.nlc.data.vs.(name)   = empty_cell;
    case 'qdc'          %% quadratic cost set
        om.qdc.data.Q.(name)    = empty_cell;
        om.qdc.data.c.(name)    = empty_cell;
        om.qdc.data.k.(name)    = empty_cell;
        om.qdc.data.vs.(name)   = empty_cell;
end
