function om = add_named_set(om, set_type, name, idx, N, varargin)
% add_named_set - Adds a named set of variables/constraints/costs to the model.
% ::
%
%   -----  PRIVATE METHOD  -----
%
%   This method is intended to be a private method, used internally by
%   the public methods ADD_VAR, ADD_LIN_CONSTRAINT, ADD_NLN_CONSTRAINT
%   ADD_QUAD_COST, ADD_NLN_COST and ADD_LEGACY_COST.
%
%   Variable Set
%       OM.ADD_NAMED_SET('var', NAME, IDX_LIST, N, V0, VL, VU, VT);
%
%   Linear Constraint Set
%       OM.ADD_NAMED_SET('lin', NAME, IDX_LIST, N, A, L, U, VARSETS);
%
%   Nonlinear Inequality Constraint Set
%       OM.ADD_NAMED_SET('nle', NAME, IDX_LIST, N, FCN, HESS, COMPUTED_BY, VARSETS);
%
%   Nonlinear Inequality Constraint Set
%       OM.ADD_NAMED_SET('nli', NAME, IDX_LIST, N, FCN, HESS, COMPUTED_BY, VARSETS);
%
%   Quadratic Cost Set
%       OM.ADD_NAMED_SET('qdc', NAME, IDX_LIST, N, CP, VARSETS);
%
%   General Nonlinear Cost Set
%       OM.ADD_NAMED_SET('nlc', NAME, IDX_LIST, N, FCN, VARSETS);
%
%   Legacy Cost Set
%       OM.ADD_NAMED_SET('cost', NAME, IDX_LIST, N, CP, VARSETS);
%
% See also opt_model, add_var, add_lin_constraint, add_nln_constraint
% add_quad_cost, add_nln_cost, add_legacy_cost.

%   MATPOWER
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% call parent method (also checks for valid type for named set)
om = add_named_set@opt_model(om, set_type, name, idx, N, varargin{:});

switch set_type
    case 'cost'         %% cost set
        %% add type-specific data for named set
        om_ff = om.cost;
        om.cost = [];

        [cp, varsets] = deal(varargin{:});

        if isempty(idx)
            om_ff.data.N.(name)  = cp.N;
            om_ff.data.Cw.(name) = cp.Cw;
            om_ff.data.vs.(name) = varsets;
            if isfield(cp, 'H')
                om_ff.data.H.(name)  = cp.H;
            end
            if isfield(cp, 'dd')
                om_ff.data.dd.(name) = cp.dd;
            end
            if isfield(cp, 'rh')
                om_ff.data.rh.(name) = cp.rh;
            end
            if isfield(cp, 'kk')
                om_ff.data.kk.(name) = cp.kk;
            end
            if isfield(cp, 'mm')
                om_ff.data.mm.(name) = cp.mm;
            end
        else
            %% calls to substruct() are relatively expensive, so we pre-build the
            %% struct for addressing cell array fields
            %% sc = substruct('.', name, '{}', idx);
            sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field

            om_ff.data.N  = subsasgn(om_ff.data.N,  sc, cp.N);
            om_ff.data.Cw = subsasgn(om_ff.data.Cw, sc, cp.Cw);
            om_ff.data.vs = subsasgn(om_ff.data.vs, sc, varsets);
            if isfield(cp, 'H')
                om_ff.data.H = subsasgn(om_ff.data.H, sc, cp.H);
            end
            if isfield(cp, 'dd')
                om_ff.data.dd = subsasgn(om_ff.data.dd, sc, cp.dd);
            end
            if isfield(cp, 'rh')
                om_ff.data.rh = subsasgn(om_ff.data.rh, sc, cp.rh);
            end
            if isfield(cp, 'kk')
                om_ff.data.kk = subsasgn(om_ff.data.kk, sc, cp.kk);
            end
            if isfield(cp, 'mm')
                om_ff.data.mm = subsasgn(om_ff.data.mm, sc, cp.mm);
            end
        end
        if ~isempty(om_ff.params)       %% clear cache of aggregated params
            om_ff.params = [];
        end
        om.cost = om_ff;
end
