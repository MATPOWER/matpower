function om = add_named_set(om, set_type, name, idx, N, varargin)
% add_named_set - Adds a named set of variables/constraints/costs to the model.
% ::
%
%   -----  PRIVATE METHOD  -----
%
%   This method is intended to be a private method, used internally by
%   the public methods ADD_VAR, ADD_LIN_CONSTRAINT, ADD_NLN_CONSTRAINT
%   ADD_QUAD_COST and ADD_NLN_COST.
%
%   Variable Set
%       OM.ADD_NAMED_SET('var', NAME, IDX_LIST, N, V0, VL, VU, VT);
%
%   Linear Constraint Set
%       OM.ADD_NAMED_SET('lin', NAME, IDX_LIST, N, A, L, U, VARSETS, TR);
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
% See also opt_model, add_var, add_lin_constraint, add_nln_constraint
% add_quad_cost, add_nln_cost.

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% call parent method (also checks for valid type for named set)
om = add_named_set@mp_idx_manager(om, set_type, name, idx, N, varargin{:});

if ~isempty(idx)
    %% calls to substruct() are relatively expensive, so we pre-build the
    %% struct for addressing cell array fields
    %% sc = substruct('.', name, '{}', idx);
    sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
end

%% add type-specific data for named set
om_ff = om.(set_type);
om.(set_type) = [];
switch set_type
    case 'var'          %% variable set
        [v0, vl, vu, vt] = deal(varargin{:});
        if isempty(idx)
            om_ff.data.v0.(name) = v0;          %% initial value
            om_ff.data.vl.(name) = vl;          %% lower bound
            om_ff.data.vu.(name) = vu;          %% upper bound
            om_ff.data.vt.(name) = vt;          %% variable type
        else
            om_ff.data.v0 = subsasgn(om_ff.data.v0, sc, v0);    %% initial value
            om_ff.data.vl = subsasgn(om_ff.data.vl, sc, vl);    %% lower bound
            om_ff.data.vu = subsasgn(om_ff.data.vu, sc, vu);    %% upper bound
            om_ff.data.vt = subsasgn(om_ff.data.vt, sc, vt);    %% variable type
        end
        if ~isempty(om_ff.params)       %% clear cache of aggregated params
            om_ff.params = [];
        end
    case 'lin'          %% linear constraint set
        [A, l, u, varsets, tr] = deal(varargin{:});
        if isempty(idx)
            om_ff.data.A.(name)  = A;
            om_ff.data.l.(name)  = l;
            om_ff.data.u.(name)  = u;
            om_ff.data.tr.(name)  = tr;
            om_ff.data.vs.(name) = varsets;
        else
            om_ff.data.A  = subsasgn(om_ff.data.A, sc, A);
            om_ff.data.l  = subsasgn(om_ff.data.l, sc, l);
            om_ff.data.u  = subsasgn(om_ff.data.u, sc, u);
            om_ff.data.tr  = subsasgn(om_ff.data.tr, sc, tr);
            om_ff.data.vs = subsasgn(om_ff.data.vs, sc, varsets);
        end
        if ~isempty(om_ff.params)       %% clear cache of aggregated params
            om_ff.params = [];
        end
    case {'nle', 'nli'} %% nonlinear constraint set
        [fcn, hess, computed_by, varsets] = deal(varargin{:});
        if isempty(idx)
            if isempty(computed_by)
                if ~isempty(fcn)
                    om_ff.data.fcn.(name)  = fcn;
                    om_ff.data.hess.(name) = hess;
                end
            else
                if isfield(om_ff.data.include, computed_by)
                    om_ff.data.include.(computed_by).name{end+1} = name;
                    om_ff.data.include.(computed_by).N(end+1) = N;
                else
                    om_ff.data.include.(computed_by).name = {name};
                    om_ff.data.include.(computed_by).N = N;
                end
            end
            om_ff.data.vs.(name) = varsets;
        else
            if ~isempty(computed_by)
                error('add_named_set: a nonlinear constraint set computed by another set is currently only implemented for simple named sets, not yet for indexed named sets');
            end
            om_ff.data.fcn  = subsasgn(om_ff.data.fcn, sc, fcn);
            om_ff.data.hess = subsasgn(om_ff.data.hess, sc, hess);
            om_ff.data.vs   = subsasgn(om_ff.data.vs, sc, varsets);
        end
    case 'qdc'          %% quadratic cost set
        [Q, c, k, varsets] = deal(varargin{:});
        if isempty(idx)
            om_ff.data.Q.(name)  = Q;
            om_ff.data.c.(name)  = c;
            om_ff.data.k.(name)  = k;
            om_ff.data.vs.(name) = varsets;
        else
            om_ff.data.Q  = subsasgn(om_ff.data.Q, sc, Q);
            om_ff.data.c  = subsasgn(om_ff.data.c, sc, c);
            om_ff.data.k  = subsasgn(om_ff.data.k, sc, k);
            om_ff.data.vs = subsasgn(om_ff.data.vs, sc, varsets);
        end
        if ~isempty(om_ff.params)       %% clear cache of aggregated params
            om_ff.params = [];
        end
    case 'nlc'          %% general nonlinear cost set
        [fcn, varsets] = deal(varargin{:});
        if isempty(idx)
            om_ff.data.fcn.(name)  = fcn;
            om_ff.data.vs.(name) = varsets;
        else
            om_ff.data.fcn  = subsasgn(om_ff.data.fcn, sc, fcn);
            om_ff.data.vs   = subsasgn(om_ff.data.vs, sc, varsets);
        end
end
om.(set_type) = om_ff;
