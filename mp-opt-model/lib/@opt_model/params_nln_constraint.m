function [N, fcn, hess, vs, include] = params_nln_constraint(om, iseq, name, idx)
% params_nln_constraint - Returns parameters for general nonlinear constraints.
% ::
%
%   N = OM.PARAMS_NLN_CONSTRAINT(ISEQ, NAME)
%   N = OM.PARAMS_NLN_CONSTRAINT(ISEQ, NAME, IDX_LIST)
%   [N, FCN] = OM.PARAMS_NLN_CONSTRAINT(...)
%   [N, FCN, HESS] = OM.PARAMS_NLN_CONSTRAINT(...)
%   [N, FCN, HESS, VS] = OM.PARAMS_NLN_CONSTRAINT(...)
%   [N, FCN, HESS, VS, INCLUDE] = OM.PARAMS_NLN_CONSTRAINT(...)
%
%   Returns the parameters N, and optionally FCN, and HESS provided when
%   the corresponding named nonlinear constraint set was added to the
%   model. Likewise for indexed named sets specified by NAME and IDX_LIST.
%   The ISEQ input should be set to 1 for equality constrainst and to 0
%   for inequality constraints.
%
%   An optional 4th output argument VS indicates the variable sets used by
%   this constraint set.

%   And, for constraint sets whose functions compute the constraints for
%   another set, an optional 5th output argument returns a struct with a
%   cell array of set names in the 'name' field and an array of
%   corresponding dimensions in the 'N' field.
%
% See also opt_model, add_nln_constraint, eval_nln_constraint.

%   MP-Opt-Model
%   Copyright (c) 2017-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% get constraint type
if iseq         %% equality constraints
    om_nlx = om.nle;
else            %% inequality constraints
    om_nlx = om.nli;
end

if nargin < 4
    idx = {};
end

if isempty(idx)
    dims = size(om_nlx.idx.i1.(name));
    if prod(dims) ~= 1
        error('opt_model.params_nln_constraint: nonlinear constraint set ''%s'' requires an IDX_LIST arg', name)
    end
    N = om_nlx.idx.N.(name);
    if nargout > 1
        if isfield(om_nlx.data.fcn, name)
            fcn = om_nlx.data.fcn.(name);
        else
            fcn = '';
        end
        if nargout > 2
            if isfield(om_nlx.data.hess, name)
                hess = om_nlx.data.hess.(name);
            else
                hess = '';
            end
            if nargout > 3
                if isfield(om_nlx.data.vs, name)
                    vs = om_nlx.data.vs.(name);
                else
                    vs = {};
                end
                if nargout > 4
                    if isfield(om_nlx.data.include, name)
                        include = om_nlx.data.include.(name);
                    else
                        include = '';
                    end
                end
            end
        end
    end
else
    %% calls to substruct() are relatively expensive, so we pre-build the
    %% structs for addressing cell and numeric array fields
    %% sn = substruct('.', name, '()', idx);
    %% sc = substruct('.', name, '{}', idx);
    sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
    sn = sc; sn(2).type = '()';                             %% num array field

    N = subsref(om_nlx.idx.N, sn);
    if nargout > 1
        if isfield(om_nlx.data.fcn, name)
            fcn = subsref(om_nlx.data.fcn, sc);
        else
            fcn = '';
        end
        if nargout > 2
            if isfield(om_nlx.data.hess, name)
                hess = subsref(om_nlx.data.hess, sc);
            else
                hess = '';
            end
            if nargout > 3
                if isfield(om_nlx.data.vs, name)
                    vs = subsref(om_nlx.data.vs, sc);
                else
                    vs = {};
                end
                if nargout > 4
                    error('opt_model.params_nln_constraint: nonlinear constraint set ''%s'' cannot return INCLUDE, since a nonlinear constraint set computed by another set is currently only implemented for simple named sets, not yet for indexed named sets', name)
                end
            end
        end
    end
end
