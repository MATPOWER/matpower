function [N, fcn, vs] = params_nln_cost(om, name, idx)
% params_nln_cost - Returns cost parameters for general nonlinear costs.
% ::
%
%   [N, FCN] = OM.PARAMS_NLN_COST(NAME)
%   [N, FCN] = OM.PARAMS_NLN_COST(NAME, IDX_LIST)
%   [N, FCN, VS] = OM.PARAMS_NLN_COST(...)
%
%   Returns the parameters N and FCN provided when the corresponding
%   named general nonlinear cost set was added to the model. Likewise
%   for indexed named sets specified by NAME and IDX_LIST.
%
%   An optional 3rd output argument VS indicates the variable sets used by
%   this cost set.
%
% See also opt_model, add_nln_cost, eval_nln_cost.

%   MP-Opt-Model
%   Copyright (c) 2017-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

nlc = om.nlc;
if nargin < 3
    idx = {};
end

if isempty(idx)
    dims = size(nlc.idx.i1.(name));
    if prod(dims) ~= 1
        error('opt_model.params_nln_cost: general nonlinear cost set ''%s'' requires an IDX_LIST arg', name)
    end
    N = nlc.idx.N.(name);
    fcn = nlc.data.fcn.(name);
    vs = nlc.data.vs.(name);
else
    %% calls to substruct() are relatively expensive, so we pre-build the
    %% structs for addressing cell and numeric array fields
    %% sn = substruct('.', name, '()', idx);
    %% sc = substruct('.', name, '{}', idx);
    sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
    sn = sc; sn(2).type = '()';                             %% num array field

    N = subsref(nlc.idx.N, sn);
    fcn = subsref(nlc.data.fcn, sc);
    vs = subsref(nlc.data.vs, sc);
end
