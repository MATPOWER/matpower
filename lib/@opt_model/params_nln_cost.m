function [N, fcn, vs] = params_nln_cost(om, name, idx)
%PARAMS_NLN_COST  Returns cost parameters for general nonlinear costs.
%   [N, FCN] = OM.PARAMS_NLN_COST(NAME)
%   [N, FCN] = OM.PARAMS_NLN_COST(NAME, IDX)
%   [N, FCN, VS] = OM.PARAMS_NLN_COST(...)
%
%   Returns the parameters N and FCN for the corresponding named general
%   nonlinear cost set. Likewise for indexed named sets specified
%   by NAME and IDX.
%
%   An optional 3rd output argument VS indicates the variable sets used by
%   this cost set.
%
%   See also OPT_MODEL, ADD_NLN_COST, EVAL_NLN_COST.

%   MATPOWER
%   Copyright (c) 2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

nlc = om.nlc;
if nargin < 3
    idx = {};
end
dims = size(nlc.idx.i1.(name));
if ~isempty(idx) || prod(dims) == 1 %% indexed, or simple named set
    %% calls to substruct() are relatively expensive, so we pre-build the
    %% structs for addressing cell and numeric array fields
    %% sn = substruct('.', name, '()', idx);
    %% sc = substruct('.', name, '{}', idx);
    sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
    sn = sc; sn(2).type = '()';                             %% num array field

    if isempty(idx)
        N = nlc.idx.N.(name);
        fcn = nlc.data.fcn.(name);
        vs = nlc.data.vs.(name);
    else
        N = subsref(nlc.idx.N, sn);
        fcn = subsref(nlc.data.fcn, sc);
        vs = subsref(nlc.data.vs, sc);
    end
else
    error('@opt_model/params_nln_cost: general nonlinear cost set ''%s'' requires an IDX arg', name)
end
