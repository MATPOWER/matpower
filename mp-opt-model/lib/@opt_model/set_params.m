function om = set_params(om, st, name, idx, params, vals)
% set_params - Modifies parameters for variable, cost or constraint in model
%
% .. note::
%    .. deprecated:: 4.3 Please use mp.sm_lin_constraint.set_params, 
%       mp.sm_nln_constraint.set_params, mp.sm_nln_cost.set_params, 
%       mp.sm_quad_cost.set_params, or mp.sm_variable.set_params instead,
%       as in ``om.lin.set_params(...)``, ``om.nle.set_params(...)``,
%       ``om.nli.set_params(...)``, ``om.nlc.set_params(...)``,
%       ``om.qdc.set_params(...)``, or ``om.var.set_params(...)``.
%
% ::
%
%   This method can be used to modify parameters for an existing variable,
%   constraint or cost in the model.
%
%   OM.SET_PARAMS(SET_TYPE, NAME, PARAMS, VALS)
%   OM.SET_PARAMS(SET_TYPE, NAME, IDX, PARAMS, VALS)
%
%   Inputs:
%       SET_TYPE : one of 'var', 'lin', 'nle', 'nli', 'nlc', 'qdc' for
%           variables, linear constraints, nonlinear equality constraints,
%           nonlinear inequality constraints, general nonlinear costs,
%           and quadratic costs, respectively
%       NAME : name of set
%       IDX : index of named set (for an indexed set)
%       PARAMS : can be one of three options:
%           1 - 'all', indicating that VALS is a cell array whose elements
%               correspond to the input parameters of the respective
%               add_*() method
%           2 - the name of a PARAM, VAL is the value of that parameter
%           3 - a cell array of PARAM names, VALS is a cell array of
%               corresponding values
%           Note: Changing the dimension of a 'var' is not allowed
%               and changing the #1 ('all') is the only option for 'nle', 'nli', and 'nlc'
%       VALS : new value or cell array of new values for PARAMS
%
%   Valid PARAM names:
%       var - N, v0, vl, vu, vt
%       lin - A, l, u, vs
%       nle - N, fcn, hess, include, vs
%       nli - N, fcn, hess, include, vs
%       nlc - N, fcn, vs
%       qdc - Q, c, k, vs
%
%   Examples:
%       om.set_params('var', 'Pg', 'v0', Pg0);
%       om.set_params('lin', 'y', {2,3}, {'l', 'u'}, {l, u});
%       om.set_params('nle', 'Pmis', 'all', {N, @fcn, @hess, vs});
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

if nargin < 6
    vals = params;
    params = idx;
    idx = {};
end

%% create default list of parameters to update based on set type & inputs
switch st
    case 'var'
        om.var.set_params(name, idx, params, vals);
    case {'lin', 'nle', 'nli', 'nlc', 'qdc'}
        om.(st).set_params(om.var, name, idx, params, vals);
    otherwise
        error('opt_model.set_params: ''%s'' is not a valid SET_TYPE', st);
end
