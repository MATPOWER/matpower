function om = add_legacy_cost(om, name, idx, varargin)
%ADD_LEGACY_COST  Adds a set of user costs to the model.
%
%   -----  DEPRECATED - Please use one of the following instead:        -----
%   -----  ADD_QUAD_COST, ADD_NLN_COST, INIT_INDEXED_NAME  -----
%
%   OM.ADD_LEGACY_COST(NAME, CP);
%   OM.ADD_LEGACY_COST(NAME, CP, VARSETS);
%   OM.ADD_LEGACY_COST(NAME, IDX_LIST, CP);
%   OM.ADD_LEGACY_COST(NAME, IDX_LIST, CP, VARSETS);
%
%   Adds a named block of user-defined costs to the model. Each set is
%   defined by the CP struct described below. All user-defined sets of
%   costs are combined together into a single set of cost parameters in
%   a single CP struct by BULD_COST_PARAMS. This full aggregate set of
%   cost parameters can be retreived from the model by GET_COST_PARAMS.
%
%   Examples:
%       cp1 = struct('N', N1, 'Cw', Cw1);
%       cp2 = struct('N', N2, 'Cw', Cw2, 'H', H, 'dd', dd, ...
%                     'rh', rh, 'kk', kk, 'mm', mm);
%       om.add_legacy_cost('usr1', cp1, {'Pg', 'Qg', 'z'});
%       om.add_legacy_cost('usr2', cp2, {'Vm', 'Pg', 'Qg', 'z'});
%
%       om.init_indexed_name('c', {2, 3});
%       for i = 1:2
%         for j = 1:3
%           om.add_legacy_cost('c', {i, j}, cp(i,j), ...);
%         end
%       end
%
%   Let x refer to the vector formed by combining the specified VARSETS,
%   and f_u(x, CP) be the cost at x corresponding to the cost parameters
%   contained in CP, where CP is a struct with the following fields:
%       N      - nw x nx sparse matrix (optional, identity matrix by default)
%       Cw     - nw x 1 vector
%       H      - nw x nw sparse matrix (optional, all zeros by default)
%       dd, mm - nw x 1 vectors (optional, all ones by default)
%       rh, kk - nw x 1 vectors (optional, all zeros by default)
%
%   These parameters are used as follows to compute f_u(x, CP)
%
%       R  = N*x - rh
%
%               /  kk(i),  R(i) < -kk(i)
%       K(i) = <   0,     -kk(i) <= R(i) <= kk(i)
%               \ -kk(i),  R(i) > kk(i)
%
%       RR = R + K
%
%       U(i) =  /  0, -kk(i) <= R(i) <= kk(i)
%               \  1, otherwise
%
%       DDL(i) = /  1, dd(i) = 1
%                \  0, otherwise
%
%       DDQ(i) = /  1, dd(i) = 2
%                \  0, otherwise
%
%       Dl = diag(mm) * diag(U) * diag(DDL)
%       Dq = diag(mm) * diag(U) * diag(DDQ)
%
%       w = (Dl + Dq * diag(RR)) * RR
%
%       f_u(x, CP) = 1/2 * w'*H*w + Cw'*w
%
%   See also OPT_MODEL, PARAMS_LEGACY_COST, EVAL_LEGACY_COST.

%   MATPOWER
%   Copyright (c) 2008-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if iscell(idx)
    cp = varargin{1};
    args = varargin(2:end);
else                            %% simple named set
    cp = idx;
    args = varargin;
    idx = {};
end

if isempty(args)
    varsets = {};
else
    varsets = args{1};
end

%% convert varsets from cell to struct array if necessary
varsets = om.varsets_cell2struct(varsets);
nv = om.varsets_len(varsets);   %% number of variables

if isfield(cp, 'N')
    [nw, nx] = size(cp.N);
else
    nw = length(cp.Cw);
    nx = nw;
    cp.N = speye(nw, nx);
end

%% check sizes
if nx ~= nv
    if nw == 0
        cp.N = sparse(nw, nx);
    else
        error('@opt_model/add_legacy_cost: number of columns in N (%d x %d) does not match\nnumber of variables (%d)\n', nw, nx, nv);
    end
end
if size(cp.Cw, 1) ~= nw
    error('@opt_model/add_legacy_cost: number of rows of Cw (%d x %d) and N (%d x %d) must match\n', size(cp.Cw), nw, nx);
end
if isfield(cp, 'H') && (size(cp.H, 1) ~= nw || size(cp.H, 2) ~= nw)
    error('@opt_model/add_legacy_cost: both dimensions of H (%d x %d) must match the number of rows in N (%d x %d)\n', size(cp.H), nw, nx);
end
if isfield(cp, 'dd') && size(cp.dd, 1) ~= nw
    error('@opt_model/add_legacy_cost: number of rows of dd (%d x %d) and N (%d x %d) must match\n', size(cp.dd), nw, nx);
end
if isfield(cp, 'rh') && size(cp.rh, 1) ~= nw
    error('@opt_model/add_legacy_cost: number of rows of rh (%d x %d) and N (%d x %d) must match\n', size(cp.rh), nw, nx);
end
if isfield(cp, 'kk') && size(cp.kk, 1) ~= nw
    error('@opt_model/add_legacy_cost: number of rows of kk (%d x %d) and N (%d x %d) must match\n', size(cp.kk), nw, nx);
end
if isfield(cp, 'mm') && size(cp.mm, 1) ~= nw
    error('@opt_model/add_legacy_cost: number of rows of mm (%d x %d) and N (%d x %d) must match\n', size(cp.mm), nw, nx);
end

%% add the legacy cost set
om.add_named_set('cost', name, idx, nw, cp, varsets);
