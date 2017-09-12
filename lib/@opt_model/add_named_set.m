function om = add_named_set(om, set_type, name, idx, N, varargin)
%ADD_NAMED_SET  Adds a named set of variables/constraints/costs to the model.
%
%   -----  PRIVATE METHOD  -----
%
%   This method is intended to be a private method, used internally by
%   the public methods ADD_VARS, ADD_LIN_CONSTRAINTS, ADD_NLN_CONSTRAINT
%   ADD_QUAD_COST, ADD_NLN_COST and ADD_LEGACY_COSTS.
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
%   See also OPT_MODEL, ADD_VARS, ADD_LIN_CONSTRAINTS, ADD_NLN_CONSTRAINT
%            ADD_QUAD_COST and ADD_NLN_COST.

%   MATPOWER
%   Copyright (c) 2008-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% check for valid type for named set
st_label = om.valid_named_set_type(set_type);
if st_label
    ff = set_type;
else
    error('@opt_model/add_named_set: ''%s'' is not a valid SET_TYPE, must be one of ''var'', ''lin'', ''nle'', ''nli'', ''cost''', set_type);
end

%% add general indexing info about this named set
if isempty(idx)     %% simple named set
    %% prevent duplicate name in set of specified type
    if isfield(om.(ff).idx.N, name)
        error('@opt_model/add_named_set: %s set named ''%s'' already exists', st_label, name);
    end

    %% add indexing info about this set
    om.(ff).idx.i1.(name)  = om.(ff).N + 1; %% starting index
    om.(ff).idx.iN.(name)  = om.(ff).N + N; %% ending index
    om.(ff).idx.N.(name)   = N;             %% number in set
    om.(ff).N  = om.(ff).idx.iN.(name);     %% number of elements of this type
    om.(ff).NS = om.(ff).NS + 1;            %% number of sets of this type
    om.(ff).order(om.(ff).NS).name = name;  %% add name to ordered list of sets
    om.(ff).order(om.(ff).NS).idx  = {};
else                %% indexed named set
    %% calls to substruct() are relatively expensive, so we pre-build the
    %% structs for addressing cell and numeric array fields
    %% sn = substruct('.', name, '()', idx);
    %% sc = substruct('.', name, '{}', idx);
    sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
    sn = sc; sn(2).type = '()';                             %% num array field

    %% prevent duplicate name in set of specified type
    if subsref(om.(ff).idx.i1, sn) ~= 0
        str = '%d'; for m = 2:length(idx), str = [str ',%d']; end
        nname = sprintf(['%s(' str, ')'], name, idx{:});
        error('@opt_model/add_named_set: %s set named ''%s'' already exists', st_label, nname);
    end

    %% add indexing info about this set
    om.(ff).idx.i1  = subsasgn(om.(ff).idx.i1, sn, om.(ff).N + 1);  %% starting index
    om.(ff).idx.iN  = subsasgn(om.(ff).idx.iN, sn, om.(ff).N + N);  %% ending index
    om.(ff).idx.N   = subsasgn(om.(ff).idx.N,  sn, N);              %% number in set
    om.(ff).N  = subsref(om.(ff).idx.iN, sn);   %% number of elements of this type
    om.(ff).NS = om.(ff).NS + 1;                %% number of sets of this type
    om.(ff).order(om.(ff).NS).name = name;      %% add name to ordered list of sets
    om.(ff).order(om.(ff).NS).idx  = idx;       %% add indexes to ordered list of sets
end

%% add type-specific data for named set
switch ff
    case 'var'          %% variable set
        [v0, vl, vu, vt] = deal(varargin{:});
        if isempty(idx)
            om.var.data.v0.(name) = v0;             %% initial value
            om.var.data.vl.(name) = vl;             %% lower bound
            om.var.data.vu.(name) = vu;             %% upper bound
            om.var.data.vt.(name) = vt;             %% variable type
        else
            om.var.data.v0 = subsasgn(om.var.data.v0, sc, v0); %% initial value
            om.var.data.vl = subsasgn(om.var.data.vl, sc, vl); %% lower bound
            om.var.data.vu = subsasgn(om.var.data.vu, sc, vu); %% upper bound
            om.var.data.vt = subsasgn(om.var.data.vt, sc, vt); %% variable type
        end
    case 'lin'          %% linear constraint set
        [A, l, u, varsets] = deal(varargin{:});
        if isempty(idx)
            om.lin.data.A.(name)  = A;
            om.lin.data.l.(name)  = l;
            om.lin.data.u.(name)  = u;
            om.lin.data.vs.(name) = varsets;
        else
            om.lin.data.A  = subsasgn(om.lin.data.A, sc, A);
            om.lin.data.l  = subsasgn(om.lin.data.l, sc, l);
            om.lin.data.u  = subsasgn(om.lin.data.u, sc, u);
            om.lin.data.vs = subsasgn(om.lin.data.vs, sc, varsets);
        end
        if ~isempty(om.lin.params)      %% clear cache of aggregated params
            om.lin.params = [];
        end
    case {'nle', 'nli'} %% nonlinear constraint set
        [fcn, hess, computed_by, varsets] = deal(varargin{:});
        if isempty(idx)
            if isempty(computed_by)
                if ~isempty(fcn)
                    om.(ff).data.fcn.(name)  = fcn;
                    om.(ff).data.hess.(name) = hess;
                end
            else
                if isfield(om.(ff).data.include, computed_by)
                    om.(ff).data.include.(computed_by).name{end+1} = name;
                    om.(ff).data.include.(computed_by).N(end+1) = N;
                else
                    om.(ff).data.include.(computed_by).name = {name};
                    om.(ff).data.include.(computed_by).N = N;
                end
            end
            om.(ff).data.vs.(name) = varsets;
        else
            if ~isempty(computed_by)
                error('add_named_set: a nonlinear constraint set computed by another set is currently only implemented for simple named sets, not yet for indexed named sets');
            end
            om.(ff).data.fcn  = subsasgn(om.(ff).data.fcn, sc, fcn);
            om.(ff).data.hess = subsasgn(om.(ff).data.hess, sc, hess);
            om.(ff).data.vs   = subsasgn(om.(ff).data.vs, sc, varsets);
        end
    case 'qdc'          %% quadratic cost set
        [Q, c, k, varsets] = deal(varargin{:});
        if isempty(idx)
            om.qdc.data.Q.(name)  = Q;
            om.qdc.data.c.(name)  = c;
            om.qdc.data.k.(name)  = k;
            om.qdc.data.vs.(name) = varsets;
        else
            om.qdc.data.Q  = subsasgn(om.qdc.data.Q, sc, Q);
            om.qdc.data.c  = subsasgn(om.qdc.data.c, sc, c);
            om.qdc.data.k  = subsasgn(om.qdc.data.k, sc, k);
            om.qdc.data.vs = subsasgn(om.qdc.data.vs, sc, varsets);
        end
        if ~isempty(om.qdc.params)      %% clear cache of aggregated params
            om.qdc.params = [];
        end
    case 'nlc'          %% general nonlinear cost set
        [fcn, varsets] = deal(varargin{:});
        if isempty(idx)
            om.nlc.data.fcn.(name)  = fcn;
            om.nlc.data.vs.(name) = varsets;
        else
            om.nlc.data.fcn  = subsasgn(om.nlc.data.fcn, sc, fcn);
            om.nlc.data.vs   = subsasgn(om.nlc.data.vs, sc, varsets);
        end
    case 'cost'         %% cost set
        [cp, varsets] = deal(varargin{:});
        if isempty(idx)
            om.cost.data.N.(name)  = cp.N;
            om.cost.data.Cw.(name) = cp.Cw;
            om.cost.data.vs.(name) = varsets;
            if isfield(cp, 'H')
                om.cost.data.H.(name)  = cp.H;
            end
            if isfield(cp, 'dd')
                om.cost.data.dd.(name) = cp.dd;
            end
            if isfield(cp, 'rh')
                om.cost.data.rh.(name) = cp.rh;
            end
            if isfield(cp, 'kk')
                om.cost.data.kk.(name) = cp.kk;
            end
            if isfield(cp, 'mm')
                om.cost.data.mm.(name) = cp.mm;
            end
        else
            om.cost.data.N  = subsasgn(om.cost.data.N,  sc, cp.N);
            om.cost.data.Cw = subsasgn(om.cost.data.Cw, sc, cp.Cw);
            om.cost.data.vs = subsasgn(om.cost.data.vs, sc, varsets);
            if isfield(cp, 'H')
                om.cost.data.H = subsasgn(om.cost.data.H, sc, cp.H);
            end
            if isfield(cp, 'dd')
                om.cost.data.dd = subsasgn(om.cost.data.dd, sc, cp.dd);
            end
            if isfield(cp, 'rh')
                om.cost.data.rh = subsasgn(om.cost.data.rh, sc, cp.rh);
            end
            if isfield(cp, 'kk')
                om.cost.data.kk = subsasgn(om.cost.data.kk, sc, cp.kk);
            end
            if isfield(cp, 'mm')
                om.cost.data.mm = subsasgn(om.cost.data.mm, sc, cp.mm);
            end
        end
        if ~isempty(om.cost.params)     %% clear cache of aggregated params
            om.cost.params = [];
        end
end
