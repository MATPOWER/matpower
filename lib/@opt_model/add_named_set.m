function om = add_named_set(om, set_type, name, idx, N, varargin)
%ADD_NAMED_SET  Adds a named set of variables/constraints/costs to the model.
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
%   See also OPT_MODEL, ADD_VAR, ADD_LIN_CONSTRAINT, ADD_NLN_CONSTRAINT
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
    om_ff = om.(ff);
else
    error('@opt_model/add_named_set: ''%s'' is not a valid SET_TYPE, must be one of ''var'', ''lin'', ''nle'', ''nli'', ''qdc'', ''nlc'', ''cost''', set_type);
end

%% add general indexing info about this named set
if isempty(idx)     %% simple named set
    %% prevent duplicate name in set of specified type
    if isfield(om_ff.idx.N, name)
        error('@opt_model/add_named_set: %s set named ''%s'' already exists', st_label, name);
    end

    %% add indexing info about this set
    om_ff.idx.i1.(name)  = om_ff.N + 1; %% starting index
    om_ff.idx.iN.(name)  = om_ff.N + N; %% ending index
    om_ff.idx.N.(name)   = N;           %% number in set
    om_ff.N  = om_ff.idx.iN.(name);     %% number of elements of this type
    om_ff.NS = om_ff.NS + 1;            %% number of sets of this type
    om_ff.order(om_ff.NS).name = name;  %% add name to ordered list of sets
    om_ff.order(om_ff.NS).idx  = {};
else                %% indexed named set
    %% calls to substruct() are relatively expensive, so we pre-build the
    %% structs for addressing cell and numeric array fields
    %% sn = substruct('.', name, '()', idx);
    %% sc = substruct('.', name, '{}', idx);
    sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
    sn = sc; sn(2).type = '()';                             %% num array field

    %% prevent duplicate name in set of specified type
    if subsref(om_ff.idx.i1, sn) ~= 0
        str = '%d'; for m = 2:length(idx), str = [str ',%d']; end
        nname = sprintf(['%s(' str, ')'], name, idx{:});
        error('@opt_model/add_named_set: %s set named ''%s'' already exists', st_label, nname);
    end

    %% add indexing info about this set
    om_ff.idx.i1  = subsasgn(om_ff.idx.i1, sn, om_ff.N + 1);    %% starting index
    om_ff.idx.iN  = subsasgn(om_ff.idx.iN, sn, om_ff.N + N);    %% ending index
    om_ff.idx.N   = subsasgn(om_ff.idx.N,  sn, N);              %% number in set
    om_ff.N  = subsref(om_ff.idx.iN, sn);   %% number of elements of this type
    om_ff.NS = om_ff.NS + 1;                %% number of sets of this type
    om_ff.order(om_ff.NS).name = name;      %% add name to ordered list of sets
    om_ff.order(om_ff.NS).idx  = idx;       %% add indices to ordered list of sets
end

%% add type-specific data for named set
switch ff
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
    case 'lin'          %% linear constraint set
        [A, l, u, varsets] = deal(varargin{:});
        if isempty(idx)
            om_ff.data.A.(name)  = A;
            om_ff.data.l.(name)  = l;
            om_ff.data.u.(name)  = u;
            om_ff.data.vs.(name) = varsets;
        else
            om_ff.data.A  = subsasgn(om_ff.data.A, sc, A);
            om_ff.data.l  = subsasgn(om_ff.data.l, sc, l);
            om_ff.data.u  = subsasgn(om_ff.data.u, sc, u);
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
    case 'cost'         %% cost set
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
end
om.(ff) = om_ff;
