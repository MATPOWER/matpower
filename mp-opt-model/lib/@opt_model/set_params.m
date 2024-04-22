function om = set_params(om, st, name, idx, params, vals)
% set_params - Modifies parameters for variable, cost or constraint in model
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
        default_params = {'N', 'v0', 'vl', 'vu', 'vt'};
    case 'lin'
        default_params = {'A', 'l', 'u', 'vs', 'tr'};
    case {'nle', 'nli'}
        default_params = {'N', 'fcn', 'hess', 'vs'};
    case 'nlc'
        default_params = {'N', 'fcn', 'vs'};
    case 'qdc'
        default_params = {'Q', 'c', 'k', 'vs'};
    otherwise
        error('opt_model.set_params: ''%s'' is not a valid SET_TYPE', st);
end

%% standardize provided arguments in cell arrays params, vals
is_all = 0;     %% flag to indicate all params for set are being replaced
if ischar(params)
    if strcmp(params, 'all')
        is_all = 1;
        np = length(vals);      %% number of parameter values provided
        params = default_params(1:np);
    else
        np = 1;                 %% number of parameter values provided
        params = {params};
        vals = {vals};
    end
else
    np = length(vals);          %% number of parameter values provided
end

if ~isempty(idx)
    %% calls to substruct() are relatively expensive, so we pre-build the
    %% struct for addressing cell array fields
    %% sc = substruct('.', name, '{}', idx);
    sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
end

switch st
    case 'var'
        %% get current parameters
        [v0, vl, vu, vt] = om.params_var(name, idx);
        N0 = om.getN('var', name, idx);
        p = struct('N', N0, 'v0', v0, 'vl', vl, 'vu', vu, 'vt', vt);    %% current parameters
        u = struct('N',  0, 'v0',  0, 'vl',  0, 'vu',  0, 'vt',  0);    %% which ones to update

        %% replace with new parameters
        for k = 1:np
            p.(params{k}) = vals{k};
            u.(params{k}) = 1;
        end
        N = p.N;

        %% set missing default params for 'all'
        if is_all
            if np < 5
                p.vt = 'C';
                u.vt = 1;               %% update vt
                if np < 4
                    p.vu = Inf(N, 1);
                    u.vu = 1;           %% update vu
                    if np < 3
                        p.vl = -Inf(N, 1);
                        u.vl = 1;       %% update vl
                        if np < 2
                            p.v0 = zeros(N, 1);
                            u.v0 = 1;   %% update v0
                        end
                    end
                end
            end
        end

        %% check consistency of parameters
        %% no dimension change
        if N ~= N0
            error('opt_model.set_params: dimension change for ''%s'' ''%s'' not allowed', st, nameidxstr(name, idx));
        end

        %% check sizes of new values of v0, vl, vu, vt
        for pn = {'v0', 'vl', 'vu', 'vt'}
            if u.(pn{1})
                nn = length(p.(pn{1}));
                if nn ~= N
                    if nn == 0
                        switch pn{1}
                            case 'v0'
                                p.(pn{1}) = zeros(N, 0);
                            case 'vl'
                                p.(pn{1}) = -Inf(N, 0);
                            case 'vu'
                                p.(pn{1}) =  Inf(N, 0);
                            case 'vt'
                                p.(pn{1}) = 'C';
                        end
                    elseif nn == 1
                        if pn{1} ~= 'vt'
                            p.(pn{1}) = p.(pn{1}) * ones(N, 1);   %% expand from scalar
                        end
                    else
                        error('opt_model.set_params: parameter ''%s'' ''%s'' ''%s'' should have length %d (or 1)', st, nameidxstr(name, idx), pn{1}, N);
                    end
                end
            end
        end

        %% assign new parameters
        if isempty(idx)     %% simple named set
            for k = 2:length(default_params)
                pn = default_params{k};     %% param name
                if u.(pn)   %% assign new val for this parameter
                    om.var.data.(pn).(name) = p.(pn);
                end
            end
        else                %% indexed named set
            for k = 2:length(default_params)
                pn = default_params{k};     %% param name
                if u.(pn)   %% assign new val for this parameter
                    om.var.data.(pn) = subsasgn(om.var.data.(pn), sc, p.(pn));
                end
            end
        end

        %% clear cached parameters
        om.var.params = [];
    case 'lin'
        %% get current parameters
        [A, l, u, vs, ~, ~, tr] = om.params_lin_constraint(name, idx);
        if tr
            [M0, N0] = size(A);
        else
            [N0, M0] = size(A);
        end
        if isempty(vs), vs = {vs}; end
        p = struct('A', A, 'l', l, 'u', u, 'vs', vs, 'tr', tr); %% current parameters
        u = struct('A', 0, 'l', 0, 'u', 0, 'vs',  0, 'tr', 0);  %% which ones to update

        %% replace with new parameters
        for k = 1:np
            p.(params{k}) = vals{k};
            u.(params{k}) = 1;
        end

        %% set missing default params for 'all'
        if p.tr
            [M, N] = size(p.A);
        else
            [N, M] = size(p.A);
        end
        if is_all
            u.A = 1;            %% always update A
            u.l = 1;            %% alwaus update l
            if np < 5
                if p.tr
                    p.tr = 0;
                    u.tr = 1;       %% update tr
                    [N, M] = deal(M, N);
                end
                if np < 4
                    p.vs = {};
                    u.vs = 1;       %% update vs
                    if np < 3
                        p.u = Inf(N, 1);
                        u.u = 1;    %% update u
                    end
                end
            end
        end

        %% check consistency of parameters
        %% no dimension change unless 'all'
        if N ~= N0 && ~is_all
            error('opt_model.set_params: dimension change for ''%s'' ''%s'' not allowed except for ''all''', st, nameidxstr(name, idx));
        end
        %% no transpose change unless providing A
        if u.tr && ~u.A
            error('opt_model.set_params: update to ''tr'' for ''%s'' ''%s'' requires update to ''A''', st, nameidxstr(name, idx));
        end

        %% check sizes of new values of l, u
        for pn = {'l', 'u'}
            if u.(pn{1})
                nn = length(p.(pn{1}));
                if nn ~= N
                    if nn == 0
                        switch pn{1}
                            case 'l'
                                p.(pn{1}) = -Inf(N, 0);
                            case 'u'
                                p.(pn{1}) =  Inf(N, 0);
                        end
                    elseif nn == 1
                        p.(pn{1}) = p.(pn{1}) * ones(N, 1);   %% expand from scalar
                    else
                        error('opt_model.set_params: parameter ''%s'' ''%s'' ''%s'' should have length %d (or 1)', st, nameidxstr(name, idx), pn{1}, N);
                    end
                end
            end
        end

        %% check consistency of A and vs
        if u.A || u.vs
            p.vs = om.varsets_cell2struct(p.vs);
            nv = om.varsets_len(p.vs);      %% number of variables
            if M ~= nv
                error('opt_model.set_params: for ''%s'' ''%s'' number of columns of ''A'' (%d) must be consistent with ''vs'' (%d)', st, nameidxstr(name, idx), M, nv);
            end
        end

        %% assign new parameters
        if isempty(idx)     %% simple named set
            for k = 1:length(default_params)
                pn = default_params{k};     %% param name
                if u.(pn)   %% assign new val for this parameter
                    om.lin.data.(pn).(name) = p.(pn);
                end
            end
        else                %% indexed named set
            for k = 1:length(default_params)
                pn = default_params{k};     %% param name
                if u.(pn)   %% assign new val for this parameter
                    om.lin.data.(pn) = subsasgn(om.lin.data.(pn), sc, p.(pn));
                end
            end
        end

        %% clear cached parameters
        om.lin.params = [];
    case 'qdc'
        %% get current parameters
        [Q, c, kk, vs] = om.params_quad_cost(name, idx);
        [MQ0, NQ0] = size(Q);
        [Mc0, Nc0] = size(c);
        Nk0 = length(kk);
        nx0 = max([MQ0 Mc0 Nk0]);
        N0 = om.getN('qdc', name, idx);
        if isempty(vs), vs = {vs}; end
        p = struct('Q', Q, 'c', c, 'k', kk, 'vs', vs);  %% current parameters
        u = struct('Q', 0, 'c', 0, 'k',  0, 'vs',  0);  %% which ones to update

        %% replace with new parameters
        for k = 1:np
            p.(params{k}) = vals{k};
            u.(params{k}) = 1;
        end

        %% set missing default params for 'all'
        [MQ, NQ] = size(p.Q);
        [Mc, Nc] = size(p.c);
        Nk = length(p.k);
        nx = max([MQ Mc Nk]);
        if NQ == 1 || (isempty(p.Q) && (Nk > 1 || k == 0))
            %% Q is a column vector (cost is element-wise, i.e. a vector)
            %% OR Q is empty and k is either a vector or zero
            N = nx;
        else            %% Q is a square matrix (cost is a scalar)
            N = 1;
        end
        if is_all
            u.Q = 1;                %% always update Q
            if np < 4
                p.vs = {};
                u.vs = 1;           %% update vs
                if np < 3
                    p.k = 0;
                    u.k = 1;        %% update k
                    if np < 2
                        p.c = [];
                        u.c = 1;    %% update c
                    end
                end
            end
        end

        %% check consistency of parameters
        %% no dimension change unless 'all'
        if (N ~= N0 || nx ~= nx0) && ~is_all
            error('opt_model.set_params: dimension change for ''%s'' ''%s'' not allowed except for ''all''', st, nameidxstr(name, idx));
        end

        %% Q and c can't both be empty
        if ~MQ && ~Mc
            error('opt_model.set_params: ''%s'' ''%s'' : ''Q'' and ''c'' cannot both be empty', st, nameidxstr(name, idx));
        end

        %% check sizes of new values of Q, c, k
        if ~isempty(p.Q) && (MQ ~= nx || (MQ ~= NQ && NQ ~= 1) )
            error('opt_model.set_params: ''%s'' ''%s'' : ''%s'' is expected to be (%d x %d)', st, nameidxstr(name, idx), 'Q', MQ, NQ);
        end
        if ~isempty(p.c) && Mc ~= nx
            error('opt_model.set_params: ''%s'' ''%s'' : ''%s'' is expected to be (%d x %d)', st, nameidxstr(name, idx), 'c', Mc, 1);
        end
        if ~isempty(p.k) && any(p.k) && Nk ~= N && Nk ~= 1
            error('opt_model.set_params: ''%s'' ''%s'' : ''%s'' is expected to be (%d x %d)', st, nameidxstr(name, idx), 'k', N, 1);
        end

        %% check consistency of Q, c, k and vs
        if u.Q || u.c || u.vs
            p.vs = om.varsets_cell2struct(p.vs);
            nv = om.varsets_len(p.vs);      %% number of variables
            if nx ~= nv
                error('opt_model.set_params: for ''%s'' ''%s'' dimensions of ''Q'', ''c'', ''k'' (%d) must be consistent with ''vs'' (%d)', st, nameidxstr(name, idx), nx, nv);
            end
        end

        %% assign new parameters
        if isempty(idx)     %% simple named set
            for k = 1:length(default_params)
                pn = default_params{k};     %% param name
                if u.(pn)   %% assign new val for this parameter
                    om.qdc.data.(pn).(name) = p.(pn);
                end
            end
        else                %% indexed named set
            for k = 1:length(default_params)
                pn = default_params{k};     %% param name
                if u.(pn)   %% assign new val for this parameter
                    om.qdc.data.(pn) = subsasgn(om.qdc.data.(pn), sc, p.(pn));
                end
            end
        end

        %% clear cached parameters
        om.qdc.params = [];
    case {'nle', 'nli'}
        %% get current parameters
        if isempty(idx)
            [N0, fcn, hess, vs, include] = om.params_nln_constraint(st(3) == 'e', name, idx);
        else
            [N0, fcn, hess, vs] = om.params_nln_constraint(st(3) == 'e', name, idx);
            include = '';
        end
        if isempty(vs), vs = {vs}; end
        p = struct('N', N0, 'fcn', fcn, 'hess', hess, 'vs', vs);    %% current parameters
        u = struct('N',  0, 'fcn',   0, 'hess',    0, 'vs',  0);    %% which ones to update

        %% replace with new parameters
        for k = 1:np
            p.(params{k}) = vals{k};
            u.(params{k}) = 1;
        end
        N = p.N;

        %% set missing default params for 'all'
        if is_all
            u.N    = 1;         %% always update N
            u.fcn  = 1;         %% alwaus update fcn
            u.hess = 1;         %% alwaus update hess
            if np < 4
                p.vs = {};
                u.vs = 1;       %% update vs
            end
        end

        %% check consistency of parameters
        %% no dimension change unless 'all'
        if N ~= N0 && ~is_all
            error('opt_model.set_params: dimension change for ''%s'' ''%s'' not allowed except for ''all''', st, nameidxstr(name, idx));
        end

        %% included constraints not yet implemented
        if ~isempty(include)
            error('opt_model.set_params: modifications for ''%s'' ''%s'' not (yet) supported since it includes evaluation of other constraints', st, nameidxstr(name, idx));
        end

        %% convert vs to struct
        if u.vs
            p.vs = om.varsets_cell2struct(p.vs);
        end

        %% assign new parameters
        if isempty(idx)     %% simple named set
            for k = 2:length(default_params)
                pn = default_params{k};     %% param name
                if u.(pn)   %% assign new val for this parameter
                    om.(st).data.(pn).(name) = p.(pn);
                end
            end
        else                %% indexed named set
            for k = 2:length(default_params)
                pn = default_params{k};     %% param name
                if u.(pn)   %% assign new val for this parameter
                    om.(st).data.(pn) = subsasgn(om.(st).data.(pn), sc, p.(pn));
                end
            end
        end
    case 'nlc'
        %% get current parameters
        [N0, fcn, vs] = om.params_nln_cost(name, idx);
        if isempty(vs), vs = {vs}; end
        p = struct('N', N0, 'fcn', fcn, 'vs', vs);  %% current parameters
        u = struct('N',  0, 'fcn',   0, 'vs',  0);  %% which ones to update

        %% replace with new parameters
        for k = 1:np
            p.(params{k}) = vals{k};
            u.(params{k}) = 1;
        end
        N = p.N;

        %% set missing default params for 'all'
        if is_all
            u.N   = 1;          %% always update N
            u.fcn = 1;          %% alwaus update fcn
            if np < 3
                p.vs = {};
                u.vs = 1;       %% update vs
            end
        end

        %% check consistency of parameters
        %% no dimension change unless 'all'
        if N ~= N0 && ~is_all
            error('opt_model.set_params: dimension change for ''%s'' ''%s'' not allowed except for ''all''', st, nameidxstr(name, idx));
        end

        %% vector valued costs not yet implemented
        if N ~= 1
            error('opt_model.set_params: vector value for ''%s'' ''%s'' not yet implemented', st, nameidxstr(name, idx));
        end

        %% assign new parameters
        if isempty(idx)     %% simple named set
            for k = 2:length(default_params)
                pn = default_params{k};     %% param name
                if u.(pn)   %% assign new val for this parameter
                    om.nlc.data.(pn).(name) = p.(pn);
                end
            end
        else                %% indexed named set
            for k = 2:length(default_params)
                pn = default_params{k};     %% param name
                if u.(pn)   %% assign new val for this parameter
                    om.nlc.data.(pn) = subsasgn(om.nlc.data.(pn), sc, p.(pn));
                end
            end
        end
    otherwise
        error('opt_model.set_params: ''%s'' is not a valid SET_TYPE', st);
end

%% update dimensions and indexing, if necessary
dN = N - N0;
if is_all && dN
    %% calls to substruct() are relatively expensive, so we pre-build the
    %% struct for addressing num array fields
    %% sn = substruct('.', name, '()', idx);
    sn = struct('type', {'.', '()'}, 'subs', {'', 1});  %% num array field
    om_ff = om.(st);
    update = 0;             %% not yet reached set being updated
    update_i1 = 0;          %% flag to indicate whether to update i1
    for k = 1:om_ff.NS
        o = om_ff.order(k);
        if ~update && strcmp(o.name, name) && isequal(o.idx, idx)
            update = 1;     %% arrived at set being updated
        end
        if update
            if isempty(o.idx)   %% simple named set
                if update_i1
                    om_ff.idx.i1.(o.name) = om_ff.idx.i1.(o.name) + dN;
                else
                    om_ff.idx.N.(o.name) = om_ff.idx.N.(o.name) + dN;
                end
                om_ff.idx.iN.(o.name) = om_ff.idx.iN.(o.name) + dN;
            else                %% indexed named set
                sn(1).subs = o.name;
                sn(2).subs = o.idx;
                if update_i1
                    v = subsref(om_ff.idx.i1, sn);
                    om_ff.idx.i1 = subsasgn(om_ff.idx.i1, sn, v + dN);
                else
                    v = subsref(om_ff.idx.N, sn);
                    om_ff.idx.N = subsasgn(om_ff.idx.N, sn, v + dN);
                end
                v = subsref(om_ff.idx.iN, sn);
                om_ff.idx.iN = subsasgn(om_ff.idx.iN, sn, v + dN);
            end
            update_i1 = 1;  %% update i1 from here on out
        end
    end
    om_ff.N = om_ff.N + dN;
    om.(st) = om_ff;
end


function str = nameidxstr(name, idx)
str = sprintf('%s%s', name, idxstr(idx));

function str = idxstr(idx)
if isempty(idx)
    str = '';
elseif length(idx) == 1
    str = sprintf('(%d)', idx{1});
else
    str = ['(' sprintf('%d', idx{1}) sprintf(',%d', idx{2:end}) ')'];
end
