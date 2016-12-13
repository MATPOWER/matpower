function om = add_costs(om, name, idx, varargin)
%ADD_COSTS  Adds a set of user costs to the model.
%   OM = ADD_COSTS(OM, NAME, CP);
%   OM = ADD_COSTS(OM, NAME, CP, VARSETS);
%   OM = ADD_COSTS(OM, NAME, DIM_LIST);
%   OM = ADD_COSTS(OM, NAME, IDX_LIST, CP);
%   OM = ADD_COSTS(OM, NAME, IDX_LIST, CP, VARSETS);
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
%       om = add_costs(om, 'usr1', cp1, {'Pg', 'Qg', 'z'});
%       om = add_costs(om, 'usr2', cp2, {'Vm', 'Pg', 'Qg', 'z'});
%
%       om = add_costs(om, 'c', {2, 3});
%       for i = 1:2
%         for j = 1:3
%           om = add_costs(om, 'c', {i, j}, cp(i,j), ...);
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
%   See also OPT_MODEL, BUILD_COST_PARAMS, GET_COST_PARAMS, COMPUTE_COST.

%   MATPOWER
%   Copyright (c) 2008-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if iscell(idx)
    if ~isempty(varargin)       %% indexed named set
        % (calls to substruct() are relatively expensive ...
        % s1 = substruct('.', name, '()', idx);
        % s2 = substruct('.', name, '{}', idx);
        % ... so replace them with these more efficient lines)
        s1 = struct('type', {'.', '()'}, 'subs', {name, idx});
        s2 = s1;
        s2(2).type = '{}';
        
        %% prevent duplicate named cost sets
        if subsref(om.cost.idx.i1, s1) ~= 0
            str = '%d'; for m = 2:length(idx), str = [str ',%d']; end
            nname = sprintf(['%s(' str, ')'], name, idx{:});
            error('@opt_model/add_costs: cost set named ''%s'' already exists', nname);
        end
        
        cp = varargin{1};
        args = varargin(2:end);
    else                        %% just setting dimensions for indexed set
        %% prevent duplicate named cost sets
        if isfield(om.cost.idx.N, name)
            error('@opt_model/add_costs: cost set named ''%s'' already exists', name);
        end
        
        cp = [];
        args = {};
    end
else                            %% simple named set
    %% prevent duplicate named cost sets
    if isfield(om.cost.idx.N, name)
        error('@opt_model/add_costs: cost set named ''%s'' already exists', name);
    end
    
    cp = idx;
    args = varargin;
    idx = {};
end

if isempty(cp)                  %% just setting dimensions for indexed set
    %% use column vector if single dimension
    if length(idx) == 1
        idx = {idx{:}, 1};
    end
    
    %% add info about this cost set
    om.cost.idx.i1.(name)  = zeros(idx{:}); %% starting index
    om.cost.idx.iN.(name)  = zeros(idx{:}); %% ending index
    om.cost.idx.N.(name)   = zeros(idx{:}); %% number of costs (nw)
    om.cost.data.N.(name)  = cell(idx{:});
    om.cost.data.Cw.(name) = cell(idx{:});
    om.cost.data.vs.(name) = cell(idx{:});
else
    if isempty(args)
        varsets = {};
    else
        varsets = args{1};
    end
    if ~isempty(varsets) && iscell(varsets)
        empty_cells = cell(1, length(varsets));
        [empty_cells{:}] = deal({});    %% empty cell arrays
        varsets = struct('name', varsets, 'idx', empty_cells);
    end
    if isfield(cp, 'N')
        [nw, nx] = size(cp.N);
    else
        nw = length(cp.Cw);
        nx = nw;
        cp.N = speye(nw, nx);
    end
    
    %% check sizes
    if isempty(varsets)
        nv = om.var.N;
    else
        nv = 0;
        s = struct('type', {'.', '()'}, 'subs', {'', 1});
        for k = 1:length(varsets)
            % (calls to substruct() are relatively expensive ...
            % s = substruct('.', varsets(k).name, '()', varsets(k).idx);
            % ... so replace it with these more efficient lines)
            s(1).subs = varsets(k).name;
            s(2).subs = varsets(k).idx;
            nv = nv + subsref(om.var.idx.N, s);
        end
    end
    if nx ~= nv
        if nw == 0
            cp.N = sparse(nw, nx);
        else
            error('@opt_model/add_costs: number of columns in N (%d x %d) does not match\nnumber of variables (%d)\n', nw, nx, nv);
        end
    end
    if size(cp.Cw, 1) ~= nw
        error('@opt_model/add_costs: number of rows of Cw (%d x %d) and N (%d x %d) must match\n', size(cp.Cw), nw, nx);
    end
    if isfield(cp, 'H') && (size(cp.H, 1) ~= nw || size(cp.H, 2) ~= nw)
        error('@opt_model/add_costs: both dimensions of H (%d x %d) must match the number of rows in N (%d x %d)\n', size(cp.H), nw, nx);
    end
    if isfield(cp, 'dd') && size(cp.dd, 1) ~= nw
        error('@opt_model/add_costs: number of rows of dd (%d x %d) and N (%d x %d) must match\n', size(cp.dd), nw, nx);
    end
    if isfield(cp, 'rh') && size(cp.rh, 1) ~= nw
        error('@opt_model/add_costs: number of rows of rh (%d x %d) and N (%d x %d) must match\n', size(cp.rh), nw, nx);
    end
    if isfield(cp, 'kk') && size(cp.kk, 1) ~= nw
        error('@opt_model/add_costs: number of rows of kk (%d x %d) and N (%d x %d) must match\n', size(cp.kk), nw, nx);
    end
    if isfield(cp, 'mm') && size(cp.mm, 1) ~= nw
        error('@opt_model/add_costs: number of rows of mm (%d x %d) and N (%d x %d) must match\n', size(cp.mm), nw, nx);
    end
    
    if isempty(idx)     %% simple named set
        %% add info about this user cost set
        om.cost.idx.i1.(name)  = om.cost.N + 1;     %% starting index
        om.cost.idx.iN.(name)  = om.cost.N + nw;    %% ending index
        om.cost.idx.N.(name)   = nw;                %% number of costs (nw)
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
        
        %% update number of vars and var sets
        om.cost.N  = om.cost.idx.iN.(name);
        om.cost.NS = om.cost.NS + 1;
        
        %% add to ordered list of var sets
        om.cost.order(om.cost.NS).name = name;
        om.cost.order(om.cost.NS).idx  = {};
    else                %% indexed named set
        %% add info about this user cost set
        om.cost.idx.i1  = subsasgn(om.cost.idx.i1, s1, om.cost.N + 1);  %% starting index
        om.cost.idx.iN  = subsasgn(om.cost.idx.iN, s1, om.cost.N + nw); %% ending index
        om.cost.idx.N   = subsasgn(om.cost.idx.N,  s1, nw);             %% number of costs (nw)
        
        om.cost.data.N  = subsasgn(om.cost.data.N,  s2, cp.N);
        om.cost.data.Cw = subsasgn(om.cost.data.Cw, s2, cp.Cw);
        om.cost.data.vs = subsasgn(om.cost.data.vs, s2, varsets);
        if isfield(cp, 'H')
            om.cost.data.H = subsasgn(om.cost.data.H, s2, cp.H);
        end
        if isfield(cp, 'dd')
            om.cost.data.dd = subsasgn(om.cost.data.dd, s2, cp.dd);
        end
        if isfield(cp, 'rh')
            om.cost.data.rh = subsasgn(om.cost.data.rh, s2, cp.rh);
        end
        if isfield(cp, 'kk')
            om.cost.data.kk = subsasgn(om.cost.data.kk, s2, cp.kk);
        end
        if isfield(cp, 'mm')
            om.cost.data.mm = subsasgn(om.cost.data.mm, s2, cp.mm);
        end
        
        %% update number of costs and cost sets
        om.cost.N  = subsref(om.cost.idx.iN, s1);
        om.cost.NS = om.cost.NS + 1;
        
        %% add to ordered list of cost sets
        om.cost.order(om.cost.NS).name = name;
        om.cost.order(om.cost.NS).idx  = idx;
    end
end
