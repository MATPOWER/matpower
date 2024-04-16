classdef (Abstract) mm_shared_opf_legacy < handle
% mp.mm_shared_opf_legacy - Mixin class for legacy optimal power flow (OPF) **math model** objects.
%
% An abstract mixin class inherited by optimal power flow (OPF) **math model**
% objects that need to handle legacy user customization mechanisms.

%   MATPOWER
%   Copyright (c) 2008-2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        cost = [];
        mpc = struct();
    end     %% properties

    methods
        function obj = def_set_types_legacy(obj)
            %

            obj.set_types = struct(...
                    'var',  'VARIABLES', ...
                    'lin',  'LINEAR CONSTRAINTS', ...
                    'nle',  'NONLIN EQ CONSTRAINTS', ...
                    'nli',  'NONLIN INEQ CONSTRAINTS', ...
                    'qdc',  'QUADRATIC COSTS', ...
                    'nlc',  'GEN NONLIN COSTS', ...
                    'cost', 'LEGACY COSTS'  ...
                );
        end

        function obj = init_set_types_legacy(obj)
            %

            %% finish initializing data structures for each type
            es = struct();  %% empty struct
            obj.cost.data = struct( ...
                'N', es, ...
                'H', es, ...
                'Cw', es, ...
                'dd', es, ...
                'rh', es, ...
                'kk', es, ...
                'mm', es, ...
                'vs', es );
            obj.cost.params = [];
        end

        function mpc = get_mpc(om)
            %

            mpc = om.mpc;
        end

        function obj = build_legacy(obj, nm, dm, mpopt)
            %

            if strcmp(obj.form_tag, 'dc') && toggle_softlims(obj.mpc, 'status')
                %% user data required by toggle_softlims
                branch_nme = nm.elements.branch;
                [Bbr, pbr] = branch_nme.get_params(1:branch_nme.nk, {'B', 'p'});
                obj.userdata.Bf = Bbr * branch_nme.C';
                obj.userdata.Pfinj = pbr;
            end

            %% execute userfcn callbacks for 'formulation' stage
            if isfield(obj.mpc, 'userfcn')
                userfcn = obj.mpc.userfcn;
            else
                userfcn = [];
            end
            obj = run_userfcn(userfcn, 'formulation', obj, mpopt);
        end

        function add_legacy_user_vars(obj, nm, dm, mpopt)
            %

            %% save data
            obj.userdata.user_vars = obj.legacy_user_var_names();

            %% add any user-defined vars
            if isfield(dm.userdata.legacy_opf_user_mods, 'z')
                z = dm.userdata.legacy_opf_user_mods.z;
                if z.nz > 0
                    obj.add_var('z', z.nz, z.z0, z.zl, z.zu);
                    obj.userdata.user_vars{end+1} = 'z';
                end
            end
        end

        function add_legacy_user_costs(obj, nm, dm, dc)
            %

            if isfield(dm.userdata.legacy_opf_user_mods, 'cost') && ...
                    dm.userdata.legacy_opf_user_mods.cost.nw
                user_cost = dm.userdata.legacy_opf_user_mods.cost;
                uv = obj.get_userdata('user_vars');
                obj.add_legacy_cost('usr', user_cost, uv);

                %% implement legacy user costs using quadratic or general non-linear costs
                cp = obj.params_legacy_cost();  %% construct/fetch the parameters
                [N, H, Cw, rh, m] = deal(cp.N, cp.H, cp.Cw, cp.rh, cp.mm);
                [nw, nx] = size(N);
                if nw
                    if any(cp.dd ~= 1) || any(cp.kk)    %% not simple quadratic form
                        if dc                           %% (includes "dead zone" or
                            if any(cp.dd ~= 1)          %%  quadratic "penalty")
                                error('mp.mm_shared_opf_legacy.add_legacy_user_costs: DC OPF can only handle legacy user-defined costs with d = 1');
                            end
                            if any(cp.kk)
                                error('mp.mm_shared_opf_legacy.add_legacy_user_costs: DC OPF can only handle legacy user-defined costs with no "dead zone", i.e. k = 0');
                            end
                        else
                            %% use general nonlinear cost to implement legacy user cost
                            legacy_cost_fcn = @(x)opf_legacy_user_cost_fcn(x, cp);
                            obj.add_nln_cost('usr', 1, legacy_cost_fcn);
                        end
                    else                                %% simple quadratic form
                        %% use a quadratic cost to implement legacy user cost
                        %% f = 1/2 * w'*H*w + Cw'*w, where w = diag(m)*(N*x - rh)
                        %% Let: MN = diag(m)*N
                        %%      MR = M * rh
                        %%      HMR  = H  * MR;
                        %%      HtMR = H' * MR;
                        %%  =>   w = MN*x - MR
                        %% f = 1/2 * (MN*x - MR)'*H*(MN*x - MR) + Cw'*(MN*x - MR)
                        %%   = 1/2 * x'*MN'*H*MN*x +
                        %%          (Cw'*MN - 1/2 * MR'*(H+H')*MN)*x +
                        %%          1/2 * MR'*H*MR - Cw'*MR
                        %%   = 1/2 * x'*Q*w + c'*x + k

                        M    = sparse(1:nw, 1:nw, m, nw, nw);
                        MN   = M * N;
                        MR   = M * rh;
                        HMR  = H  * MR;
                        HtMR = H' * MR;
                        Q = MN' * H * MN;
                        c = full(MN' * (Cw - 1/2*(HMR+HtMR)));
                        k = (1/2 * HtMR - Cw)' * MR;
                        obj.add_quad_cost('usr', Q, c, k);
                    end
                end
            end
        end

        function add_legacy_user_constraints(obj, nm, dm, mpopt)
            %

            %% user-defined linear constraints
            if isfield(dm.userdata.legacy_opf_user_mods, 'lin')
                lin = dm.userdata.legacy_opf_user_mods.lin;
                if lin.nlin
                    uv = obj.get_userdata('user_vars');
                    obj.add_lin_constraint('usr', lin.A, lin.l, lin.u, uv);
                end
            end
        end

        function add_legacy_user_constraints_ac(obj, nm, dm, mpopt)
            %

            obj.add_legacy_user_constraints(nm, dm, mpopt);

            if ~isempty(dm.userdata.legacy_opf_user_mods)
                uc = dm.userdata.legacy_opf_user_mods.nlc;
                for k = 1:length(uc)
                    obj.add_nln_constraint(uc{k}{:});
                end
            end
        end

        function om = add_legacy_cost(om, name, idx, varargin)
            % add_legacy_cost - Add a set of user costs to the model
            % ::
            %
            %   mm.add_legacy_cost(name, cp)
            %   mm.add_legacy_cost(name, idx, varsets)
            %   mm.add_legacy_cost(name, idx_list, cp)
            %   mm.add_legacy_cost(name, idx_list, cp, varsets)

            %ADD_LEGACY_COST  Adds a set of user costs to the model.
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
                    error('mp.mm_shared_opf_legacy.add_legacy_cost: number of columns in N (%d x %d) does not match\nnumber of variables (%d)\n', nw, nx, nv);
                end
            end
            if size(cp.Cw, 1) ~= nw
                error('mp.mm_shared_opf_legacy.add_legacy_cost: number of rows of Cw (%d x %d) and N (%d x %d) must match\n', size(cp.Cw), nw, nx);
            end
            if isfield(cp, 'H') && (size(cp.H, 1) ~= nw || size(cp.H, 2) ~= nw)
                error('mp.mm_shared_opf_legacy.add_legacy_cost: both dimensions of H (%d x %d) must match the number of rows in N (%d x %d)\n', size(cp.H), nw, nx);
            end
            if isfield(cp, 'dd') && size(cp.dd, 1) ~= nw
                error('mp.mm_shared_opf_legacy.add_legacy_cost: number of rows of dd (%d x %d) and N (%d x %d) must match\n', size(cp.dd), nw, nx);
            end
            if isfield(cp, 'rh') && size(cp.rh, 1) ~= nw
                error('mp.mm_shared_opf_legacy.add_legacy_cost: number of rows of rh (%d x %d) and N (%d x %d) must match\n', size(cp.rh), nw, nx);
            end
            if isfield(cp, 'kk') && size(cp.kk, 1) ~= nw
                error('mp.mm_shared_opf_legacy.add_legacy_cost: number of rows of kk (%d x %d) and N (%d x %d) must match\n', size(cp.kk), nw, nx);
            end
            if isfield(cp, 'mm') && size(cp.mm, 1) ~= nw
                error('mp.mm_shared_opf_legacy.add_legacy_cost: number of rows of mm (%d x %d) and N (%d x %d) must match\n', size(cp.mm), nw, nx);
            end

            %% add the legacy cost set
            om.add_named_set('cost', name, idx, nw, cp, varsets);
        end

        function [f, df, d2f] = eval_legacy_cost(om, x, name, idx)
            % eval_legacy_cost - Evaluate individual or full set of legacy user costs.
            % ::
            %
            %   f = mm.eval_legacy_cost(x ...)
            %   [f, df] = mm.eval_legacy_cost(x ...)
            %   [f, df, d2f] = mm.eval_legacy_cost(x ...)
            %   [f, df, d2f] = mm.eval_legacy_cost(x, name)
            %   [f, df, d2f] = mm.eval_legacy_cost(x, name, idx_list)

            %EVAL_LEGACY_COST  Evaluates individual or full set of legacy user costs.
            %   F = OM.EVAL_LEGACY_COST(X ...)
            %   [F, DF] = OM.EVAL_LEGACY_COST(X ...)
            %   [F, DF, D2F] = OM.EVAL_LEGACY_COST(X ...)
            %   [F, DF, D2F] = OM.EVAL_LEGACY_COST(X, NAME)
            %   [F, DF, D2F] = OM.EVAL_LEGACY_COST(X, NAME, IDX_LIST)
            %   Evaluates an individual named set or the full set of legacy user
            %   costs and their derivatives for a given value of the optimization vector
            %   X, based on costs added by ADD_LEGACY_COST.
            %
            %   Example:
            %       [f, df, d2f] = om.eval_legacy_cost(x)
            %       [f, df, d2f] = om.eval_legacy_cost(x, name)
            %       [f, df, d2f] = om.eval_legacy_cost(x, name, idx)
            %
            %   See also OPT_MODEL, ADD_LEGACY_COST, PARAMS_LEGACY_COST.

            if om.cost.N
                done = 0;

                %% collect cost parameters
                if nargin < 3                       %% full set
                    [cp, vs] = om.params_legacy_cost();
                elseif nargin < 4 || isempty(idx)   %% name, no idx provided
                    dims = size(om.cost.idx.i1.(name));
                    if prod(dims) == 1              %% simple named set
                        [cp, vs] = om.params_legacy_cost(name);
                    elseif nargout == 1             %% indexing required, recurse
                        f = 0;          %% initialize cumulative cost
                        idx = num2cell(ones(size(dims))); %% initialize idx
                        while ~done     %% call eval_legacy_cost() recursively
                            f = f + om.eval_legacy_cost(x, name, idx);

                            %% increment idx
                            D = length(dims);
                            idx{D} = idx{D} + 1;    %% increment last dimension
                            for d = D:-1:2          %% increment next dimension, if necessary
                                if idx{d} > dims(d)
                                    idx{d} = 1;
                                    idx{d-1} = idx{d-1} + 1;
                                end
                            end
                            if idx{1} > dims(1)     %% check if done
                                done = 1;
                            end
                        end
                    else
                        error('mp.mm_shared_opf_legacy.eval_legacy_cost: legacy cost set ''%s'' requires an IDX_LIST arg when requesting DF output', name)
                    end
                else                                %% indexed named set
                    [cp, vs] = om.params_legacy_cost(name, idx);
                end

                if ~done
                    %% assemble appropriately-sized x vector
                    xx = om.varsets_x(x, vs, 'vector');

                    %% compute function & derivatives
                    if nargout == 1
                        f = opf_legacy_user_cost_fcn(xx, cp);
                    elseif nargout == 2
                        [f, df] = opf_legacy_user_cost_fcn(xx, cp);
                    else    %% nargout == 3
                        [f, df, d2f] = opf_legacy_user_cost_fcn(xx, cp);
                    end
                end
            else
                f = 0;
                if nargout > 1
                    df = [];
                    if nargout > 2
                        d2f = [];
                    end
                end
            end
        end

        function [cp, vs, i1, iN] = params_legacy_cost(om, name, idx)
            % params_legacy_cost - Return cost parameters for legacy user-defined costs.
            % ::
            %
            %   cp = mm.params_legacy_cost()
            %   cp = mm.params_legacy_cost(name)
            %   cp = mm.params_legacy_cost(name, idx)
            %   [cp, vs] = mm.params_legacy_cost(...)
            %   [cp, vs, i1, iN] = mm.params_legacy_cost(...)

            %PARAMS_LEGACY_COST  Returns cost parameters for legacy user-defined costs.
            %   CP = OM.PARAMS_LEGACY_COST()
            %   CP = OM.PARAMS_LEGACY_COST(NAME)
            %   CP = OM.PARAMS_LEGACY_COST(NAME, IDX_LIST)
            %   [CP, VS] = OM.PARAMS_LEGACY_COST(...)
            %   [CP, VS, I1, IN] = OM.PARAMS_LEGACY_COST(...)
            %
            %   With no input parameters, it assembles and returns the parameters
            %   for the aggregate legacy user-defined cost from all legacy cost sets
            %   added using ADD_LEGACY_COST. The values of these parameters are cached
            %   for subsequent calls. The parameters are contained in the struct CP,
            %   described below.
            %
            %   If a NAME is provided then it simply returns parameter struct CP for the
            %   corresponding named set. Likewise for indexed named sets specified
            %   by NAME and IDX_LIST.
            %
            %   An optional 2nd output argument VS indicates the variable sets used by
            %   this cost set. The size of CP.N will be consistent with VS.
            %
            %   If NAME is provided, optional 3rd and 4th output arguments I1 and IN
            %   indicate the starting and ending row indices of the corresponding
            %   cost set in the full aggregate cost matrix CP.N.
            %
            %   Let X refer to the vector formed by combining the corresponding varsets VS,
            %   and F_U(X, CP) be the cost at X corresponding to the cost parameters
            %   contained in CP, where CP is a struct with the following fields:
            %       N      - nw x nx sparse matrix (optional, identity matrix by default)
            %       Cw     - nw x 1 vector
            %       H      - nw x nw sparse matrix (optional, all zeros by default)
            %       dd, mm - nw x 1 vectors (optional, all ones by default)
            %       rh, kk - nw x 1 vectors (optional, all zeros by default)
            %
            %   These parameters are used as follows to compute F_U(X, CP)
            %
            %       R  = N*X - rh
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
            %       F_U(X, CP) = 1/2 * w'*H*w + Cw'*w
            %
            %   See also OPT_MODEL, ADD_LEGACY_COST, EVAL_LEGACY_COST.

            if nargin > 1       %% individual set
                if nargin < 3
                    idx = {};
                end
                if isempty(idx)                 %% name, no index provided
                    if numel(om.cost.idx.i1.(name)) == 1    %% simple named set
                        N  = om.cost.data.N.( name);
                        Cw = om.cost.data.Cw.(name);
                        [nw, nx] = size(N);
                        e1 = ones(nw, 1);
                        e0 = zeros(nw, 1);
                        cp = struct(...
                            'N',    N, ...
                            'Cw',   Cw, ...
                            'H',    sparse(nw, nw), ... %% default => no quadratic term
                            'dd',   e1, ...             %% default => linear
                            'rh',   e0, ...             %% default => no shift
                            'kk',   e0, ...             %% default => no dead zone
                            'mm',   e1 );               %% default => no scaling
                        if isfield(om.cost.data.H, name) && ~isempty(om.cost.data.H.(name))
                            cp.H = om.cost.data.H.(name);
                        end
                        if isfield(om.cost.data.dd, name) && ~isempty(om.cost.data.dd.(name))
                            cp.dd = om.cost.data.dd.(name);
                        end
                        if isfield(om.cost.data.rh, name) && ~isempty(om.cost.data.rh.(name))
                            cp.rh = om.cost.data.rh.(name);
                        end
                        if isfield(om.cost.data.kk, name) && ~isempty(om.cost.data.kk.(name))
                            cp.kk = om.cost.data.kk.(name);
                        end
                        if isfield(om.cost.data.mm, name) && ~isempty(om.cost.data.mm.(name))
                            cp.mm = om.cost.data.mm.(name);
                        end
                        if nargout > 1
                            vs = om.cost.data.vs.(name);
                            if nargout > 3
                                i1 = om.cost.idx.i1.(name);     %% starting row index
                                iN = om.cost.idx.iN.(name);     %% ending row index
                            end
                        end
                    else                                    %% indexing required
                        error('mp.mm_shared_opf_legacy.params_legacy_cost: legacy cost set ''%s'' requires an IDX_LIST arg', name);
                    end
                else                            %% indexed named set
                    % (calls to substruct() are relatively expensive ...
                    % s = substruct('.', name, '{}', idx);
                    % ... so replace it with these more efficient lines)
                    sc = struct('type', {'.', '{}'}, 'subs', {name, idx});
                    N  = subsref(om.cost.data.N,  sc);
                    Cw = subsref(om.cost.data.Cw, sc);
                    [nw, nx] = size(N);
                    e1 = ones(nw, 1);
                    e0 = zeros(nw, 1);
                    cp = struct(...
                        'N',    N, ...
                        'Cw',   Cw, ...
                        'H',    sparse(nw, nw), ... %% default => no quadratic term
                        'dd',   e1, ...             %% default => linear
                        'rh',   e0, ...             %% default => no shift
                        'kk',   e0, ...             %% default => no dead zone
                        'mm',   e1 );               %% default => no scaling
                    if isfield(om.cost.data.H, name) && ~isempty(subsref(om.cost.data.H, sc))
                        cp.H = subsref(om.cost.data.H,  sc);
                    end
                    if isfield(om.cost.data.dd, name) && ~isempty(subsref(om.cost.data.dd, sc))
                        cp.dd = subsref(om.cost.data.dd, sc);
                    end
                    if isfield(om.cost.data.rh, name) && ~isempty(subsref(om.cost.data.rh, sc))
                        cp.rh = subsref(om.cost.data.rh, sc);
                    end
                    if isfield(om.cost.data.kk, name) && ~isempty(subsref(om.cost.data.kk, sc))
                        cp.kk = subsref(om.cost.data.kk, sc);
                    end
                    if isfield(om.cost.data.mm, name) && ~isempty(subsref(om.cost.data.mm, sc))
                        cp.mm = subsref(om.cost.data.mm, sc);
                    end
                    if nargout > 1
                        vs = subsref(om.cost.data.vs, sc);
                        if nargout > 3
                            sn = sc; sn(2).type = '()';         %% num array field
                            i1 = subsref(om.cost.idx.i1, sn);   %% starting row index
                            iN = subsref(om.cost.idx.iN, sn);   %% ending row index
                        end
                    end
                end
            else                %% aggregate
                cache = om.cost.params;
                if isempty(cache)       %% build the aggregate
                    nx = om.var.N;          %% number of variables
                    nw = om.cost.N;         %% number of cost rows
                    Nt = sparse(nx, nw);    %% use N transpose for speed
                    Cw = zeros(nw, 1);
                    H = sparse(nw, nw);     %% default => no quadratic term
                    dd = ones(nw, 1);       %% default => linear
                    rh = Cw;                %% default => no shift
                    kk = Cw;                %% default => no dead zone
                    mm = dd;                %% default => no scaling

                    %% fill in each piece
                    for k = 1:om.cost.NS
                        name = om.cost.order(k).name;
                        idx  = om.cost.order(k).idx;
                        [cp, vs, i1, iN] = om.params_legacy_cost(name, idx);
                        [mk, nk] = size(cp.N);      %% size of Nk
                        if mk
                            Nkt_full = sparse(nx, nw);
                            if isempty(vs)
                                if nk == nx     %% full size
                                    Nkt_full(:, i1:iN) = cp.N';
                                else            %% vars added since adding this cost set
                                    Nk_all_cols = sparse(mk, nx);
                                    Nk_all_cols(:, 1:nk) = cp.N;
                                    Nkt_full(:, i1:iN) = Nk_all_cols';
                                end
                            else
                                jj = om.varsets_idx(vs);    %% indices for var set
                                Nk_all_cols = sparse(mk, nx);
                                Nk_all_cols(:, jj) = cp.N;
                                Nkt_full(:, i1:iN) = Nk_all_cols';
                            end
                            Nt = Nt + Nkt_full;
                            Cw(i1:iN) = cp.Cw;
                            H_all_cols = sparse(mk, nw);
                            H_all_cols(:, i1:iN) = cp.H';
                            H(:, i1:iN) = H_all_cols';
                            dd(i1:iN) = cp.dd;
                            rh(i1:iN) = cp.rh;
                            kk(i1:iN) = cp.kk;
                            mm(i1:iN) = cp.mm;
                        end
                    end

                    %% cache aggregated parameters
                    om.cost.params = struct( ...
                        'N', Nt', 'Cw', Cw, 'H', H, 'dd', dd, 'rh', rh, 'kk', kk, 'mm', mm );
                    cp = om.cost.params;
                else                    %% return cached values
                    cp = cache;
                end
                if nargout > 1
                    vs = {};
                end
            end
        end
    end     %% methods

    methods (Access=protected)
        function om = add_named_set_legacy(om, set_type, name, idx, N, varargin)
            %ADD_NAMED_SET  Adds a named set of variables/constraints/costs to the model.
            %
            %   -----  PRIVATE METHOD  -----
            %
            %   Legacy Cost Set
            %       OM.ADD_NAMED_SET_LEGACY('cost', NAME, IDX_LIST, N, CP, VARSETS);
            %
            %   See also OPT_MODEL, ADD_VAR, ADD_LIN_CONSTRAINT, ADD_NLN_CONSTRAINT
            %            ADD_QUAD_COST, ADD_NLN_COST and ADD_LEGACY_COST.

            switch set_type
                case 'cost'         %% cost set
                    %% add type-specific data for named set
                    om_ff = om.cost;
                    om.cost = [];

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
                        %% calls to substruct() are relatively expensive, so we pre-build the
                        %% struct for addressing cell array fields
                        %% sc = substruct('.', name, '{}', idx);
                        sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field

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
                    om.cost = om_ff;
            end
        end
    end     %% methods (Access=protected)
end         %% classdef
