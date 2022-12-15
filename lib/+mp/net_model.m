classdef (Abstract) net_model < mp.nm_element & mp.element_container & mp_idx_manager% & mp.form
%MP.NET_MODEL  Abstract base class for MATPOWER network model
%   Explicitly a subclass of MP.NM_ELEMENT, MP_IDX_MANAGER and
%   MP.ELEMENT_CONTAINER, and implicitly assumed to be a subclass of MP.FORM
%   as well.

%   MATPOWER
%   Copyright (c) 2019-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        the_np = 0;     %% number of ports
        the_nz = 0;     %% number of non-voltage states
        nv = 0;         %% total number of (real) v variables
        node = [];
        port = [];
        state = [];
    end

    methods
        function name = name(obj)
            name = 'network';
        end

        function np = np(obj)
            np = obj.the_np;    %% number of ports
        end

        function nz = nz(obj)
            nz = obj.the_nz;    %% number of (possibly complex) non-voltage states
        end

        function obj = build(obj, dm)
            %% Due to a bug related to inheritance in constructors in
            %% Octave 5.2 and earlier (https://savannah.gnu.org/bugs/?52614),
            %% INIT_SET_TYPES() cannot be called directly in the
            %% MP_IDX_MANAGER constructor, as desired.
            %%
            %% WORKAROUND:  Initialize MP_IDX_MANAGER fields here, if needed,
            %%              after object construction, but before object use.
            if isempty(obj.node)        %% only if not already initialized
                obj.init_set_types();
            end

            %% create element objects for each class with data
            obj.nk = 1;
            obj.elements = mp.mapped_array();
            for c = obj.element_classes
                nme = c{1}();       %% element constructor
                if nme.count(dm)
                    obj.elements.add_elements(nme, nme.name);
                    obj.the_np = obj.the_np + nme.np * nme.nk;  %% number of ports
                    obj.the_nz = obj.the_nz + nme.nz * nme.nk;  %% number of z_ vars
                end
            end

            if obj.np ~= 0      %% skip for empty model
                %% create nodes and node voltage state variables
                obj.add_nodes(obj, dm);

                %% create non-voltage states and corresponding state variables
                obj.add_states(obj, dm);

                %% build params
                obj.build_params(obj, dm);
            end
        end

        function obj = add_nodes(obj, nm, dm)
            %% each element adds its nodes
            for k = 1:length(obj.elements)
                if obj.elements{k}.nn   %% element has nodes to create
                    obj.elements{k}.add_nodes(obj, dm);
                end
            end

            %% if network has its own nodes to add, do that now
            if obj.nn
                add_nodes@mp.nm_element(obj, nm, dm);
            end

            %% add voltage variables for each node
            obj.add_vvars(obj, dm);
        end

        function obj = add_states(obj, nm, dm)
            %% each element adds its states
            for k = 1:length(obj.elements)
                if obj.elements{k}.nz   %% element has states to create
                    obj.elements{k}.add_states(obj, dm);
                end
            end

%             %% if network has its own states to add, do that now
%             %% (this doesn't work because nz includes states from all elements)
%             if obj.nz
%                 add_states@mp.nm_element(obj, nm, dm);
%             end

            %% add state variables for each node
            obj.add_zvars(obj, dm);
        end

        function obj = build_params(obj, nm, dm)
            %% each element builds parameters, aggregate incidence matrices
            C = {};
            D = {};
            for k = 1:length(obj.elements)
                nme = obj.elements{k};
                nme.build_params(obj, dm);
                C = horzcat(C, {nme.C});
                D = horzcat(D, {nme.D});

                %% add ports (to track indexing of all network ports)
                if nme.np > 1
                    obj.init_indexed_name('port', nme.name, {nme.np});
                    for p = 1:nme.np
                        obj.add_port(nme.name, {p}, nme.nk);
                    end
                elseif nme.np
                    obj.add_port(nme.name, nme.nk);
                end
            end
            obj.C = horzcat(C{:});
            obj.D = horzcat(D{:});
        end

        function M = stack_matrix_params(obj, name, vnotz)
            nn = obj.getN('node');
            if vnotz
                nc = nn;
            else
                nc = obj.nz;
            end
            ii = {};
            jj = {};
            ss = {};
            last_i = 0;
            last_j = 0;
            for k = 1:length(obj.elements)
                nme = obj.elements{k};
                Mk = nme.(name);
                if ~isempty(Mk)
                    [i, j, s] = find(Mk);
                    ii = horzcat(ii, i + last_i);
                    jj = horzcat(jj, j + last_j);
                    ss = horzcat(ss, s);
                end
                m = nme.nk * nme.np;        %% total number of ports for class
                if vnotz
                    n = m;
                else
                    n = nme.nk * nme.nz;    %% total number of states for class
                end
                last_i = last_i + m;
                last_j = last_j + n;
            end

            M = sparse(vertcat(ii{:}), vertcat(jj{:}), vertcat(ss{:}), last_i, last_j);
        end

        function v = stack_vector_params(obj, name)
            nn = obj.getN('node');
            vv = {};

            for k = 1:length(obj.elements)
                nme = obj.elements{k};
                vk = nme.(name);
                if isempty(vk)
                    vk = zeros(nme.nk * nme.np, 1);
                end
                vv = horzcat(vv, vk);
            end
            v = vertcat(vv{:});
        end

        function obj = add_vvars(obj, nm, dm, idx)
            for k = 1:length(obj.node.order)
                nme = obj.elements.(obj.node.order(k).name);
                nme.add_vvars(obj, dm, obj.node.order(k).idx);
            end
            for vtype = obj.model_vvars
                obj.nv = obj.nv + obj.getN(vtype{1});
            end
        end

        function obj = add_zvars(obj, nm, dm, idx)
            for k = 1:length(obj.state.order)
                nme = obj.elements.(obj.state.order(k).name);
                nme.add_zvars(obj, dm, obj.state.order(k).idx);
            end
        end

        %%-----  mp_idx_manager methods  -----
        function obj = def_set_types(obj)
            obj.set_types = struct(...
                    'node', 'NODES', ...
                    'state', 'STATES', ...
                    'port', 'PORTS' ...
                );
        end

        function obj = init_set_types(obj)
            %% call parent to create base data structures for each type
            init_set_types@mp_idx_manager(obj);

            %% finish initializing data structures for each type
            es = struct();                          %% empty struct
                                                    %% variable types
            for vtype = horzcat(obj.model_vvars(), obj.model_zvars())
                assert(isfield(obj.set_types, vtype{1}), ...
                    'var type ''%'' is missing from def_set_types()', vtype{1});
                obj.(vtype{1}).data = struct( ...
                    'v0', es, ...
                    'vl', es, ...
                    'vu', es, ...
                    'vt', es );
            end
        end

        function display(obj)
            %DISPLAY  Displays the object.
            %   Called when semicolon is omitted at the command-line. Displays the details
            %   of the variables, constraints, costs included in the model.
            %
            %   See also MP_IDX_MANAGER.

            %% Due to a bug related to inheritance in constructors in
            %% Octave 5.2 and earlier (https://savannah.gnu.org/bugs/?52614),
            %% init_set_types() cannot be called directly in the
            %% MP_IDX_MANAGER constructor, as desired.
            %%
            %% WORKAROUND:  Initialize MP_IDX_MANAGER fields here, if needed,
            %%              after object construction, but before object use.
            if isempty(obj.node)        %% only if not already initialized
                obj.init_set_types();
            end

            %% base element info
            display@mp.nm_element(obj)
            fprintf('\n');

            %% nodes and states
            obj.display_set('node');
            obj.display_set('port');
            obj.display_set('state');

            %% variables
            vvars = obj.model_vvars();
            zvars = obj.model_zvars();
            for k = 1:length(vvars)
                obj.display_set(vvars{k});
            end
            for k = 1:length(zvars)
                obj.display_set(zvars{k});
            end

            %% elements
            model_params = obj.model_params();
            fprintf('ELEMENTS\n')
            fprintf('========\n')
            fprintf('  name                  N      np    nz    class, param(m,n))\n');
            fprintf(' ----------------   --------  ----  ----  --------------------\n');
            for k = 1:length(obj.elements)
                nme = obj.elements{k};
                fprintf('  %-13s %11d %5d %5d    %s', nme.name, nme.nk, nme.np, nme.nz, class(nme));

                for j = 1:length(model_params)
                    pn = model_params{j};   %% parameter name
                    if ~isempty(nme.(pn))
                        [m, n] = size(nme.(pn));
                        fprintf(', %s(%d,%d)', pn, m, n);
                    end
                end
                fprintf('\n');
            %     nme
            end

            %% user data
            fields = fieldnames(obj.userdata);
            if ~isempty(fields)
                fprintf('\nUSER DATA\n')
                fprintf('=========\n')
                fprintf('  name                         size       class\n');
                fprintf(' ------------------------   -----------  --------------------\n');
                for k = 1:length(fields)
                    f = obj.userdata.(fields{k});
                    [m, n] = size(f);
                    fprintf('  %-24s  %5dx%-5d   %s\n', fields{k}, m, n, class(f));
                end
            end
        end

        function add_node(obj, name, idx, N)
            %   obj.add_node(name, N)
            %   obj.add_node(name, idx_list, N)
            if ~iscell(idx)
                N = idx;
                idx = {};
            end

            %% add the named node set
            obj.add_named_set('node', name, idx, N);
        end

        function add_port(obj, name, idx, N)
            %   obj.add_port(name, N)
            %   obj.add_port(name, idx_list, N)
            if ~iscell(idx)
                N = idx;
                idx = {};
            end

            %% add the named state set
            obj.add_named_set('port', name, idx, N);
        end

        function add_state(obj, name, idx, N)
            %   obj.add_state(name, N)
            %   obj.add_state(name, idx_list, N)
            if ~iscell(idx)
                N = idx;
                idx = {};
            end

            %% add the named state set
            obj.add_named_set('state', name, idx, N);
        end

        function s = set_type_idx_map(obj, set_type, idxs, dm, group_by_name)
            %%  s = obj.set_type_idx_map(set_type, idxs)
            %%  s = obj.set_type_idx_map(set_type, idxs, dm)
            %%  s = obj.set_type_idx_map(set_type, idxs, dm, group_by_name)
            %%  output struct s has fields:
            %%      name - name of named set
            %%      idx - optional cell array of indices of named set
            %%      i - index of element in named set
            %%      e - external index (i.e. corresp. row in data model)
            %%      ID - external ID (i.e. corresp. element ID in data model)
            %%      j - (only for group_by_name == 1), corresponding index of
            %%          set type (i.e. node, port or state)

            if nargin < 5
                group_by_name = 0;
                if nargin < 4
                    dm = [];
                end
            end

            s = set_type_idx_map@mp_idx_manager(obj, ...
                                            set_type, idxs, group_by_name);
            n = length(s(:));
            if nargin > 3 && ~isempty(dm)
                for k = 1:n
                    dme = dm.elements.(s(k).name);
                    s(k).e = zeros(size(s(k).i));
                    s(k).e(:) = dme.on(s(k).i);         %% external index
                    s(k).ID = s(k).e;
                    if ~isempty(dme.ID)
                        s(k).ID(:) = dme.ID(s(k).e);    %% external ID
                    end
                end
            end
        end

        function label = set_type_label(obj, set_type, idxs, dm)
            %%  label = obj.set_type_label(set_type, idxs)
            %%  label = obj.set_type_label(set_type, idxs, dm)
            label = cell(size(idxs));
            ID = cell(size(idxs));
            if nargin > 3
                s = obj.set_type_idx_map(set_type, idxs, dm);
                [ID{1:length(idxs(:))}] = deal(s.ID);
            else
                s = obj.set_type_idx_map(set_type, idxs);
                [ID{1:length(idxs(:))}] = deal(s.i);
            end
            for k = 1:length(idxs(:))
                if isempty(s(k).idx)
                    label{k} = sprintf('%s %d', s(k).name, ID{k});
                else
                    if length(s(k).idx) <= 1
                        idxstr = sprintf('%d', s(k).idx{1});
                    else
                        idxstr = [sprintf('%d', s(k).idx{1}) sprintf(',%d', s(k).idx{2:end})];
                    end
                    label{k} = sprintf('%s(%s) %d', s(k).name, idxstr, ID{k});
                end
            end
            if isscalar(idxs)               %% return scalar
                label = label{1};
            end
        end

        function add_var(obj, vtype, name, idx, varargin)
            %   obj.add_var(vtype, name, N, v0, vl, vu, vt)
            %   obj.add_var(vtype, name, N, v0, vl, vu)
            %   obj.add_var(vtype, name, N, v0, vl)
            %   obj.add_var(vtype, name, N, v0)
            %   obj.add_var(vtype, name, N)
            %   obj.add_var(vtype, name, idx_list, N, v0, vl, vu, vt)
            %   obj.add_var(vtype, name, idx_list, N, v0, vl, vu)
            %   obj.add_var(vtype, name, idx_list, N, v0, vl)
            %   obj.add_var(vtype, name, idx_list, N, v0)
            %   obj.add_var(vtype, name, idx_list, N)
            if iscell(idx)
                N = varargin{1};
                args = varargin(2:end);
            else
                N = idx;
                idx = {};
                args = varargin;
            end
            nargs = length(args);

            v0 = []; vl = []; vu = []; vt = [];
            if nargs >= 1
                v0 = args{1};
                if N > 1 && length(v0) == 1         %% expand from scalar as needed
                    v0 = v0 * ones(N, 1);
                end
                if nargs >= 2
                    vl = args{2};
                    if N > 1 && length(vl) == 1     %% expand from scalar as needed
                        vl = vl * ones(N, 1);
                    end
                    if nargs >= 3
                        vu = args{3};
                        if N > 1 && length(vu) == 1 %% expand from scalar as needed
                            vu = vu * ones(N, 1);
                        end
                        if nargs >= 4
                            vt = args{4};
                        end
                    end
                end
            end
            if isempty(v0)
                v0 = zeros(N, 1);   %% init to zero by default
            end
            if isempty(vl)
                vl = -Inf(N, 1);    %% unbounded below by default
            end
            if isempty(vu)
                vu = Inf(N, 1);     %% unbounded above by default
            end
            if isempty(vt) && N > 0
                vt = 'C';           %% all continuous by default
            end

            %% add the named variable set
            obj.add_named_set(vtype, name, idx, N);

            %% add type-specific data for var (v0, vl, vu, vt)
            d = obj.(vtype).data;
            obj.(vtype).data = [];
            if isempty(idx)
                d.v0.(name) = v0;       %% initial value
                d.vl.(name) = vl;       %% lower bound
                d.vu.(name) = vu;       %% upper bound
                d.vt.(name) = vt;       %% variable type
            else
                %% calls to substruct() are relatively expensive, so we
                %% pre-build the struct for addressing cell array fields
                %% sc = substruct('.', name, '{}', idx);
                sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
                d.v0 = subsasgn(d.v0, sc, v0);      %% initial value
                d.vl = subsasgn(d.vl, sc, vl);      %% lower bound
                d.vu = subsasgn(d.vu, sc, vu);      %% upper bound
                d.vt = subsasgn(d.vt, sc, vt);      %% variable type
            end
            obj.(vtype).data = d;
        end

        function [v0, vl, vu, vt] = params_var(obj, vtype, name, idx)
            %PARAMS_VAR  Returns initial value, lower bound and upper bound for opt variables.
            %   [V0, VL, VU] = OBJ.PARAMS_VAR(VTYPE)
            %   [V0, VL, VU] = OBJ.PARAMS_VAR(VTYPE, NAME)
            %   [V0, VL, VU] = OBJ.PARAMS_VAR(VTYPE, NAME, IDX_LIST)
            %   [V0, VL, VU, VT] = OBJ.PARAMS_VAR(...)
            %   Returns the initial value V0, lower bound VL and upper bound VU for
            %   the full optimization variable vector, or for a specific named or named
            %   and indexed variable set. Optionally also returns a corresponding char
            %   vector VT of variable types, where 'C', 'I' and 'B' represent continuous
            %   integer and binary variables, respectively.
            %
            %   Examples:
            %       [x0, xmin, xmax] = obj.params_var();
            %       [pg0, pg_lb, pg_ub] = obj.params_var('pg');
            %       [zij0, zij_lb, zij_ub, ztype] = obj.params_var('z', {i, j});
            %
            %   See also OPT_MODEL/PARAMS_VAR.

            if nargout > 3
                have_vt = 1;
            else
                have_vt = 0;
            end
            var = obj.(vtype);
            if nargin < 3
                v0 = []; vl = []; vu = []; vt = char([]);
                %% calls to substruct() are relatively expensive, so we pre-build the
                %% structs for addressing cell and numeric array fields, updating only
                %% the subscripts before use
                sc = struct('type', {'.', '{}'}, 'subs', {'', 1});  %% cell array field
                for k = 1:var.NS
                    name = var.order(k).name;
                    idx = var.order(k).idx;
                    if isempty(idx)
                        v0 = [ v0; var.data.v0.(name) ];
                        vl = [ vl; var.data.vl.(name) ];
                        vu = [ vu; var.data.vu.(name) ];
                        if have_vt
                            N = var.idx.N.(name);
                            vt0 = var.data.vt.(name);
                            if isscalar(vt0) && N > 1
                                vt = [ vt char(vt0 * ones(1, N)) ];
                            else
                                vt = [ vt vt0 ];
                            end
                        end
                    else
                        % (calls to substruct() are relatively expensive ...
                        % sc = substruct('.', name, '{}', idx);
                        % ... so replace it with these more efficient lines)
                        sc(1).subs = name;
                        sc(2).subs = idx;
                        v0 = [ v0; subsref(var.data.v0, sc) ];
                        vl = [ vl; subsref(var.data.vl, sc) ];
                        vu = [ vu; subsref(var.data.vu, sc) ];
                        if have_vt
                            % (calls to substruct() are relatively expensive ...
                            % sn = substruct('.', name, '()', idx);
                            % ... so replace it with these more efficient lines)
                            sn = sc; sn(2).type = '()';
                            N = subsref(var.idx.N, sn);
                            vt0 = subsref(var.data.vt, sc);
                            if isscalar(vt0) && N > 1
                                vt = [ vt char(vt0 * ones(1, N)) ];
                            else
                                if ~isempty(vt0)
                                    vt = [ vt vt0 ];
                                end
                            end
                        end
                    end
                end
            else
                if isfield(var.idx.N, name)
                    if nargin < 4 || isempty(idx)
                        v0 = var.data.v0.(name);
                        vl = var.data.vl.(name);
                        vu = var.data.vu.(name);
                        if have_vt
                            N = var.idx.N.(name);
                            vt0 = var.data.vt.(name);
                            if isscalar(vt0) && N > 1
                                vt = char(vt0 * ones(1, N));
                            else
                                vt = vt0;
                            end
                        end
                    else
                        % (calls to substruct() are relatively expensive ...
                        % sc = substruct('.', name, '{}', idx);
                        % ... so replace it with these more efficient lines)
                        sc = struct('type', {'.', '{}'}, 'subs', {name, idx});
                        v0 = subsref(var.data.v0, sc);
                        vl = subsref(var.data.vl, sc);
                        vu = subsref(var.data.vu, sc);
                        if have_vt
                            % (calls to substruct() are relatively expensive ...
                            % sn = substruct('.', name, '()', idx);
                            % ... so replace it with these more efficient lines)
                            sn = sc; sn(2).type = '()';
                            N = subsref(var.idx.N, sn);
                            vt0 = subsref(var.data.vt, sc);
                            if isscalar(vt0) && N > 1
                                vt = char(vt0 * ones(1, N));
                            else
                                vt = vt0;
                            end
                        end
                    end
                else
                    v0 = [];
                    vl = [];
                    vu = [];
                    if have_vt
                        vt = [];
                    end
                end
            end
        end

        function varargout = get_node_idx(obj, name)
            %% [i1 iN] = obj.get_node_idx(name)
            %% nidx = obj.get_node_idx(name), where
            %%      nidx = [i1:iN]' or
            %%      nidx = {[i1(1):iN(1)]', ..., [i1(n):iN(n)]'}
            [varargout{1:nargout}] = obj.get_set_type_idx('node', name);
        end

        function varargout = get_port_idx(obj, name)
            %% [i1 iN] = obj.get_port_idx(name)
            %% pidx = obj.get_port_idx(name), where
            %%      pidx = [i1:iN]' or
            %%      pidx = {[i1(1):iN(1)]', ..., [i1(n):iN(n)]'}
            [varargout{1:nargout}] = obj.get_set_type_idx('port', name);
        end

        function varargout = get_state_idx(obj, name)
            %% [i1 iN] = obj.get_state_idx(name)
            %% sidx = obj.get_state_idx(name), where
            %%      sidx = [i1:iN]' or
            %%      sidx = {[i1(1):iN(1)]', ..., [i1(n):iN(n)]'}
            [varargout{1:nargout}] = obj.get_set_type_idx('state', name);
        end

        function [ref, pv, pq, by_elm] = node_types(obj, nm, dm, idx, skip_ensure_ref)
            %%          ntv          = obj.node_types(nm, dm)
            %%         [ntv, by_elm] = obj.node_types(nm, dm)
            %% [ref, pv, pq]         = obj.node_types(nm, dm)
            %% [ref, pv, pq, by_elm] = obj.node_types(nm, dm)
            %%  where by_elm(k) is a struct for the k-th node-creating element
            %%  type, with fields:
            %%      'name' - name of corresponding node-creating element type
            %%      'ntv' - node type vector (if by_elm is 2nd output arg)
            %%      'ref'/'pv'/'pq' - index vectors into elements of
            %%          corresponding node-creating element type  (if by_elm
            %%          is 4th output arg)
            if nargin < 5
                skip_ensure_ref = 0;
            end

            %% empty cell array for initialization
            c0 = cell(length(obj.node.order), 1);

            %% get node types from each node-creating NME
            if nargout > 2
                %% initialize outputs
                [rr, vv, qq] = deal(c0);
                if nargout > 3
                    by_elm = struct('name', c0, 'ref', c0, 'pv', c0, 'pq', c0);
                end

                %% get node types for each node-creating element type
                for k = 1:obj.node.NS
                    name = obj.node.order(k).name;
                    idx = obj.node.order(k).idx;
                    if isempty(idx)
                        i0 = obj.node.idx.i1.(name) - 1;
                    else
                        sn = struct('type', {'.', '()'}, 'subs', {name, idx});
                        i0 = subsref(obj.node.idx.i1, sn) - 1;
                    end
                    nme = obj.elements.(name);
                    [rr{k}, vv{k}, qq{k}] = nme.node_types(obj, dm, idx);
                    if nargout > 3
                        by_elm(k).name = name;
                        by_elm(k).ref = rr{k};
                        by_elm(k).pv  = vv{k};
                        by_elm(k).pq  = qq{k};
                    end
                    [rr{k}, vv{k}, qq{k}] = deal(rr{k}+i0, vv{k}+i0, qq{k}+i0);
                end
                [ref, pv, pq] = ...     %% concatenate into single vectors
                    deal(vertcat(rr{:}), vertcat(vv{:}), vertcat(qq{:}));
                if ~skip_ensure_ref
                    [ref, pv, pq] = obj.ensure_ref_node(dm, ref, pv, pq);
                end
            else
                %% initialize outputs
                tt = c0;
                if nargout > 3
                    by_elm = struct('name', c0, 'ntv', c0);
                end

                %% get node types for each node-creating element type
                for k = 1:length(obj.node.order)
                    name = obj.node.order(k).name;
                    idx = obj.node.order(k).idx;
                    nme = obj.elements.(name);
                    tt{k} = nme.node_types(obj, dm, idx);
                    if nargout > 1
                        by_elm(k).name = name;
                        by_elm(k).ntv = tt{k};
                    end
                end
                ref = vertcat(tt{:});       %% concatenate into a single vector
                if ~skip_ensure_ref
                    ref = obj.ensure_ref_node(dm, ref);
                end
            end
        end

        function [ref, pv, pq] = ensure_ref_node(obj, dm, ref, pv, pq)
            %% [ref, pv, pq] = obj.ensure_ref_node(dm, ref, pv, pq)
            %% ntv = obj.ensure_ref_node(dm, ntv)
            if nargout > 1      %% ref, pv, pq
                if isempty(ref)
                    if isempty(pv)
                        error('mp.net_model/ensure_ref_node: must have at least one REF or PV node');
                    end
                    obj.set_node_type_ref(dm, pv(1));
                    ref = pv(1);
                    pv(1) = [];
                end
            else                %% ntv
                ntv = ref;
                ref = find(ntv == mp.NODE_TYPE.REF);    %% ref node indices
                if isempty(ref)
                    pv = find(ntv == mp.NODE_TYPE.PV);  %% PV node indices
                    if isempty(pv)
                        error('mp.net_model/ensure_ref_node: must have at least one REF or PV node');
                    end
                    obj.set_node_type_ref(dm, pv(1));

                    %% update node type vector, skipping ensure ref node step
                    ntv = obj.node_types(obj, dm, {}, 1);
                end
                ref = ntv;
            end
        end

        function set_node_type_ref(obj, dm, idx)
            %% obj.set_node_type_ref(dm, idx)
            %% idx is internal nme index
            s = obj.set_type_idx_map('node', idx, dm, 1);
            for k = 1:length(s)
                nme = obj.elements.(s(k).name);
                nme.set_node_type_ref(obj, dm, s(k).i);
            end
        end

        function set_node_type_pv(obj, dm, idx)
            %% obj.set_node_type_pv(dm, idx)
            %% idx is internal nme index
            s = obj.set_type_idx_map('node', idx, dm, 1);
            for k = 1:length(s)
                nme = obj.elements.(s(k).name);
                nme.set_node_type_pv(obj, dm, s(k).i);
            end
        end

        function set_node_type_pq(obj, dm, idx)
            %% obj.set_node_type_pq(dm, idx)
            %% idx is internal nme index
            s = obj.set_type_idx_map('node', idx, dm, 1);
            for k = 1:length(s)
                nme = obj.elements.(s(k).name);
                nme.set_node_type_pq(obj, dm, s(k).i);
            end
        end
    end     %% methods

    methods (Access = private)
        function [i1, iN] = get_set_type_idx(obj, node_state, name)
            %% [i1 iN] = obj.get_set_type_idx(node_state, name)
            %% stidx = obj.get_set_type_idx(node_state, name), where
            %%      stidx = [i1:iN]' or
            %%      stidx = {[i1(1):iN(1)]', ..., [i1(n):iN(n)]'}
            idx = obj.get_idx(node_state);
            i1 = idx.i1.(name);
            iN = idx.iN.(name);
            if nargout == 1
                N = length(i1);
                if N == 1
                    i1 = (i1:iN)';
                else
                    t = cell(1, N);
                    for k = 1:N
                        t{k} = (i1(k):iN(k))';
                    end
                    i1 = t;
                end
            end
        end
    end     %% methods (private)
end         %% classdef
