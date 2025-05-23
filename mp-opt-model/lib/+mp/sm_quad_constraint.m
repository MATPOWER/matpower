classdef sm_quad_constraint < mp.set_manager_opt_model
% mp.sm_quad_constraint -  MP Set Manager class for quadratic constraints.
% ::
%
%   qcn = mp.sm_quad_constraint()
%   qcn = mp.sm_quad_constraint(label)
%
% MP Set Manager class for quadratic constraints of the form
%
% .. math::
%    :label: eq_qcn_form
%
%    \li \le \frac{1}{2}\trans{\x} \QQ_i \x + \b_i \x \le \ui, \ \ \ 
%    i = 1,2,..., n_q
%
% Manages quadratic constraint sets and their indexing.
%
% The parameters defining the set of constraints are
% :math:`\QQ_1, \dots, \QQ_{n_q}, \Bb, \l, \u`, where :math:`\b_i` in
% :eq:`eq_qcn_form` is row :math:`i` of matrix :math:`\Bb`, :math:`\li`
% and :math:`\ui` are the :math:`i`-th elements of vectors :math:`\l` and
% :math:`\u`, respectively, and :math:`n_q` is the number of quadratic
% constraints.
%
% Let :math:`\{\rmat{A}\}_{\times n}` denote the set :math:`\{\rmat{A},
% \rmat{A}, \dots, \rmat{A}\}` where :math:`\rmat{A}` is repeated :math:`n`
% times, let :math:`\mathcal{A} = \{\rmat{A}_i\}_{i=1}^n` be a set of :math:`n`
% matrices or vectors indexed by :math:`i`, and let
% :math:`\textrm{diag}(\rmat{A}) = \trans{\left[a_{11} \ a_{22} \ \dots \
% a_{nn}\right]}` denote the matrix-to-vector diagonal operator. If we also let
% :math:`\diag{\mathcal{A}}` denote a block diagonal matrix with the set
% :math:`\mathcal{A}` of matrices or vectors on the block diagonal, then
% :eq:`eq_qcn_form` can be expressed in matrix form as:
%
% .. math::
%    :label: eq_qcn_matform
%
%    \l &\le \frac{1}{2}\textrm{diag}\left(\trans{\Xblk}
%       \Qblk \Xblk\right) + \Bb \x \le \u \\
%    \l &\le \frac{1}{2}\textrm{diag}\left(\trans{\diag{\{\x\}_{\times n_q}}}
%       \diag{\{\QQ_i\}_{i=1}^{n_q}} \diag{\{\x\}_{\times n_q}}\right) + \Bb \x
%       \le \u,
%
% where
%
% .. math::
%    :label: eq_qcn_matdef
%
%    \Xblk = \diag{\{\x\}_{\times n_q}} \\
%    \Qblk = \diag{\{\QQ_i\}_{i=1}^{n_q}}.
%
% By convention, ``qcn`` is the variable name used for mp.sm_quad_constraint
% objects.
%
% mp.sm_quad_constraint Properties:
%   * cache - struct for caching aggregated parameters for the set
%
% mp.sm_quad_constraint Methods:
%   * sm_quad_constraint - constructor
%   * add - add a subset of quadratic constraints, with parameters
%     :math:`\QQ_1, \dots, \QQ_{n_q}, \Bb, \l, \u`
%   * params - build and return quadratic constraint parameters
%     :math:`\QQ_1, \dots, \QQ_{n_q}, \Bb, \l, \u`
%   * set_params - modify quadratic constraint parameter data
%   * eval - evaluate individual or full set of quadratic constraints
%   * display_soln - display solution values for quadratic constraints
%   * get_soln - fetch solution values for specific named/indexed subsets
%   * parse_soln - parse solution for quadratic constraints
%   * blkprod2vertcat - compute all products of two block diagonal matrices
%
% See also mp.set_manager, mp.set_manager_opt_model.

%   MP-Opt-Model
%   Copyright (c) 2019-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

    properties
        % struct for caching aggregated parameters for quadratic constraints
        cache = [];
    end     %% properties

    methods
        function obj = sm_quad_constraint(varargin)
            % Constructor.
            % ::
            %
            %   qcn = mp.sm_quad_constraint(label)

            es = struct();  %% empty struct
            obj@mp.set_manager_opt_model(varargin{:});
            obj.data = struct( ...
                'Q', es, ...
                'B', es, ...
                'l', es, ...
                'u', es, ...
                'vs', es );
        end

        function obj = add(obj, var, name, idx, varargin)
            % Add a subset of quadratic constraints, with parameters :math:`\QQ_1, \dots, \QQ_{n_q}, \Bb, \l, \u`.
            % ::
            %
            %   qcn.add(var, name, Q, B, l, u);
            %   qcn.add(var, name, Q, B, l, u, vs);
            %
            %   qcn.add(var, name, idx_list, Q, B, l, u);
            %   qcn.add(var, name, idx_list, Q, B, l, u, vs);
            %
            % Add a named, and possibly indexed, subset of :math:`n_q` quadratic
            % constraints to the set, of the form in :eq:`eq_qcn_form` and
            % :eq:`eq_qcn_matform`, where :math:`\x` is a vector made up of the
            % variables specified in the optional ``vs`` *(in the order given)*.
            % This allows the set of :math:`\QQ_i` matrices and :math:`\Bb`
            % matrix to be defined in terms of only the relevant variables
            % without the need to manually create a lot of properly located zero
            % rows and/or columns.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   name (char array) : name of subset/block of constraints to add
            %   idx_list (cell array) : *(optional)* index list for subset/block
            %       of constraints to add (for an indexed subset)
            %   Q (cell vector): :math:`n_q \times 1` cell array of sparse
            %       quadratic coefficient matrices :math:`\QQ_i`, where the
            %       size of each element must be consistent with :math:`\x`.
            %   B (double) : possibly sparse matrix :math:`\Bb` of linear
            %       coefficients, with row vector :math:`\b_i` corresponding
            %       to quadratic constraint :math:`i`.
            %   l (double) : *(optional, default =* ``-Inf`` *)* constraint
            %       left-hand side vector :math:`\l`, or scalar which is
            %       expanded to a vector
            %   u (double) : *(optional, default =* ``Inf`` *)* constraint
            %       right-hand side vector :math:`\u`, or scalar which is
            %       expanded to a vector
            %   vs (cell or struct array) : *(optional, default* ``{}`` *)*
            %       variable set defining vector :math:`\x` for this
            %       constraint subset; can be either a cell array of names of
            %       variable subsets, or a struct array of ``name``, ``idx``
            %       pairs of indexed named subsets of variables; order of
            %       ``vs`` determines order of blocks in :math:`\x`; if
            %       empty, :math:`\x` is assumed to be the full variable vector
            %
            % Examples::
            %
            %   qcn.add(var, 'my_quad', Q1, B1, l1, u1, {'my_var1', 'my_var2'});
            %
            %   qcn.init_indexed_name('my_set', {3, 2})
            %   for i = 1:3
            %       for j = 1:2
            %           qcn.add(var, 'my_set', {i, j}, Q{i,j}, B{i,j}, l{i,j}, ...);
            %       end
            %   end
            %
            % See also params, set_params, eval.

            %% set up default args
            if isfield(obj.idx.N, name)     %% indexed named set
                Q = varargin{1};
                args = varargin(2:end);
            else                            %% simple named set
                Q = idx;
                idx = {};
                args = varargin;
            end
            nargs = length(args);

            %% prepare data
            B = []; l = []; u = []; vs = {};
            if nargs >= 1
                B = args{1};
                if nargs >= 2
                    l = args{2};
                    if nargs >= 3
                        u = args{3};
                        if nargs >= 4
                            vs = args{4};
                        end
                    end
                end
            end

            %% convert varsets from cell to struct array if necessary
            vs = mp.sm_variable.varsets_cell2struct(vs);
            nv = var.varsets_len(vs);   %% number of variables

            %% get sizes
            [MQ, NQ] = size(Q);
            [MB, NB] = size(B);
            if ~isempty(Q)
                N = MQ;
            elseif ~isempty(B)
                N = MB;
            elseif nv
                error('mp.sm_quad_constraint.add: Q and B cannot both be empty');
            end

            %% check parameters and assign defaults
            if ~isempty(Q)
                if NQ ~= 1
                    error('mp.sm_quad_constraint.add: Q (%d x %d) must be a column cell array (or empty)', MQ, NQ);
                end
                for k = 1:N
                    [MQi, NQi] = size(Q{k});
                    if MQi ~= nv || NQi ~= nv
                        error('mp.sm_quad_constraint.add: Q{%d} (%d x %d) be a square matrix matching the number of variables (%d)', k, MQi, NQi, nv);
                    end
                end
            end
            if ~isempty(B)
                if MB ~= N
                    error('mp.sm_quad_constraint.add: number of rows in Q (%d) and B (%d) must match', N, MB);
                end
                if NB ~= nv
                    error('mp.sm_quad_constraint.add: number of columns in B (%d) must match the number of variables (%d)', NB, nv);
                end
            else
                B = sparse(N, nv);
            end
            if isempty(l)
                l = -Inf(N, 1);
            elseif N > 1 && numel(l) == 1
                l = l*ones(N, 1);
            elseif numel(l) ~= N
                error('mp.sm_quad_constraint.add: l (%d x %d) must (%d x 1), a scalar, or empty', size(l, 1), size(l, 2), N);
            end
            if isempty(u)
                u = Inf(N, 1);
            elseif N > 1 && numel(u) == 1
                u = u*ones(N, 1);
            elseif numel(u) ~= N
                error('mp.sm_quad_constraint.add: u (%d x %d) must (%d x 1), a scalar, or empty', size(u, 1), size(u, 2), N);
            end

            %% call parent to handle standard indexing
            if isempty(idx)
                add@mp.set_manager_opt_model(obj, name, N, args{:});
            else
                add@mp.set_manager_opt_model(obj, name, idx, N, args{:});
            end

            %% assign data (wgv: this replaces the old add_named_set method)
            if isempty(idx)
                obj.data.Q.(name)  = Q;
                obj.data.B.(name)  = B;
                obj.data.l.(name)  = l;
                obj.data.u.(name)  = u;
                obj.data.vs.(name) = vs;
            else
                % (calls to substruct() are relatively expensive ...
                % s = substruct('.', name, '{}', idx);
                % ... so replace it with these more efficient lines)
                sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
                obj.data.Q  = subsasgn(obj.data.Q, sc, Q);
                obj.data.B  = subsasgn(obj.data.B, sc, B);
                obj.data.l  = subsasgn(obj.data.l, sc, l);
                obj.data.u  = subsasgn(obj.data.u, sc, u);
                obj.data.vs = subsasgn(obj.data.vs, sc, vs);
            end
            if ~isempty(obj.cache)  %% clear cache of aggregated params
                obj.cache = [];
            end
        end

        function [Qblk, B, l, u, vs, i1, iN] = params(obj, var, name, idx, isblk)
            % Build and return quadratic constraint parameters :math:`\Qblk, \Bb, \l, \u`.
            % ::
            %
            %   [Qblk, B, l, u] = qcn.params(var)
            %   [Qblk, B, l, u] = qcn.params(var, name)
            %   [Qblk, B, l, u] = qcn.params(var, name, idx_list)
            %   [Qblk, B, l, u, vs] = qcn.params(...)
            %   [Qblk, B, l, u, vs, i1, iN] = qcn.params(name ...)
            %
            % With no input parameters, it assembles and returns the parameters
            % for the aggregate quadratic constraints from all quadratic
            % constraint sets added using add(). The values of these parameters
            % are cached for subsequent calls. The parameters are
            % :math:`\Qblk`, :math:`\Bb`, :math:`\l`, and :math:`\u`, where the
            % quadratic constraint is of the form in :eq:`eq_qcn_matform`.
            %
            % If ``isblk`` is 1, ``Qblk`` is the block diagonal matrix
            % :math:`\Qblk` of :eq:`eq_qcn_matdef`. If ``isblk`` is 0, ``Qblk``
            % is a vertical cell array containing the individual quadratic
            % coefficient matrices :math:`\QQ_i` for all quadratic constraint
            % sets. If a name or name and index list are provided, then it
            % simply returns the parameters for the corresponding named set, as
            % a block diagonal matrix or cell array, depending on the value of
            % ``isblk``. It can also optionally return the variable sets used by
            % this constraint set (the size of :math:`\QQ_i` and :math:`\Bb`
            % will be consistent with this variable set), and the starting and
            % ending row indices of the subset within the full aggregate
            % constraint set.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   name (char array) : *(optional)* name of subset
            %   idx_list (cell array) : *(optional)* index list for subset
            %   isblk : *(optional default = 0)* determines form of ``Qblk``
            %       output, 1 for block diagonal matrix :math:`\Qblk`, 0 for
            %       individual :math:`\QQ_i` matrices stacked vertically in
            %       separate elements of a cell array
            %
            % Outputs:
            %   Qblk (double or cell array) : constraint coefficient matrix
            %       :math:`\Qblk` or cell array of individual :math:`\QQ_i`
            %       depending on the value of input ``isblk``
            %   B (double) : linear term contraint coefficient matrix
            %       :math:`\Bb`
            %   l (double) : constraint left-hand side vector :math:`\l`
            %   u (double) : constraint right-hand side vector :math:`\u`
            %   vs (struct array) : variable set, ``name``, ``idx`` pairs
            %       specifying the set of variables defining vector :math:`\x`
            %       for this constraint subset; order of ``vs`` determines
            %       order of blocks in :math:`\x`
            %   i1 (integer) : index of 1st row of specified subset in full set
            %   iN (integer) : index of last row of specified subset in full set
            %
            % Examples::
            %
            %   [Qblk, B, l, u] = qcn.params(var)
            %   [Qblk, B, l, u] = qcn.params(var, 'my_set')
            %
            % See also add.

            if nargin > 2       %% individual set
                if nargin < 5
                    isblk = 0;
                    if nargin < 4
                        idx = {};
                    end
                end
                if isempty(idx)                 %% name, no index provided
                    if numel(obj.idx.i1.(name)) == 1     %% simple named set
                        Qblk = obj.data.Q.(name);   %% cell array matrices
                        if isblk                    %% block-diagonal matrix
                            if issparse(Qblk{1})    %% sparse matrix
                                Qblk = blkdiag(Qblk{:});
                            else                    %% dense [row col var] matrix
                                % Pending ...
                            end
                        end
                        B = obj.data.B.(name);
                        l = obj.data.l.(name);
                        u = obj.data.u.(name);
                        if nargout > 4
                            vs = obj.data.vs.(name);
                            if nargout > 6
                                i1 = obj.idx.i1.(name);      %% starting row index
                                iN = obj.idx.iN.(name);      %% ending row index
                            end
                        end
                    else                                    %% indexing required
                        error('mp.sm_lin_constraint.params: quadratic constraint set ''%s'' requires an IDX_LIST arg', name);
                    end
                else                            %% indexed named set
                    % (calls to substruct() are relatively expensive ...
                    % s = substruct('.', name, '{}', idx);
                    % ... so replace it with these more efficient lines)
                    sc = struct('type', {'.', '{}'}, 'subs', {name, idx});
                    Qblk  = subsref(obj.data.Q, sc);
                    B     = subsref(obj.data.B, sc);
                    l     = subsref(obj.data.l, sc);
                    u     = subsref(obj.data.u, sc);
                    if isblk
                        Qblk = blkdiag(Qblk{:});
                    end
                    if nargout > 3
                        vs = subsref(obj.data.vs, sc);
                        if nargout > 5
                            sn = sc; sn(2).type = '()';     %% num array field
                            i1 = subsref(obj.idx.i1, sn);   %% starting row index
                            iN = subsref(obj.idx.iN, sn);   %% ending row index
                        end
                    end
                end
            else                %% aggregate
                cache = obj.cache;
                if isempty(cache)
                    % Initialize parameters for aggregate
                    nx = var.N;              %% number of full set of variables
                    nquad = obj.N;           %% number of quadratic constraints
                    Qblk = cell(nquad, 1);   %% cell array of quadratic constraints
                    B = sparse(nquad, nx);   %% matrix of linear components of quadratic constraints
                    u = Inf(nquad, 1);       %% upper bound
                    l = -u;                  %% lower bound

                    for j = 1:obj.NS   %% For each set of quadratic constraints
                        name = obj.order(j).name;
                        idx  = obj.order(j).idx;
                        [Qj, Bj, lj, uj, vsj, i1, iN] = obj.params(var, name, idx);

                        if isempty(vsj)   % full nx vars
                            Qblk(i1:iN) = Qj;
                            B(i1:iN,:) = Bj;
                        else              % subset of vars
                            n = obj.get_N(name, idx); % number of constraints in current set
                            jj = var.varsets_idx(vsj); % indices for var set
                            rowid = repmat(jj', length(jj), 1);
                            colid = repmat(jj, length(jj), 1); colid = colid(:);
                            rowid = repmat(rowid, n, 1);
                            colid = repmat(colid, n, 1);
                            vals = cellfun(@(x)(x(:)), Qj, 'UniformOutput', false);
                            Qaux = mat2cell([rowid colid cell2mat(vals)], length(jj)^2*ones(n,1));
                            Qblk(i1:iN) = cellfun(@(x)(sparse(x(:,1), x(:,2), x(:,3), nx, nx)), Qaux, 'UniformOutput', false);
                            Baux = sparse(n, length(jj));
                            Baux(:, 1:length(jj)) = Bj;
                            B(i1:iN, jj) = Baux;
                        end
                        l(i1:iN,:) = lj;
                        u(i1:iN,:) = uj;
                    end

                    %% cache aggregated parameters
                    obj.cache = struct('B', B, 'l', l, 'u', u);
                    obj.cache.Qblk = Qblk;
                else
                    Qblk = cache.Qblk;
                    B = cache.B;
                    l = cache.l;
                    u = cache.u;
                end
                if nargout > 3
                    vs = {};
                end
            end
        end

        function obj = set_params(obj, var, name, idx, params, vals)
            % Modify quadratic constraint parameter data.
            % ::
            %
            %   qcn.set_params(var, name, params, vals)
            %   qcn.set_params(var, name, idx_list, params, vals)
            %
            % This method can be used to modify parameters for an existing
            % subset of quadratic constraints.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   name (char array) : name of subset/block of quadratic
            %       constraints to modify
            %   idx_list (cell array) : *(optional)* index list for subset/block
            %       of quadratic constraints to modify (for an indexed subset)
            %   params : can be one of three options:
            %
            %       - ``'all'`` - indicates that ``vals`` is a cell array
            %         whose elements correspond to the input parameters of
            %         the add() method
            %       - name of a parameter - ``val`` is the value of that
            %         parameter
            %       - cell array of parameter names - ``vals`` is a cell array
            %         of corresponding values
            %   vals : new value or cell array of new values corresponding to
            %       ``params``
            %
            % Valid parameter names are ``Q``, ``B``, ``l``, ``u``, ``vs``.
            %
            % Examples::
            %
            %   qcn.set_params(var, 'y', {2,3}, {'l', 'u'}, {l, u});
            %   qcn.set_params(var, 'Pmis', 'all', {Q, B, l, u, vs});
            %
            % See also add, params.

            if nargin < 6
                vals = params;
                params = idx;
                idx = {};
            end

            %% create default list of parameters to update based on set type & inputs
            default_params = {'Q', 'B', 'l', 'u', 'vs'};

            %% standardize provided arguments in cell arrays params, vals
            [is_all, np, params, vals] = ...
                obj.set_params_std_args(default_params, params, vals);

            %% get current parameters
            [Q, B, l , u, vs] = obj.params(var, name, idx);
            MB0 = size(B, 1);
            MQ0 = size(Q, 1);
            if isempty(vs), vs = {vs}; end
            p = struct('Q', {Q}, 'B', B, 'l', l, 'u', u, 'vs', vs); %% current parameters
            u = struct('Q', 0, 'B', 0, 'l', 0, 'u', 0, 'vs',  0);   %% which ones to update

            %% replace with new parameters
            for j = 1:np
                p.(params{j}) = vals{j};
                u.(params{j}) = 1;
            end

            %% set missing default params for 'all'
            [MB, NB] = size(p.B);
            [MQ, NQ] = size(p.Q);
            if is_all
                u.Q = 1;            %% always update Q
                u.B = 1;            %% always update B
                u.l = 1;            %% always update l
                if np < 5
                    p.vs = {};
                    u.vs = 1;       %% update vs
                    if np < 4
                        p.u = Inf(MB, 1);
                        u.u = 1;    %% update u
                    end
                end
            end

            %% check consistency of parameters
            % Q must be a MQ x 1 cell array
            if NQ ~= 1
                error('mp.sm_quad_constraint.set_params: ''Q'' must be column vector cell array');
            end

            %% no dimension change unless 'all'
            if (MB ~= MB0 && ~is_all) || (MQ ~= MQ0 && ~is_all)
                error('mp.sm_quad_constraint.set_params: dimension change for ''%s'' not allowed except for ''all''', obj.nameidxstr(name, idx));
            end

            %% check sizes of new values of l, u
            for pn = {'l', 'u'}
                if u.(pn{1})
                    nn = length(p.(pn{1}));
                    if nn ~= MB
                        if nn == 0
                            switch pn{1}
                                case 'l'
                                    p.(pn{1}) = -Inf(MB, 0);
                                case 'u'
                                    p.(pn{1}) =  Inf(MB, 0);
                            end
                        elseif nn == 1
                            p.(pn{1}) = p.(pn{1}) * ones(MB, 1);   %% expand from scalar
                        else
                            error('mp.sm_quad_constraint.set_params: parameter ''%s'' ''%s'' should have length %d (or 1)', obj.nameidxstr(name, idx), pn{1}, NB);
                        end
                    end
                end
            end

            %% check consistency of Q and vs
            p.vs = mp.sm_variable.varsets_cell2struct(p.vs);
            nv = var.varsets_len(p.vs);     %% number of variables
            if u.Q
                if iscell(p.Q)
                    sz = cellfun(@(x)(size(x)), p.Q, 'UniformOutput', false);
                else
                    error('mp.sm_quad_constraint.set_params: parameter ''Q'' must be a %d x 1 cell array', MQ);
                end
                if sum(sum(cell2mat(sz)))/(2*MQ) ~= nv
                    error('mp.sm_quad_constraint.set_params: for ''%s'' number of columns of ''Q'' (%d) must be consistent with ''vs'' (%d)', obj.nameidxstr(name, idx), sum(sum(cell2mat(sz)))/(2*MQ) , nv);
                end
            end

            %% check consistency of B and vs
            if u.B || u.vs
                if NB ~= nv
                    error('mp.sm_quad_constraint.set_params: for ''%s'' number of columns of ''B'' (%d) must be consistent with ''vs'' (%d)', obj.nameidxstr(name, idx), MB, nv);
                end
            end

            %% assign new parameters
            if isempty(idx)     %% simple named set
                for k = 1:length(default_params)
                    pn = default_params{k};     %% param name
                    if u.(pn)   %% assign new val for this parameter
                        obj.data.(pn).(name) = p.(pn);
                    end
                end
            else                %% indexed named set
                %% calls to substruct() are relatively expensive, so we pre-build the
                %% struct for addressing cell array fields
                %% sc = substruct('.', name, '{}', idx);
                sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
                for k = 1:length(default_params)
                    pn = default_params{k};     %% param name
                    if u.(pn)   %% assign new val for this parameter
                        obj.data.(pn) = subsasgn(obj.data.(pn), sc, p.(pn));
                    end
                end
            end

            %% clear cached parameters
            obj.cache = [];

            %% update dimensions and indexing, if necessary
            dMB = MB - MB0;
            if is_all && dMB
                obj.set_params_update_dims(dMB, name, idx);
            end
        end

        function [g_u, J, g, l_g] = eval(obj, var, x, name, idx)
            % Evaluate individual or full set of quadratic constraints.
            % ::
            %
            %   g_u = qcn.eval(var, x)
            %   g_u = qcn.eval(var, x, name)
            %   g_u = qcn.eval(var, x, name, idx_list)
            %   [g_u, J] = qcn.eval(...)
            %   [g_u, J, g] = qcn.eval(...)
            %   [g_u, J, g, l_g] = qcn.eval(...)
            %
            % For a given value of the variable vector :math:`\x`, this method
            % evaluates the quadratic constraints for an individual subset, if
            % name or name and index list are provided, otherwise, for the full
            % set of constraints.
            %
            % The constraints are of the form
            %
            % .. math::
            %    :label: eq_qcn_eval_1
            %
            %    \li \le g_i(\x) \le \ui, \ \ \ i = 1,2,..., n_q
            %
            % or in matrix form for a subset or full set of constraints
            %
            % .. math::
            %    :label: eq_qcn_eval_all
            %
            %    \l \le \g(\x) \le \u
            %
            % where
            %
            % .. math::
            %    :label: eq_qcn_eval_g1
            %
            %    g_i(\x) = \frac{1}{2}\trans{\x} \QQ_i \x + \b_i \x
            %
            % and
            %
            % .. math::
            %    :label: eq_qcn_eval_g
            %
            %    \g(\x) &= \frac{1}{2}\textrm{diag}\left(\trans{\Xblk}
            %       \Qblk \Xblk\right) + \Bb \x \\
            %       &=\frac{1}{2}\textrm{diag}\left(\trans{
            %           \diag{\{\x\}_{\times n_q}}}\diag{\{\QQ_i\}_{i=1}^{n_q}}
            %           \diag{\{\x\}_{\times n_q}}\right) + \Bb \x
            %
            % Returns :math:`\g(\x) - \u`, and optionally the jacobian
            % :math:`\rmat{J}`, the contraint function :math:`\g(\x)`, and
            % :math:`\l - \g(\x)`.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   x (double) : full :math:`n_x \times 1` variable vector :math:`\x`
            %   name (char array) : name of subset/block of quadratic constraints
            %       to evaluate
            %   idx_list (cell array) : *(optional)* index list for subset/block
            %       of quadratic constraints to evaluate (for an indexed subset)
            %
            % Outputs:
            %   g_u (double) : value of :math:`\g(\x) - \u`
            %   J (double) : *(optional)* constraint Jacobian
            %       :math:`\rmat{J} = \der{\g}{\x}`
            %   g (double) : *(optional)* value of :math:`\g(\x)`
            %   l_g (double) : *(optional)* value of :math:`\l - \g(\x)`
            %
            % Examples::
            %
            %   g_u = qcn.eval(var, x)
            %   [g_u, J] = qcn.eval(var, x)
            %   [g_u, J, g] = qcn.eval(var, x)
            %   [g_u, J, g, l_g] = qcn.eval(var, x)
            %   [g_u, J, g, l_g] = qcn.eval(var, x, 'my_set')
            %   [g_u, J, g, l_g] = qcn.eval(var, x, 'my_set', {3,2})
            %
            % See also add, params.

            if obj.N
                %% collect constraint parameters
                if nargin < 4                       %% full set
                    [Qblk_aux, B, l, u, vs] = obj.params(var);
                    Qblk = blkdiag(Qblk_aux{:});
                    Nq = obj.N;
                elseif nargin < 5 || isempty(idx)   %% name, no idx provided
                    dims = size(obj.idx.i1.(name));
                    if prod(dims) == 1              %% simple named set
                        [Qblk, B, l, u, vs] = obj.params(var, name, {}, 1);
                        Nq = obj.get_N(name);
                    else
                        error('mp.sm_quad_constraint.eval: quadratic constraint set ''%s'' requires an IDX_LIST arg', name)
                    end
                else                                %% indexed named set
                    [Qblk, B, l, u, vs] = obj.params(var, name, idx, 1);
                    Nq = obj.get_N(name, idx);
                end

                %% assemble block diagonal matrix from x vector
                x = var.varsets_x(x, vs, 'vector');
                xx = mat2cell(repmat(sparse(x'), Nq, 1), ones(Nq,1));
                blkx = blkdiag(xx{:});

                %% Compute quadratic constraints
                if isempty(B)
                    gg = 1/2 * diag(blkx * Qblk * blkx');
                    g_u = gg - u;
                else
                    gg = 1/2 * diag(blkx * Qblk * blkx') + B * x;
                    g_u = gg - u;
                end

                if nargout > 1 %% Jacobian is requested
                    Qx = obj.blkprod2vertcat(blkx, Qblk, length(x));
                    J = Qx + B;
                    if nargout > 2
                        g = gg;
                        if nargout > 3
                            l_g = l - gg;
                        end
                    end
                end
            else
                g_u = [];
                if nargout > 1
                    J = [];
                    if nargout > 2
                        g = [];
                        if nargout > 3
                            l_g = [];
                        end
                    end
                end
            end
        end

        function obj = display_soln(obj, var, soln, varargin)
            % Display solution values for quadratic constraints.
            % ::
            %
            %   qcn.display_soln(var, soln)
            %   qcn.display_soln(var, soln, name)
            %   qcn.display_soln(var, soln, name, idx_list)
            %   qcn.display_soln(var, soln, fid)
            %   qcn.display_soln(var, soln, fid, name)
            %   qcn.display_soln(var, soln, fid, name, idx_list)
            %
            % Displays the solution values for all quadratic constraints
            % (default) or an individual named or named/indexed subset.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   soln (struct) : full solution struct with these fields
            %       (among others):
            %
            %           - ``eflag`` - exit flag, 1 = success, 0 or negative =
            %             solver-specific failure code
            %           - ``x`` - variable values
            %           - ``lambda`` - constraint shadow prices, struct with
            %             fields:
            %
            %               - ``eqnonlin`` - nonlinear equality constraints
            %               - ``ineqnonlin`` - nonlinear inequality constraints
            %               - ``mu_l`` - linear constraint lower bounds
            %               - ``mu_u`` - linear constraint upper bounds
            %               - ``mu_lq`` - quadratic constraint lower bounds
            %               - ``mu_uq`` - quadratic constraint upper bounds
            %               - ``lower`` - variable lower bounds
            %               - ``upper`` - variable upper bounds
            %   fid (fileID) : fileID of open file to write to (default is
            %       1 for standard output)
            %   name (char array) : *(optional)* name of individual subset
            %   idx_list (cell array) : *(optional)* indices of individual
            %       subset

            [fid, name, idx, idxs, hdr1] = ...
                obj.display_soln_std_args(varargin{:});

            if obj.N
                [Q, B, vl, vu] = obj.params(var);
                Nq = length(vl);
                Qblk = blkdiag(Q{:});
                xx = mat2cell(repmat(sparse(soln.x'), Nq, 1), ones(Nq,1));
                blkx = blkdiag(xx{:});
                v = diag(1/2 * blkx * Qblk * blkx' + B * soln.x);
                if isempty(soln.lambda)
                    mu_l = NaN(size(v));
                    mu_u = mu_l;
                else
                    mu_l = soln.lambda.mu_lq;
                    mu_u = soln.lambda.mu_uq;
                end

                %% print header rows
                hdr2 = {'   mu_lb     lb       val      ub      mu_ub', ...
                        ' -------- -------- -------- -------- --------' };
                obj.display_soln_print_headers(fid, hdr1, hdr2);

                %% print data
                none = '- ';
                for k = 1:length(idxs)
                    obj.display_soln_print_row(fid, idxs(k));

                    if isnan(mu_l(idxs(k)))
                        mu_lb = sprintf( ' ');
                    elseif abs(mu_l(idxs(k))) < obj.mu_thresh()
                        mu_lb = sprintf(none);
                    else
                        mu_lb = obj.sprintf_num(8, mu_l(idxs(k)));
                    end
                    if isnan(mu_u(idxs(k)))
                        mu_ub = sprintf( ' ');
                    elseif abs(mu_u(idxs(k))) < obj.mu_thresh()
                        mu_ub = sprintf(none);
                    else
                        mu_ub = obj.sprintf_num(8, mu_u(idxs(k)));
                    end
                    if vl(idxs(k)) < -obj.num_inf()
                        lb = sprintf(none);
                    else
                        lb = obj.sprintf_num(8, vl(idxs(k)));
                    end
                    if vu(idxs(k)) > obj.num_inf()
                        ub = sprintf(none);
                    else
                        ub = obj.sprintf_num(8, vu(idxs(k)));
                    end
                    fprintf(fid, '%9s%9s%9s%9s%9s\n', ...
                        mu_lb, lb, obj.sprintf_num(8, v(idxs(k))), ub, mu_ub);
                end

                %% print footer rows
                fprintf(fid, '%s\n', [hdr1{2} hdr2{2}]);
                fprintf(fid, '%7s %-28s%9s%9s%9s%9s%9s\n', '', 'Min', ...
                    obj.sprintf_num(8, min(mu_l)), obj.sprintf_num(8, min(vl)), ...
                    obj.sprintf_num(8, min(v)), ...
                    obj.sprintf_num(8, min(vu)), obj.sprintf_num(8, min(mu_u)));
                fprintf(fid, '%7s %-28s%9s%9s%9s%9s%9s\n', '', 'Max', ...
                    obj.sprintf_num(8, max(mu_l)), obj.sprintf_num(8, max(vl)), ...
                    obj.sprintf_num(8, max(v)), ...
                    obj.sprintf_num(8, max(vu)), obj.sprintf_num(8, max(mu_u)));
                fprintf(fid, '\n');
            end
        end

        function varargout = get_soln(obj, var, soln, varargin)
            % Fetch solution values for specific named/indexed subsets.
            % ::
            %
            %   vals = qcn.get_soln(var, soln, name)
            %   vals = qcn.get_soln(var, soln, name, idx_list)
            %   vals = qcn.get_soln(var, soln, tags, name)
            %   vals = qcn.get_soln(var, soln, tags, name, idx_list)
            %
            % Returns named/indexed quadratic constraint results for a solved
            % model, evaluated at the solution found.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   soln (struct) : full solution struct with these fields
            %       (among others):
            %
            %           - ``eflag`` - exit flag, 1 = success, 0 or negative =
            %             solver-specific failure code
            %           - ``x`` - variable values
            %           - ``lambda`` - constraint shadow prices, struct with
            %             fields:
            %
            %               - ``eqnonlin`` - nonlinear equality constraints
            %               - ``ineqnonlin`` - nonlinear inequality constraints
            %               - ``mu_l`` - linear constraint lower bounds
            %               - ``mu_u`` - linear constraint upper bounds
            %               - ``mu_lq`` - quadratic constraint lower bounds
            %               - ``mu_uq`` - quadratic constraint upper bounds
            %               - ``lower`` - variable lower bounds
            %               - ``upper`` - variable upper bounds
            %   tags (char array or cell array of char arrays) : names of
            %       desired outputs, default is ``{'g', 'mu_l', 'mu_u'}`` with
            %       valid values:
            %
            %           - ``'g'`` - 2 element cell array with constraint values
            %             :math:`\g(\x) - \u` and :math:`\l - \g(\x)`,
            %             respectively, where :math:`\g(\x)` is the quadratic
            %             form shown in :eq:`eq_qcn_eval_g` for the quadratic
            %             constraints specified
            %           - ``'g_u'`` or ``'f'`` - constraint values
            %             :math:`\g(\x) - \u`
            %           - ``'l_g'`` - constraint values :math:`\l - \g(\x)`
            %           - ``'mu_l'`` - shadow price on :math:`\l - \g(\x)`
            %           - ``'mu_u'`` - shadow price on :math:`\g(\x) - \u`
            %   name (char array) : name of the subset
            %   idx_list (cell array) : *(optional)* indices of the subset
            %
            % Outputs:
            %     : Variable number of outputs corresponding to ``tags`` input.
            %       If ``tags`` is empty or not specified, the calling context
            %       will define the number of outputs, returned in order of
            %       default tags.
            %
            % Example::
            %
            %     [g, mu_l, mu_u] = qcn.get_soln(var, soln, 'flow');
            %     mu_lq_Pmis_5_3 = qcn.get_soln(var, soln, 'mu_l', 'Pmis', {5,3});
            %
            % For a complete set of solution values, using the parse_soln()
            % method may be more efficient.
            %
            % See also parse_soln.

            %% input arg handling
            [tags, name, idx, N, i1, iN] = obj.get_soln_std_args(varargin{:});

            %% get outputs
            varargout = cell(1, nargout);
            if N && ~isempty(soln.eflag)
                if any(ismember({'g', 'g_u', 'l_g'}, tags(1:nargout)))
                    g = cell(1,4);
                    [g{:}] = obj.eval(var, soln.x, name, idx);
                end
                for k = 1:nargout
                    switch tags{k}
                        case 'g'
                            varargout{k} = g;
                        case 'g_u'
                            varargout{k} = g{1};
                        case 'l_g'
                            varargout{k} = g{4};
                        case 'f'
                            varargout{k} = soln.f(i1:iN);
                        case 'mu_l'
                            varargout{k} = soln.lambda.mu_lq(i1:iN);
                        case 'mu_u'
                            varargout{k} = soln.lambda.mu_uq(i1:iN);
                        otherwise
                            error('mp.sm_quad_constraint.get_soln: unknown tag ''%s''', tags{k});
                    end
                end
            end     %% if N
        end

        function ps = parse_soln(obj, soln, stash)
            % Parse solution for quadratic constraints.
            % ::
            %
            %   ps = qcn.parse_soln(soln)
            %
            % Parse a full solution struct into parts corresponding to
            % individual quadratic constraint subsets.
            %
            % Input:
            %   soln (struct) : full solution struct with these fields
            %       (among others):
            %
            %           - ``x`` - variable values
            %           - ``lambda`` - constraint shadow prices, struct with
            %             fields:
            %
            %               - ``eqnonlin`` - nonlinear equality constraints
            %               - ``ineqnonlin`` - nonlinear inequality constraints
            %               - ``mu_l`` - linear constraint lower bounds
            %               - ``mu_u`` - linear constraint upper bounds
            %               - ``mu_lq`` - quadratic constraint lower bounds
            %               - ``mu_uq`` - quadratic constraint upper bounds
            %               - ``lower`` - variable lower bounds
            %               - ``upper`` - variable upper bounds
            %   stash (boolean) : if true, store return value in :attr:`soln`
            %       property
            %
            % Output:
            %   ps (struct) : parsed solution, struct where each field listed
            %       below is a struct whos names are the names of the relevant
            %       linear constraint subsets and values are scalars for named
            %       sets, arrays for named/indexed sets:
            %
            %           - ``mu_l`` - constraint lower bound shadow prices
            %           - ``mu_u`` - constraint upper bound shadow prices

            ps = [];
            if obj.get_N()
                if isfield(soln.lambda, 'mu_lq')
                    if isfield(soln.lambda, 'mu_uq')
                        params = struct('src', {soln.lambda.mu_lq, soln.lambda.mu_uq}, ...
                                        'dst', {'mu_l', 'mu_u'});
                    else
                        params = struct('src', soln.lambda.mu_lq, 'dst', 'mu_l');
                    end
                else
                    if isfield(soln.lambda, 'mu_uq')
                        params = struct('src', soln.lambda.mu_uq, 'dst', 'mu_u');
                    else
                        params = [];
                    end
                end
                if ~isempty(params)
                    ps = obj.parse_soln_fields(params);
                end
            end

            if nargin > 2 && stash
                obj.soln = ps;
            end
        end
    end     %% methods

    methods (Static)
        function M = blkprod2vertcat(blk1, blk2, n)
            % Compute all products of two block diagonal matrices.
            % ::
            %
            %   M = qcn.blkprod2vertcat(blk1, blk2, n)
            %
            % Take two block diagonal matrices and return a matrix formed
            % by stacking vertically the result of all products of the two.
            % Here the matrices are assumed to have compatible sizes. That
            % is, the number of columns of the set of matrices used to form
            % the first block diagonal matrix must be the same as the number
            % of rows of the matrices used to form the second block diagonal
            % matrix.
            %
            % Inputs:
            %   blk1 (double) : block diagonal matrix formed from a set of
            %       matrices,  each of size :math:`m_1 \times n`
            %   blk2 (double) : block diagonal matrix formed from a set of
            %       matrices, each of size :math:`n \times m_2`
            %   n (integer) : compatible dimension of the two block diagonal
            %       matrices used to find the number of blocks to be multiplied
            %
            % Outputs:
            %   M : :math:`m_1 n \times m_2` matrix holding a vertical
            %       stack of the resulting products between the two block
            %       diagonal matrices
            %
            % Examples::
            %
            %   xx = blkdiag(repmat(x, m, 1));      % x in R^n, xx in R^(mxn)
            %   QQ = blkdiag(Q{:});                 % Q is a cell array of m matrices in R^(nxn)
            %   M = qcn.blkprod2vertcat(xx, QQ, n)  % M is a n x n matrix
            %
            % See also eval.

            [row1, col1] = size(blk1);
            [row2, col2] = size(blk2);

            if (col1) ~= (row2)
                error('sm_quad_constraint.blkprod2vertcat: number of columns of elements in blk1 (%d) do not match the number of rows of elements in blk2 (%d) \n', col1/n, row2/n);
            end

            blkprod = blk1 * blk2;    % Product of blocks
            N = col1/n;               % Number of blocks in blkprod

            % id_block = sparse(logical(ones(row1/N, col2/N)));
            % id_blkprod = mat2cell(repmat(id_block, N, 1), (row1/N)*ones(N,1));
            % id_blkprod = blkdiag(id_blkprod{:});
            %
            % M = blkprod(id_blkprod);
            % M = reshape(M, N, [])';

            %                 rows                  cols
            id_start = [(0:N-1)*(row1/N)+1 ; (0:N-1)*(col2/N)+1 ];
            id_end   = [   (1:N)*(row1/N)  ;    (1:N)*(col2/N)  ];

            id_start = mat2cell(id_start(:), 2*ones(N,1));
            id_end   = mat2cell(id_end(:), 2*ones(N,1));

            M = cellfun(@(x,y)(blkprod(x(1):y(1), x(2):y(2))), id_start, id_end, 'UniformOutput', false);

            M = cell2mat(M);
        end
    end

    methods (Access=protected)
        function default_tags = get_soln_default_tags(obj)
            % Return default tags for get_soln().
            % ::
            %
            %   default_tags = nln.get_soln_default_tags()
            %
            % Output:
            %   default_tags (cell array) : tags defining the default outputs
            %       of get_soln(), namely ``{'g', 'lam', 'dg'}``
            %
            % See also get_soln.

            default_tags = {'g', 'mu_l', 'mu_u'};
        end
    end     %% methods (Access=protected)
end         %% classdef