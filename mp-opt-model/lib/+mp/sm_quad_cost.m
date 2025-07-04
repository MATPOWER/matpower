classdef sm_quad_cost < mp.set_manager_opt_model
% mp.sm_quad_cost -  MP Set Manager class for quadratic costs.
% ::
%
%   qdc = mp.sm_quad_cost()
%   qdc = mp.sm_quad_cost(label)
%
% MP Set Manager class for quadratic costs. The costs take one of two forms,
% either a scalar cost function of the form
%
% .. math:: f(\x) = \frac{1}{2}\trans{\x} \Hh \x + \trans{\c} \x + \param{k}
%   :label: eq_qdc_form_1
%
% or a vector cost function of the form
%
% .. math:: \rvec{f}(\x) = \frac{1}{2} \diag{\param{\h}} \x^2 + \diag{\c} \x + \k
%   :label: eq_qdc_form_2
%
% where :math:`\x` is an :math:`n_x \times 1` vector, and the corresponding
% coefficient parameters are conformable.
%
% Manages cost parameters :math:`\Hh, \c, \param{k}` or :math:`\param{\h}, \c, \k`,
% along with indexing.
%
% By convention, ``qdc`` is the variable name used for mp.sm_quad_cost objects.
%
% mp.sm_quad_cost Properties:
%   * cache - struct for caching aggregated parameters for the set
%
% mp.sm_quad_cost Methods:
%   * sm_quad_cost - constructor
%   * add - add a subset of quadratic costs
%   * params - build and return cost parameters :math:`\Hh, \c, \param{k}` or :math:`\param{\h}, \c, \k`
%   * set_params - modify quadratic cost parameter data
%   * eval - evaluate individual or full set of quadratic costs
%   * display_soln - display solution values for quadratic costs
%   * get_soln - fetch solution values for specific named/indexed subsets
%
% See also mp.set_manager, mp.set_manager_opt_model.

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

    properties
        % struct for caching aggregated parameters for quadratic costs
        cache = [];
    end     %% properties

    methods
        function obj = sm_quad_cost(varargin)
            % Constructor.
            % ::
            %
            %   qdc = mp.sm_quad_cost(label)

            es = struct();  %% empty struct
            obj@mp.set_manager_opt_model(varargin{:});
            if isempty(fieldnames(obj.data))
                obj.data = struct( ...
                    'H', es, ...
                    'c', es, ...
                    'k', es, ...
                    'vs', es );
            end
        end

        function new_obj = copy(obj, varargin)
            % Duplicate the object.
            % ::
            %
            %   new_sm = sm.copy()
            %   new_sm = sm.copy(new_class)
            %
            % Input:
            %   new_class (char array) : *(default = same class)* name of class
            %       to use for new object
            %
            % Make a shallow copy of the object by copying each of the
            % top-level properties. Parent method makes the copy, then
            % the ``H`` parameter is renamed to ``Q`` if the new class is
            % mp.sm_quad_cost_legacy.

            new_obj = copy@mp.set_manager_opt_model(obj, varargin{:});
            if isa(new_obj, 'mp.sm_quad_cost_legacy')
                new_obj.data.Q = new_obj.data.H;
                new_obj.data = rmfield(new_obj.data, 'H');
            end
        end

        function obj = add(obj, var, name, idx, varargin)
            % Add a subset of quadratic costs.
            % ::
            %
            %   qdc.add(var, name, H, c);
            %   qdc.add(var, name, H, c, k);
            %   qdc.add(var, name, H, c, k, vs);
            %
            %   qdc.add(var, name, idx_list, H, c);
            %   qdc.add(var, name, idx_list, H, c, k);
            %   qdc.add(var, name, idx_list, H, c, k, vs);
            %
            % Add a named, and possibly indexed, subset of quadratic costs
            % of form :eq:`eq_qdc_form_1` or :eq:`eq_qdc_form_2` to the set
            % where :math:`\x` is an :math:`n_x \times 1` vector made up of the
            % variables specified in the optional ``vs`` *(in the order
            % given)*. This allows the :math:`\Hh`, :math:`\param{\h}`, :math:`\c`,
            % and/or :math:`\k` parameters to be defined in terms of only the
            % relevant variables without the need to manually create a lot of
            % properly located zero rows/columns.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   name (char array) : name of subset/block of costs to add
            %   idx_list (cell array) : *(optional)* index list for subset/block
            %       of costs to add (for an indexed subset)
            %   H (double) : *(optional, default = all zeros)* quadratic cost
            %       coefficient, :math:`n_x \times n_x` matrix :math:`\Hh`, or
            %       :math:`n_x \times 1` vector :math:`\param{\h}`
            %   c (double) : *(optional, default = all zeros)* linear cost
            %       coefficient, :math:`n_x \times 1` vector :math:`\c`
            %   k (double) : *(optional, default = 0)* constant cost term,
            %       scalar :math:`\param{k}`, or :math:`n_x \times 1` vector
            %       :math:`\k`; scalar `k` is expanded to a vector if `H` is
            %       a vector
            %   vs (cell or struct array) : *(optional, default* ``{}`` *)*
            %       variable set defining vector :math:`\x` for this
            %       cost subset; can be either a cell array of names of
            %       variable subsets, or a struct array of ``name``, ``idx``
            %       pairs of indexed named subsets of variables; order of
            %       ``vs`` determines order of blocks in :math:`\x`; if
            %       empty, :math:`\x` is assumed to be the full variable vector
            %
            % Examples::
            %
            %   qdc.add(var, 'quad_cost1', H1, c1, 0);
            %   qdc.add(var, 'lin_cost2',  [], c2, k2, {'Vm', 'Pg', 'z'});
            %
            %   qdc.init_indexed_name('c', {2, 3});
            %   for i = 1:2
            %     for j = 1:3
            %       qdc.add(var, 'c', {i, j}, H{i,j}, ...);
            %     end
            %   end
            %
            % See also params, set_params, eval.

            %% set up default args
            if iscell(idx)          %% indexed named set
                H = varargin{1};
                args = varargin(2:end);
            else                    %% simple named set
                H = idx;
                idx = {};
                args = varargin;
            end
            nargs = length(args);

            %% prepare data
            c = []; k = []; vs = {};
            if nargs >= 1
                c = args{1};
                if nargs >= 2
                    k = args{2};
                    if nargs >= 3
                        vs = args{3};
                    end
                end
            end
            if isempty(k)
                if isempty(H)
                    k = zeros(size(c));
                else
                    k = 0;
                end
            end

            %% convert varsets from cell to struct array if necessary
            vs = mp.sm_variable.varsets_cell2struct(vs);
            nv = var.varsets_len(vs);   %% number of variables

            %% check sizes
            [MH, NH] = size(H);
            [Mc, Nc] = size(c);
            [Mk, Nk] = size(k);
            if MH
                if NH ~= MH && NH ~= 1
                    error('mp.sm_quad_cost.add: H (%d x %d) must be square or a column vector (or empty)', MH, NH);
                end
            end
            if Mc && Nc ~= 1
                error('mp.sm_quad_cost.add: c (%d x %d) must be a column vector (or empty)', Mc, Nc);
            end
            if MH
                if Mc && Mc ~= MH
                    error('mp.sm_quad_cost.add: dimensions of H (%d x %d) and c (%d x %d) are not compatible', MH, NH, Mc, Nc);
                end
                nx = MH;
            else
                if nv && ~Mc
                    error('mp.sm_quad_cost.add: H and c cannot both be empty');
                end
                nx = Mc;
            end
            if nx ~= nv
                error('mp.sm_quad_cost.add: dimensions of H (%d x %d) and c (%d x %d) do not match\nnumber of variables (%d)\n', MH, NH, Mc, Nc, nv);
            end

            %% size of named cost set
            if isempty(H)
                if Mk == 0 || Mk == Mc
                    N = nx;     %% cost is element-wise, i.e. a vector
                elseif isscalar(k)
                    N = 1;      %% cost is scalar
                else
                    error('mp.sm_quad_cost.add: k must be a scalar or (%d x 1) vector', nx);
                end
            elseif NH == 1
                N = nx;     %% cost is element-wise, i.e. a vector
                if isscalar(k)
                    k = k * ones(N, 1); %% expand k to vector if necessary
                end
            else
                N = 1;      %% cost is scalar
                if ~isscalar(k)
                    error('mp.sm_quad_cost.add: k must be a scalar when H is a matrix');
                end
            end

            %% call parent to handle standard indexing
            if isempty(idx)
                add@mp.set_manager_opt_model(obj, name, N, args{:});
            else
                add@mp.set_manager_opt_model(obj, name, idx, N, args{:});
            end

            %% assign data
            if isempty(idx)
                obj.data.H.(name)  = H;
                obj.data.c.(name)  = c;
                obj.data.k.(name)  = k;
                obj.data.vs.(name) = vs;
            else
                %% calls to substruct() are relatively expensive, so we pre-build the
                %% struct for addressing cell array fields
                %% sc = substruct('.', name, '{}', idx);
                sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
                obj.data.H  = subsasgn(obj.data.H, sc, H);
                obj.data.c  = subsasgn(obj.data.c, sc, c);
                obj.data.k  = subsasgn(obj.data.k, sc, k);
                obj.data.vs = subsasgn(obj.data.vs, sc, vs);
            end
            if ~isempty(obj.cache)  %% clear cache of aggregated params
                obj.cache = [];
            end
        end

        function [H, c, K, vs] = params(obj, var, name, idx)
            % Build and return quadratic cost parameters :math:`\Hh, \c, \param{k}` or :math:`\param{\h}, \c, \k`.
            % ::
            %
            %   [H, c] = qdc.params(var)
            %   [H, c] = qdc.params(var, name)
            %   [H, c] = qdc.params(var, name, idx_list)
            %   [H, c, k] = qdc.params(...)
            %   [H, c, k, vs] = qdc.params(...)
            %
            % With no input parameters, it assembles and returns the parameters
            % for the aggregate quadratic cost from all quadratic cost sets
            % added using add(). The values of these parameters are cached
            % for subsequent calls. The parameters are :math:`\Hh, \c`, and
            % optionally :math:`\param{k}` for the scalar cost form in
            % :eq:`eq_qdc_form_1`, or :math:`\param{\h}, \c`, and optionally
            % :math:`\k` for the vector cost form in :eq:`eq_qdc_form_2`.
            %
            % If a name or name and index list are provided, then it simply
            % returns the parameters for the corresponding set. It can also
            % optionally return the variable sets used by this cost set
            % (the dimensions of :math:`\Hh, \param{\h}, \c, \k` will be consistent
            % with this variable set).
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   name (char array) : *(optional)* name of subset
            %   idx_list (cell array) : *(optional)* index list for subset
            %
            % Outputs:
            %   H (double) : quadratic cost coefficient matrix :math:`\Hh` or
            %       vector :math:`\param{\h}`
            %   c (double) : linear cost coefficient vector :math:`\c`
            %   k (double) : constant cost term, scalar :math:`\param{k}`, or vector
            %       :math:`\k`
            %   vs (struct array) : variable set, ``name``, ``idx`` pairs
            %       specifying the set of variables defining vector :math:`\x`
            %       for this constraint subset; order of ``vs`` determines
            %       order of blocks in :math:`\x`
            %
            % See also add.

            if nargin > 2       %% individual set
                if nargin < 4
                    idx = {};
                end
                if isempty(idx)                 %% name, no index provided
                    if numel(obj.idx.i1.(name)) == 1    %% simple named set
                        H = obj.data.H.(name);
                        c = obj.data.c.(name);
                        K = obj.data.k.(name);
                        if nargout > 3
                            vs = obj.data.vs.(name);
                        end
                    else                                    %% indexing required
                        error('mp.sm_quad_cost.params: quadratic cost set ''%s'' requires an IDX_LIST arg', name);
                    end
                else                            %% indexed named set
                    % (calls to substruct() are relatively expensive ...
                    % s = substruct('.', name, '{}', idx);
                    % ... so replace it with these more efficient lines)
                    sc = struct('type', {'.', '{}'}, 'subs', {name, idx});
                    H = subsref(obj.data.H, sc);
                    c = subsref(obj.data.c, sc);
                    K = subsref(obj.data.k, sc);
                    if nargout > 3
                        vs = subsref(obj.data.vs, sc);
                    end
                end
            else                %% aggregate
                cache = obj.cache;
                if isempty(cache)       %% build the aggregate
                    nx = var.N;         %% number of variables
                    if obj.NS < 25 || obj.NS < 100 && nx < 300
                        %% METHOD 1: Add sparse matrices (original method)
                        Ht = sparse(nx, nx);    %% transpose of quadratic coefficients
                        c = zeros(nx, 1);       %% linear coefficients
                        K = 0;                  %% constant term
                        for k = 1:obj.NS
                            name = obj.order(k).name;
                            idx  = obj.order(k).idx;
                            N = obj.get_N(name, idx);
                            [Hk, ck, kk, vs] = obj.params(var, name, idx);
                            haveH = ~isempty(Hk);
                            havec = ~isempty(ck);
                            nk = max(size(Hk, 1), size(ck, 1));     %% size of Hk and/or ck
                            if isempty(vs)
                                if nk == nx     %% full size
                                    if size(Hk, 2) == 1     %% Hk is a column vector
                                        Hkt_full = spdiags(Hk, 0, nx, nx);
                                    elseif haveH            %% Hk is a matrix
                                        Hkt_full = Hk';
                                    end
                                    if havec
                                        ck_full = ck;
                                    end
                                else            %% vars added since adding this cost set
                                    if size(Hk, 2) == 1     %% Hk is a column vector
                                        Hkt_full = sparse(1:nk, 1:nk, Hk, nx, nx);
                                    elseif haveH            %% Hk is a matrix
                                        Hk_all_cols = sparse(nk, nx);
                                        Hk_all_cols(:, 1:nk) = Hk;
                                        Hkt_full(:, 1:nk) = Hk_all_cols';
                                    end
                                    if havec
                                        ck_full = zeros(nx, 1);
                                        ck_full(1:nk) = ck;
                                    end
                                end
                            else
                                jj = var.varsets_idx(vs);   %% indices for var set
                                if size(Hk, 2) == 1     %% Hk is a column vector
                                    Hkt_full = sparse(jj, jj, Hk, nx, nx);
                                elseif haveH            %% Hk is a matrix
                                    Hk_all_cols = sparse(nk, nx);
                                    Hk_all_cols(:, jj) = Hk;
                                    Hkt_full = sparse(nx, nx);
                                    Hkt_full(:, jj) = Hk_all_cols';
                                end
                                if havec
                                    ck_full = zeros(nx, 1);
                                    ck_full(jj) = ck;
                                end
                            end
                            if haveH
                                Ht = Ht + Hkt_full;
                            end
                            if havec
                                c = c + ck_full;
                            end
                            K = K + sum(kk);
                        end
                        H = Ht';
                   else
                        %% METHOD 2: construct using single call to sparse()
                        H_ijv = cell(obj.NS, 3);    %% indices/values to construct H
                        c_ijv = cell(obj.NS, 3);    %% indices/values to construct c
                        K = 0;                  %% constant term
                        for k = 1:obj.NS
                            name = obj.order(k).name;
                            idx  = obj.order(k).idx;
                            N = obj.get_N(name, idx);
                            [Hk, ck, kk, vs] = obj.params(var, name, idx);
                            haveH = ~isempty(Hk);
                            havec = ~isempty(ck);
                            if haveH
                                [i, j, v] = find(Hk);
                            end
                            if havec
                                [ic, jc, vc] = find(ck);
                            end
                            if isempty(vs)
                                if size(Hk, 2) == 1     %% Hk is a column vector
                                    H_ijv(k, :) = {i, i, v};
                                elseif haveH            %% Hk is a matrix
                                    H_ijv(k, :) = {i, j, v};
                                end
                                if havec
                                    c_ijv(k, :) = {ic, jc, vc};
                                end
                            else
                                jj = var.varsets_idx(vs)';  %% indices for var set
                                if size(Hk, 2) == 1     %% Hk is a column vector
                                    H_ijv(k, :) = {jj(i), jj(i), v};
                                elseif haveH            %% Hk is a matrix
                                    H_ijv(k, :) = {jj(i), jj(j), v};
                                end
                                if havec
                                    c_ijv(k, :) = {jj(ic), jc, vc};
                                end
                            end
                            K = K + sum(kk);
                        end
                        H = sparse( vertcat(H_ijv{:,1}), ...
                                    vertcat(H_ijv{:,2}), ...
                                    vertcat(H_ijv{:,3}), nx, nx);
                        c = accumarray(vertcat(c_ijv{:,1}), vertcat(c_ijv{:,3}), [nx 1]);
                    end

                    %% cache aggregated parameters
                    obj.cache = struct('H', H, 'c', c, 'k', K);
                else                    %% return cached values
                    H = cache.H;
                    c = cache.c;
                    K = cache.k;
                end
                if nargout > 3
                    vs = {};
                end
            end
        end

        function obj = set_params(obj, var, name, idx, params, vals)
            % Modify quadratic cost parameter data.
            % ::
            %
            %   qdc.set_params(var, name, params, vals)
            %   qdc.set_params(var, name, idx_list, params, vals)
            %
            % This method can be used to modify parameters for an existing
            % subset of quadratic costs.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   name (char array) : name of subset/block of quadratic costs to
            %       modify
            %   idx_list (cell array) : *(optional)* index list for subset/block
            %       of quadratic costs to modify (for an indexed subset)
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
            % Valid parameter names are ``H``, ``c``, ``k``, ``vs``.
            %
            % Examples::
            %
            %   qdc.set_params(var, 'y', {2,3}, {'c'}, {c});
            %   qdc.set_params(var, 'Pg', 'all', {H, c, k, vs});
            %
            % See also add, params.

            if nargin < 6
                vals = params;
                params = idx;
                idx = {};
            end

            %% create default list of parameters to update based on set type & inputs
            default_params = {'H', 'c', 'k', 'vs'};

            %% standardize provided arguments in cell arrays params, vals
            [is_all, np, params, vals] = ...
                obj.set_params_std_args(default_params, params, vals);

            %% get current parameters
            [H, c, kk, vs] = obj.params(var, name, idx);
            [MH0, NH0] = size(H);
            [Mc0, Nc0] = size(c);
            Nk0 = length(kk);
            nx0 = max([MH0 Mc0 Nk0]);
            N0 = obj.get_N(name, idx);
            if isempty(vs), vs = {vs}; end
            p = struct('H', H, 'c', c, 'k', kk, 'vs', vs);  %% current parameters
            u = struct('H', 0, 'c', 0, 'k',  0, 'vs',  0);  %% which ones to update

            %% replace with new parameters
            for k = 1:np
                p.(params{k}) = vals{k};
                u.(params{k}) = 1;
            end

            %% set missing default params for 'all'
            [MH, NH] = size(p.H);
            [Mc, Nc] = size(p.c);
            Nk = length(p.k);
            nx = max([MH Mc Nk]);
            if NH == 1 || (isempty(p.H) && (Nk > 1 || k == 0))
                %% H is a column vector (cost is element-wise, i.e. a vector)
                %% OR H is empty and k is either a vector or zero
                N = nx;
            else            %% H is a square matrix (cost is a scalar)
                N = 1;
            end
            if is_all
                u.H = 1;                %% always update H
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
                error('mp.sm_quad_cost.set_params: dimension change for ''%s'' not allowed except for ''all''', obj.nameidxstr(name, idx));
            end

            %% H and c can't both be empty
            if ~MH && ~Mc
                error('mp.sm_quad_cost.set_params: ''%s'' : ''H'' and ''c'' cannot both be empty', obj.nameidxstr(name, idx));
            end

            %% check sizes of new values of H, c, k
            if ~isempty(p.H) && (MH ~= nx || (MH ~= NH && NH ~= 1) )
                error('mp.sm_quad_cost.set_params: ''%s'' : ''%s'' is expected to be (%d x %d)', obj.nameidxstr(name, idx), 'H', MH, NH);
            end
            if ~isempty(p.c) && Mc ~= nx
                error('mp.sm_quad_cost.set_params: ''%s'' : ''%s'' is expected to be (%d x %d)', obj.nameidxstr(name, idx), 'c', Mc, 1);
            end
            if ~isempty(p.k) && any(p.k) && Nk ~= N && Nk ~= 1
                error('mp.sm_quad_cost.set_params: ''%s'' : ''%s'' is expected to be (%d x %d)', obj.nameidxstr(name, idx), 'k', N, 1);
            end

            %% check consistency of H, c, k and vs
            if u.H || u.c || u.vs
                p.vs = mp.sm_variable.varsets_cell2struct(p.vs);
                nv = var.varsets_len(p.vs);     %% number of variables
                if nx ~= nv
                    error('mp.sm_quad_cost.set_params: for ''%s'' dimensions of ''H'', ''c'', ''k'' (%d) must be consistent with ''vs'' (%d)', obj.nameidxstr(name, idx), nx, nv);
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
            dN = N - N0;
            if is_all && dN
                obj.set_params_update_dims(dN, name, idx);
            end
        end

        function [f, df, d2f] = eval(obj, var, x, name, idx)
            % Evaluate individual or full set of quadratic costs.
            % ::
            %
            %   f = qdc.eval(var, x ...)
            %   [f, df] = qdc.eval(var, x ...)
            %   [f, df, d2f] = qdc.eval(var, x ...)
            %   [f, df, d2f] = qdc.eval(var, x, name)
            %   [f, df, d2f] = qdc.eval(var, x, name, idx_list)
            %
            % For a given value of the variable vector :math:`\x`, this method
            % evaluates the quadratic cost function and optionally its
            % derivatives for an individual subset, if name or name and
            % index list are provided, otherwise, for the full set of costs.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   x (double) : full :math:`n_x \times 1` variable vector :math:`\x`
            %   name (char array) : name of subset/block of quadratic costs
            %       to evaluate
            %   idx_list (cell array) : *(optional)* index list for subset/block
            %       of quadratic costs to evaluate (for an indexed subset)
            %
            % Outputs:
            %   f (double) : scalar cost :math:`f(\x)` for :eq:`eq_qdc_form_1`,
            %       or vector cost :math:`\rvec{f}(\x)` for :eq:`eq_qdc_form_2`
            %   df (double) : *(optional)* cost first derivatives, gradient of
            %       scalar cost :math:`\der{f}{\x}` for :eq:`eq_qdc_form_1`,
            %       or vector of cost first derivatives for :eq:`eq_qdc_form_2`
            %   d2f (double) : *(optional)* second derivative of costs,
            %       :math:`\der{^2 f}{\x^2}` for :eq:`eq_qdc_form_1`, or vector
            %       of cost second derivatives for :eq:`eq_qdc_form_2`
            %
            % See also add, params.

            if obj.N
                done = 0;

                %% collect cost parameters
                if nargin < 4                       %% full set
                    [H, c, k, vs] = obj.params(var);
                    N = 1;
                elseif nargin < 5 || isempty(idx)   %% name, no idx provided
                    dims = size(obj.idx.i1.(name));
                    if prod(dims) == 1              %% simple named set
                        [H, c, k, vs] = obj.params(var, name);
                        N = obj.get_N(name);
                    elseif nargout == 1             %% indexing required, recurse
                        f = 0;          %% initialize cumulative cost
                        idx = num2cell(ones(size(dims))); %% initialize idx
                        while ~done     %% call eval() recursively
                            f = f + sum(obj.eval(var, x, name, idx));

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
                        error('mp.sm_quad_cost.eval: quadratic cost set ''%s'' requires an IDX_LIST arg when requesting DF output', name)
                    end
                else                                %% indexed named set
                    [H, c, k, vs] = obj.params(var, name, idx);
                    N = obj.get_N(name, idx);
                end

                if ~done
                    %% assemble appropriately-sized x vector
                    xx = var.varsets_x(x, vs, 'vector');

                    %% compute/assemble f
                    if N == 1               %% f is scalar (H is matrix, k is scalar)
                        f = k;                  %% start with k term
                        if ~isempty(c)
                            f = f + c'*xx;      %% add c term
                        end
                        if ~isempty(H)          %% add H term
                            f = f + (xx'*H*xx)/2;
                        end
                    else                    %% f is vector (H is vector or empty, k is vector)
                        if isempty(c)           %% H, k terms only
                            f = (H .* xx.^2)/2 + k;
                        else
                            if isempty(H)       %% c, k terms only
                                f = c .* xx + k;
                            else                %% H, c, k terms
                                f = (H .* xx.^2)/2 + c .* xx + k;
                            end
                        end
                    end

                    if nargout > 1
                        %% compute/assemble df
                        if ~isempty(c)
                            df = c;             %% start with c term
                        else
                            df = 0;             %% start with nothing
                        end
                        if ~isempty(H)
                            if N == 1       %% f is scalar (H is matrix, k is scalar)
                                df = df + H*xx;     %% add H term
                            else            %% f is vector (H is vector or empty, k is vector or scalar)
                                df = df + H.*xx;    %% add H term
                            end
                        end

                        %% assemble d2f
                        if nargout > 2
                            if isempty(H)
                                nx = length(xx);
                                if N == 1   %% f is scalar (H is matrix, k is scalar)
                                    d2f = sparse(nx, nx);
                                else        %% f is vector (H is vector or empty, k is vector or scalar)
                                    d2f = sparse(nx, 1);
                                end
                            else
                                d2f = H;
                            end
                        end
                    end     %% nargout > 1
                end         %% ~done
            else
                f = 0;
                if nargout > 1
            %         nx = length(x);
            %         df = zeros(nx, 1);
                    df = [];
                    if nargout > 2
            %             d2f = sparse(nx, nx);
                        d2f = [];
                    end
                end
            end
        end

        function obj = display_soln(obj, var, soln, varargin)
            % Display solution values for quadratic costs.
            % ::
            %
            %   qdc.display_soln(var, soln)
            %   qdc.display_soln(var, soln, name)
            %   qdc.display_soln(var, soln, name, idx_list)
            %   qdc.display_soln(var, soln, fid)
            %   qdc.display_soln(var, soln, fid, name)
            %   qdc.display_soln(var, soln, fid, name, idx_list)
            %
            % Displays the solution values for all quadratic costs (default)
            % or an individual named or named/indexed subset.
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
            %               - ``lower`` - variable lower bounds
            %               - ``upper`` - variable upper bounds
            %   fid (fileID) : fileID of open file to write to (default is
            %       1 for standard output)
            %   name (char array) : *(optional)* name of individual subset
            %   idx_list (cell array) : *(optional)* indices of individual
            %       subset

            [fid, name, idx, idxs, hdr1] = obj.display_soln_std_args(varargin{:});

            if obj.N
                c = [];
                c_k = [];
                c_lin = [];
                c_avg = [];
                vv = var.idx;
                for k = 1:length(obj.order)
                    n = obj.order(k).name;
                    i = obj.order(k).idx;
                    [HH, cc, kk, vs] = obj.params(var, n, i);
                    xx = var.varsets_x(soln.x, vs, 'vector');
                    c_total = obj.eval(var, soln.x, n, i);
                    len = length(c_total);
                    if len == 1
                        c_constant = kk;
                        if isempty(cc)
                            c_linear = 0;
                        else
                            c_linear = cc' * xx;
                        end
                        if abs(sum(xx)) > 1e-9
                            c_average = c_total / sum(xx);
                        else
                            c_average = NaN;
                        end
                    else
                        if isempty(cc)
                            c_linear = 0;
                        else
                            c_linear = cc .* xx;
                        end
                        c_average = c_total ./ xx;
                        c_average(isinf(c_average)) = NaN;
                        if isscalar(kk)
                            c_constant = ones(len, 1)*kk/len;
                        else
                            c_constant = kk;
                        end
                    end
                    if sum(sum(HH)) == 0 && sum(kk) == 0 && len == length(cc)
                        c_average = cc;
                    end
                    c = [c; c_total];
                    c_lin = [c_lin; c_linear];
                    c_k = [c_k; c_constant];
                    c_avg = [c_avg; c_average];
                end
                c_quad = c - c_lin - c_k;

                %% print header rows
                hdr2 = {'   cost  =  quad    linear  constant  average', ...
                        ' -------- -------- -------- -------- --------' };
                obj.display_soln_print_headers(fid, hdr1, hdr2);

                %% print data
                for k = 1:length(idxs)
                    obj.display_soln_print_row(fid, idxs(k));

                    if isnan(c_avg(idxs(k)))
                        cc_avg = '- ';
                    else
                        cc_avg = obj.sprintf_num(8, c_avg(idxs(k)));
                    end
                    fprintf(fid, '%9s%9s%9s%9s%9s\n', ...
                        obj.sprintf_num(8, c(idxs(k))), ...
                        obj.sprintf_num(8, c_quad(idxs(k))), ...
                        obj.sprintf_num(8, c_lin(idxs(k))), ...
                        obj.sprintf_num(8, c_k(idxs(k))), ...
                        cc_avg);
                end

                %% print footer rows
                fprintf(fid, '%s\n', [hdr1{2} hdr2{2}]);
                fprintf(fid, '%7s %-28s%9s%9s%9s%9s%9s\n', '', ...
                    'Sum of Displayed Costs', ...
                    obj.sprintf_num(8, sum(c(idxs))), ...
                    obj.sprintf_num(8, sum(c_quad(idxs))), ...
                    obj.sprintf_num(8, sum(c_lin(idxs))), ...
                    obj.sprintf_num(8, sum(c_k(idxs))), '');
                fprintf(fid, '\n');
            end
        end

        function varargout = get_soln(obj, var, soln, varargin)
            % Fetch solution values for specific named/indexed subsets.
            % ::
            %
            %   vals = qdc.get_soln(var, soln, name)
            %   vals = qdc.get_soln(var, soln, name, idx_list)
            %   vals = qdc.get_soln(var, soln, tags, name)
            %   vals = qdc.get_soln(var, soln, tags, name, idx_list)
            %
            % Returns named/indexed quadratic cost results for a solved
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
            %               - ``lower`` - variable lower bounds
            %               - ``upper`` - variable upper bounds
            %   tags (char array or cell array of char arrays) : names of
            %       desired outputs, default is ``{'f', 'df', 'd2f'}`` with
            %       valid values:
            %
            %           - ``'f'`` - scalar cost :math:`f(\x)` for
            %             :eq:`eq_qdc_form_1`, or vector cost
            %             :math:`\rvec{f}(\x)` for :eq:`eq_qdc_form_2`
            %           - ``'df'`` - cost first derivatives, gradient of
            %             scalar cost :math:`\der{f}{\x}` for
            %             :eq:`eq_qdc_form_1`, or vector of cost first
            %             derivatives for :eq:`eq_qdc_form_2`
            %           - ``'d2f'`` - second derivative of costs,
            %             :math:`\der{^2 f}{\x^2}` for :eq:`eq_qdc_form_1`, or
            %             vector of cost second derivatives for
            %             :eq:`eq_qdc_form_2`
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
            %     [f, df, d2f] = qdc.get_soln(var, soln, 'gen');
            %     df_Pg_2_4 = qdc.get_soln(var, soln, 'df', 'Pg', {2,4});

            % Add the below once parse_soln() is implemented:
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
                if ismember('d2f', tags(1:nargout))
                    [f, df, d2f] = obj.eval(var, soln.x, name, idx);
                elseif ismember('df', tags(1:nargout))
                    [f, df] = obj.eval(var, soln.x, name, idx);
                else
                    f = obj.eval(var, soln.x, name, idx);
                end
                for k = 1:nargout
                    switch tags{k}
                        case 'f'
                            varargout{k} = f;
                        case 'df'
                            varargout{k} = df;
                        case 'd2f'
                            varargout{k} = d2f;
                        otherwise
                            error('mp.sm_quad_cost.get_soln: unknown tag ''%s''', tags{k});
                    end
                end
            end     %% if N
        end
    end     %% methods

    methods (Access=protected)
        function default_tags = get_soln_default_tags(obj)
            % Return default tags for get_soln().
            % ::
            %
            %   default_tags = qdc.get_soln_default_tags()
            %
            % Output:
            %   default_tags (cell array) : tags defining the default outputs
            %       of get_soln(), namely ``{'f', 'df', 'd2f'}``
            %
            % See also get_soln.

            default_tags = {'f', 'df', 'd2f'};
        end
    end     %% methods (Access=protected)
end         %% classdef
