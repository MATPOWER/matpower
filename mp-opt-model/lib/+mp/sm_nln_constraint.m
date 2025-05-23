classdef sm_nln_constraint < mp.set_manager_opt_model
% mp.sm_nln_constraint -  MP Set Manager class for nonlinear constraints.
% ::
%
%   nln = mp.sm_nln_constraint()
%   nln = mp.sm_nln_constraint(label)
%
% MP Set Manager class for nonlinear constraints
%
% .. math:: \g(\x) = 0
%
% or
%
% .. math:: \g(\x) \le 0
%
% Manages nonlinear constraint sets and their indexing.
%
% By convention, ``nln`` (general), ``nle`` (equality), or ``nli`` (inequality)
% are the variable names used for mp.sm_nln_constraint objects.
%
% mp.sm_nln_constraint Methods:
%   * sm_nln_constraint - constructor
%   * add - add a subset of nonlinear constraints
%   * params - return nonlinear constraint parameters
%   * set_params - modify nonlinear constraint parameter data
%   * eval - evaluate individual or full set of nonlinear constraints
%   * eval_hess - evaluate "Hessian" for full set of nonlinear constraints
%   * display_soln - display solution values for nonlinear constraints
%   * get_soln - fetch solution values for specific named/indexed subsets
%   * parse_soln - parse solution for nonlinear constraints
%
% See also mp.set_manager, mp.set_manager_opt_model.

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

    methods
        function obj = sm_nln_constraint(varargin)
            % Constructor.
            % ::
            %
            %   nln = mp.sm_nln_constraint(label)

            obj@mp.set_manager_opt_model(varargin{:});
            obj.data = struct( ...
                'fcn', [], ...
                'hess', [], ...
                'include', [], ...
                'vs', struct() );
        end

        function obj = add(obj, var, name, idx, varargin)
            % Add a subset of nonlinear constraints.
            % ::
            %
            %   nln.add(var, name, N, fcn, hess);
            %   nln.add(var, name, N, fcn, hess, vs);
            %
            %   nln.add(var, name, idx_list, N, fcn, hess);
            %   nln.add(var, name, idx_list, N, fcn, hess, vs);
            %
            % Add a named, and possibly indexed, subset of nonlinear
            % constraints  :math:`\g(\x) = 0` or :math:`\g(\x) \le 0` to the
            % set, where :math:`\x` is an :math:`n_x \times 1` vector made up
            % of the variables specified in the optional ``vs`` *(in the order
            % given)*. This allows the constraint function to be defined in
            % terms of only the relevant variables without the need to manually
            % account for the locations of other variable sets.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   name (char array or cell array) : name(s) of subset/block of
            %       constraints to add
            %   idx_list (cell array) : *(optional)* index list for subset/block
            %       of constraints to add (for an indexed subset)
            %   N (integer) : number of constraints in the set, i.e. the
            %       dimension :math:`n` of constraint function :math:`\g(\x)`
            %   fcn (function handle) : handle to function that evaluates the
            %       :math:`n \times 1` constraint vector :math:`\g(\x)`, and
            %       optionally the corresponding :math:`n \times n_x` Jacobian
            %       :math:`\g_\x = \der{\g}{\x}`
            %   hess (function handle) : handle to function that evaluates the
            %       constraint "Hessian" :math:`\g_{\x\x}(\lam)` as
            %       described below
            %   vs (cell or struct array) : *(optional, default* ``{}`` *)*
            %       variable set defining vector :math:`\x` for this
            %       constraint subset; can be either a cell array of names of
            %       variable subsets, or a struct array of ``name``, ``idx``
            %       pairs of indexed named subsets of variables; order of
            %       ``vs`` determines order of blocks in :math:`\x`; if
            %       empty, :math:`\x` is assumed to be the full variable vector
            %
            % **Constraint Function Implmentation** : ``fcn``
            %
            % For a constraint function :math:`\g(\x)`, the ``fcn`` input
            % should point to a function with the following interface::
            %
            %   g = fcn(x)
            %   [g, dg] = fcn(x)
            %
            % Input:
            %   x (double or cell array of doubles) : variable :math:`\x`
            %       corresponding to the full variable vector, or to a set of
            %       sub-vectors defined by the variable set provided in ``vs``
            %
            % Outputs:
            %   g (double) : :math:`n \times 1` constraint vector :math:`\g(\x)`
            %   dg (double) : :math:`n \times n_x` constraint Jacobian
            %       :math:`\g_\x = \der{\g}{\x}`
            %
            %       .. note:: The ``dg`` return value is the transpose of what
            %           is expected from an input function for nlps_master and
            %           friends.
            %
            % **Hessian Function Implmentation** : ``hess``
            %
            % Similarly, the ``hess`` input should point to a function that
            % returns an :math:`n_x \times n_x` matrix constraint "Hessian"
            % matrix :math:`\g_{\x\x}(\lam)` for a given set of multipliers,
            % with the following interface::
            %
            %   d2g = hess(x, lam)
            %
            % Inputs:
            %   x (double or cell array of doubles) : variable :math:`\x`
            %       corresponding to the full variable vector, or to a set of
            %       sub-vectors defined by the variable set provided in ``vs``
            %   lam (double) : :math:`n \times 1` vector of multipliers
            %       :math:`\lam`
            %
            % Output:
            %   d2g (double) : :math:`n_x \times n_x` "Hessian" matrix
            %       :math:`\g_{\x\x}(\lam) = \der{}{\x}(\trans{\g_\x} \lam)`
            %
            % For both functions, the input argument ``x`` takes one of two
            % forms. If the constraint set is added with varset input ``vs``
            % empty or missing, then ``x`` will be the full variable vector.
            % Otherwise it will be a cell array of vectors corresponding to
            % the variable sets specified in ``vs``.
            %
            % **Special Case**
            %
            % For simple (not indexed) named sets, ``name`` can be a cell
            % array of constraint set names, in which case ``N`` is a vector
            % specifying the number of constraints in each corresponding set.
            % ``fcn`` and ``hess`` are each still a single function handle, but
            % the values computed by each correspond to the entire stacked
            % collection of constraint sets together, as if they were a single
            % set.
            %
            % Likewise, if ``fcn`` or ``hess`` are empty, it indicates a
            % placeholder in the indexing for a constraint set whose
            % implementation is included in another constraint set. This
            % functionality is only intended to be used internally to handle
            % constraint/gradient and Hessian functions that compute the
            % values for more than one constraint set simultaneously.
            %
            % Examples::
            %
            %   % nonlinear equality constraint with constraint/gradient and
            %   % Hessian evaluation functions provided
            %   nle.add('var, Qmis', nb, 1, fcn, hess);
            %
            %   % nonlinear inequality constraints with indexed named set 'S(i,j)'
            %   nli.init_indexed_name('S', {2, 3});
            %   for i = 1:2
            %     for j = 1:3
            %       nli.add('S', {i, j}, N{i,j}, 0, ...);
            %     end
            %   end
            %
            % See also params, set_params, eval.

            %% set up default args
            if iscell(idx)          %% indexed named set
                N = varargin{1};
                args = varargin(2:end);
            else                    %% simple named set
                N = idx;
                idx = {};
                args = varargin;
            end
            nargs = length(args);

            %% prepare data
            [fcn, hess] = deal(args{1:2});
            if nargs >= 3
                vs = args{3};
            else
                vs = {};
            end

            %% call parent to handle standard indexing
            if iscell(name)
                if length(name) ~= length(N)
                    error('mp.sm_nln_constraint.add: dimensions of NAME and N must match');
                end
                for k = 1:length(name)
                    add@mp.set_manager_opt_model(obj, name{k}, N(k), args{:});
                end
            else
                add@mp.set_manager_opt_model(obj, name, idx, N, args{:});
            end

            %% convert varsets from cell to struct array if necessary
            vs = mp.sm_variable.varsets_cell2struct(vs);

            %% assign data
            if iscell(name)
                computed_by = name{1};
                obj.data.fcn.(computed_by)  = fcn;
                obj.data.hess.(computed_by) = hess;
                obj.data.vs.(computed_by) = vs;
                obj.data.include.(computed_by).name = {};
                obj.data.include.(computed_by).N = [];
                for k = 2:length(name)
                    obj.data.include.(computed_by).name{end+1} = name{k};
                    obj.data.include.(computed_by).N(end+1) = N(k);
                end
            elseif isempty(idx)
                obj.data.fcn.(name)  = fcn;
                obj.data.hess.(name) = hess;
                obj.data.vs.(name) = vs;
            else
                %% calls to substruct() are relatively expensive, so we pre-build the
                %% struct for addressing cell array fields
                %% sc = substruct('.', name, '{}', idx);
                sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
                obj.data.fcn  = subsasgn(obj.data.fcn, sc, fcn);
                obj.data.hess = subsasgn(obj.data.hess, sc, hess);
                obj.data.vs   = subsasgn(obj.data.vs, sc, vs);
            end
        end

        function [N, fcn, hess, vs, include] = params(obj, var, name, idx)
            % Return nonlinear constraint parameters.
            % ::
            %
            %   N = nln.params(var, name)
            %   N = nln.params(var, name, idx_list)
            %   [N, fcn] = nln.params(...)
            %   [N, fcn, hess] = nln.params(...)
            %   [N, fcn, hess, vs] = nln.params(...)
            %   [N, fcn, hess, vs, include] = nln.params(...)
            %
            % Returns the parameters for the general nonlinear constraint subset
            % corresponding to the name or name and index list provided. The
            % parameters are those supplied via the add method when the
            % subset was added.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   name (char array) : name of subset
            %   idx_list (cell array) : *(optional)* index list for subset
            %
            % Outputs:
            %   N (integer) : number of constraints in the set, i.e. the
            %       dimension :math:`n` of constraint function :math:`\g(\x)`
            %   fcn (function handle) : handle to function that evaluates the
            %       :math:`n \times 1` constraint vector :math:`\g(\x)`, and
            %       optionally the corresponding :math:`n \times n_x` Jacobian
            %       :math:`\g_\x = \der{\g}{\x}`
            %   hess (function handle) : handle to function that evaluates the
            %       constraint "Hessian" :math:`\g_{\x\x}(\lam)`; see
            %       add() for details
            %   vs (struct array) : variable set, ``name``, ``idx`` pairs
            %       specifying the set of variables defining vector :math:`\x`
            %       for this cost subset; order of ``vs`` determines
            %       order of blocks in :math:`\x`
            %   include (struct) : *(optional)* for constraint sets whose
            %       functions compute the constraints for additional sets,
            %       struct with fields:
            %
            %       - ``name`` - cell array of additional set names computed
            %       - ``N`` - array of corresponding dimensions
            %
            % See also add, eval.

            if nargin < 4
                idx = {};
            end

            if isempty(idx)
                if ~isscalar(obj.idx.N.(name))
                    error('mp.sm_nln_constraint.params: nonlinear constraint set ''%s'' requires an IDX_LIST arg', name);
                end
                N = obj.idx.N.(name);
                if nargout > 1
                    if isfield(obj.data.fcn, name)
                        fcn = obj.data.fcn.(name);
                    else
                        fcn = '';
                    end
                    if nargout > 2
                        if isfield(obj.data.hess, name)
                            hess = obj.data.hess.(name);
                        else
                            hess = '';
                        end
                        if nargout > 3
                            if isfield(obj.data.vs, name)
                                vs = obj.data.vs.(name);
                            else
                                vs = {};
                            end
                            if nargout > 4
                                if isfield(obj.data.include, name)
                                    include = obj.data.include.(name);
                                else
                                    include = '';
                                end
                            end
                        end
                    end
                end
            else
                %% calls to substruct() are relatively expensive, so we pre-build the
                %% structs for addressing cell and numeric array fields
                %% sn = substruct('.', name, '()', idx);
                %% sc = substruct('.', name, '{}', idx);
                sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
                sn = sc; sn(2).type = '()';                             %% num array field

                N = subsref(obj.idx.N, sn);
                if nargout > 1
                    if isfield(obj.data.fcn, name)
                        fcn = subsref(obj.data.fcn, sc);
                    else
                        fcn = '';
                    end
                    if nargout > 2
                        if isfield(obj.data.hess, name)
                            hess = subsref(obj.data.hess, sc);
                        else
                            hess = '';
                        end
                        if nargout > 3
                            if isfield(obj.data.vs, name)
                                vs = subsref(obj.data.vs, sc);
                            else
                                vs = {};
                            end
                            if nargout > 4
                                error('mp.sm_nln_constraint.params: nonlinear constraint set ''%s'' cannot return INCLUDE, since a nonlinear constraint set computed by another set is currently only implemented for simple named sets, not yet for indexed named sets', name)
                            end
                        end
                    end
                end
            end
        end

        function obj = set_params(obj, var, name, idx, params, vals)
            % Modify nonlinear constraint parameter data.
            % ::
            %
            %   nln.set_params(var, name, params, vals)
            %   nln.set_params(var, name, idx_list, params, vals)
            %
            % This method can be used to modify parameters for an existing
            % subset of nonlinear constraints.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   name (char array) : name of subset/block of nonlinear
            %       constraints to modify
            %   idx_list (cell array) : *(optional)* index list for subset/block
            %       of nonlinear constraints to modify (for an indexed subset)
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
            % Valid parameter names are ``N``, ``fcn``, ``hess``, ``vs``.
            %
            % Examples::
            %
            %   nln.set_params(var, 'y', {2,3}, {'fcn', 'hess'}, {fcn, hess});
            %   nln.set_params(var, 'Pmis', 'all', {N, fcn, hess, vs});
            %
            % See also add, params.

            if nargin < 6
                vals = params;
                params = idx;
                idx = {};
            end

            %% create default list of parameters to update based on set type & inputs
            default_params = {'N', 'fcn', 'hess', 'vs'};

            %% standardize provided arguments in cell arrays params, vals
            [is_all, np, params, vals] = ...
                obj.set_params_std_args(default_params, params, vals);

            %% get current parameters
            if isempty(idx)
                [N0, fcn, hess, vs, include] = obj.params(var, name, idx);
            else
                [N0, fcn, hess, vs] = obj.params(var, name, idx);
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
                error('mp.sm_nln_constraint.set_params: dimension change for ''%s'' not allowed except for ''all''', obj.nameidxstr(name, idx));
            end

            %% included constraints not yet implemented
            if ~isempty(include)
                error('mp.sm_nln_constraint.set_params: modifications for ''%s'' not (yet) supported since it includes evaluation of other constraints', obj.nameidxstr(name, idx));
            end

            %% convert vs to struct
            if u.vs
                p.vs = mp.sm_variable.varsets_cell2struct(p.vs);
            end

            %% assign new parameters
            if isempty(idx)     %% simple named set
                for k = 2:length(default_params)
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
                for k = 2:length(default_params)
                    pn = default_params{k};     %% param name
                    if u.(pn)   %% assign new val for this parameter
                        obj.data.(pn) = subsasgn(obj.data.(pn), sc, p.(pn));
                    end
                end
            end

            %% update dimensions and indexing, if necessary
            dN = N - N0;
            if is_all && dN
                obj.set_params_update_dims(dN, name, idx);
            end
        end

        function [g, dg] = eval(obj, var, x, name, idx)
            % Evaluate individual or full set of nonlinear constraints.
            % ::
            %
            %   g = nln.eval(var, x)
            %   g = nln.eval(var, x, name)
            %   g = nln.eval(var, x, name, idx_list)
            %   [g, dg] = nln.eval(...)
            %
            % For a given value of the variable vector :math:`\x`, this method
            % evaluates the nonlinear constraint function and optionally its
            % derivatives for an individual subset, if name or name and index
            % list are provided, otherwise, for the full set of constraints.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   x (double) : full :math:`n_x \times 1` variable vector :math:`\x`
            %   name (char array) : name of subset/block of nonlinear
            %       constraints to evaluate
            %   idx_list (cell array) : *(optional)* index list for subset/block
            %       of nonlinear constraints to evaluate (for an indexed subset)
            %
            % Outputs:
            %   g (double) : :math:`n \times 1` constraint vector :math:`\g(\x)`
            %   dg (double) : :math:`n \times n_x` constraint Jacobian
            %       :math:`\g_\x = \der{\g}{\x}`
            %
            % See also add, params.

            %% initialize
            if nargin < 5
                idx = {};
            end

            if nargin < 4                       %% full set
                %% initialize g, dg
                g = NaN(obj.N, 1);
                if nargout > 1
                    dg = sparse(0, var.N);      %% build gradient by stacking
                end

                %% calls to substruct() are relatively expensive, so we pre-build the
                %% structs for addressing cell and numeric array fields, updating only
                %% the subscripts before use
                sc = struct('type', {'.', '{}'}, 'subs', {'', 1});  %% cell array field
                sn = sc; sn(2).type = '()';                         %% num array field

                %% fill in each piece
                for k = 1:obj.NS
                    name = obj.order(k).name;
                    idx  = obj.order(k).idx;
                    if isempty(idx)
                        if ~isfield(obj.data.fcn, name)
                            continue;   %% skip, there is no function handle stored here,
                                        %% the function value for this named set was included
                                        %% in the value computed by a previous named set
                        end
                        N = obj.idx.N.(name);   %% number of constraint functions
                                                %% evaluated for this named set
                        if isfield(obj.data.include, name)
                            N = N + sum(obj.data.include.(name).N);
                        end
                    else
                        % (calls to substruct() are relatively expensive ...
                        % sn = substruct('.', name, '()', idx);
                        % sc = substruct('.', name, '{}', idx);
                        % ... so replace them with these more efficient lines)
                        sn(1).subs = name;
                        sn(2).subs = idx;
                        sc(1).subs = name;
                        sc(2).subs = idx;
                        N = subsref(obj.idx.N, sn);
                    end
                    if N                                %% non-zero number of rows
                        if isempty(idx)
                            fcn = obj.data.fcn.(name);      %% fcn for kth constraint set
                            i1 = obj.idx.i1.(name);         %% starting row index
                            iN = i1 + N - 1;                %% ending row index
                            vs = obj.data.vs.(name);        %% var sets
                        else
                            fcn = subsref(obj.data.fcn, sc);%% fcn for kth constraint set
                            i1 = subsref(obj.idx.i1, sn);   %% starting row index
                            iN = subsref(obj.idx.iN, sn);   %% ending row index
                            vs = subsref(obj.data.vs, sc);  %% var sets
                        end
                        xx = var.varsets_x(x, vs);
                        if nargout < 2
                            gk = fcn(xx);           %% evaluate kth constraint w/o gradient
                        else
                            [gk, dgk] = fcn(xx);    %% evaluate kth constraint and gradient

                            if isempty(vs)          %% all rows of x
                                if size(dgk, 2) == var.N
                                    dg = [dg; dgk];
                                else                %% must have added vars since adding
                                                    %% this constraint set
                                    dg(i1:iN, 1:size(dgk, 2)) = dgk;
                                end
                            else                    %% selected rows of x
                                jj = var.varsets_idx(vs);   %% column indices for var set
                                dgi = sparse(N, var.N);
                                dgi(:, jj) = dgk;
                                dg = [dg; dgi];
                            end
                        end
                        g(i1:iN) = gk;          %% assign kth constraint
                    end
                end
            else                                %% individual named set
                if isempty(idx) && prod(size(obj.idx.i1.(name))) ~= 1
                    error('mp.sm_nln_constraint.eval: nonlinear constraint set ''%s'' requires an IDX_LIST arg', name)
                end
                [N, fcn, hess, vs] = obj.params(var, name, idx);
                xx = var.varsets_x(x, vs);
                if nargout < 2
                    g = fcn(xx);
                else
                    [g, dg] = fcn(xx);
                end
            end
        end

        function d2g = eval_hess(obj, var, x, lam, name, idx)
            % Evaluate "Hessian" for full set of nonlinear constraints.
            % ::
            %
            %   d2g = nln.eval_hess(var, x, lam)
            %
            %   % not yet implemented
            %   d2g = nln.eval_hess(var, x, lam, name)
            %   d2g = nln.eval_hess(var, x, lam, name, idx_list)
            %
            % For a given value of the variable vector :math:`\x`, this method
            % evaluates the nonlinear constraint "Hessian" for the full set of
            % constraints.
            %
            % Instead of evaluating the full 3 dimensional Hessian, it
            % actually evaluates the Jacobian of the vector formed by
            % multiplying the transpose of the constraint Jacobian by a
            % vector :math:`\lam` of multipliers.
            %
            % .. note:: Evaluation of Hessian for individual subsets not yet
            %   implemented.
            %
            % Inputs:
            %   var (mp.sm_variable) : corresponding mp.sm_variable object
            %   x (double) : full :math:`n_x \times 1` variable vector :math:`\x`
            %   lam (double) : :math:`n \times 1` vector of multipliers
            %       :math:`\lam`
            %   name (char array) : *(optional, and not yet implemented)* name
            %       of subset/block of nonlinear constraints to evaluate
            %   idx_list (cell array) : *(optional, and not yet implemente)*
            %       index list for subset/block of nonlinear constraints to
            %       evaluate (for an indexed subset)
            %
            % Outputs:
            %   d2g (double) : :math:`n_x \times n_x` "Hessian" matrix
            %       :math:`\g_{\x\x}(\lam) = \der{}{\x}(\trans{\g_\x} \lam)`
            %
            % See also add, params.

            %% initialize d2g (use transpose for speed on older versions of MATLAB)
            d2gt = sparse(var.N, var.N);

            %% calls to substruct() are relatively expensive, so we pre-build the
            %% structs for addressing cell and numeric array fields, updating only
            %% the subscripts before use
            sc = struct('type', {'.', '{}'}, 'subs', {'', 1});  %% cell array field
            sn = sc; sn(2).type = '()';                         %% num array field

            %% fill in each piece
            for k = 1:obj.NS
                name = obj.order(k).name;
                idx  = obj.order(k).idx;
                if isempty(idx)
                    if ~isfield(obj.data.hess, name)
                        continue;   %% skip, there is no function handle stored here,
                                    %% the function value for this named set was included
                                    %% in the value computed by a previous named set
                    end
                    N = obj.idx.N.(name);   %% number of constraint functions
                                            %% evaluated for this named set
                    if isfield(obj.data.include, name)
                        N = N + sum(obj.data.include.(name).N);
                    end
                else
                    % (calls to substruct() are relatively expensive ...
                    % sn = substruct('.', name, '()', idx);
                    % sc = substruct('.', name, '{}', idx);
                    % ... so replace them with these more efficient lines)
                    sn(1).subs = name;
                    sn(2).subs = idx;
                    sc(1).subs = name;
                    sc(2).subs = idx;
                    N = subsref(obj.idx.N, sn);
                end
                if N                                %% non-zero number of rows
                    if isempty(idx)
                        d2g_fcn = obj.data.hess.(name); %% Hessian fcn for kth constraint set
                        i1 = obj.idx.i1.(name);         %% starting row index
                        iN = i1 + N - 1;                %% ending row index
                        vs = obj.data.vs.(name);        %% var sets
                    else
                        d2g_fcn = subsref(obj.data.hess, sc);   %% Hessian fcn for kth constraint set
                        i1 = subsref(obj.idx.i1, sn);   %% starting row index
                        iN = subsref(obj.idx.iN, sn);   %% ending row index
                        vs = subsref(obj.data.vs, sc);  %% var sets
                    end
                    xx = var.varsets_x(x, vs);
                    d2gk = d2g_fcn(xx, lam(i1:iN));     %% evaluate kth Hessian

                    nk = size(d2gk, 2);
                    if isempty(vs)          %% all rows of x
                        if nk == var.N
                            d2gkt_full = d2gk';
                        else                %% must have added vars since adding
                                            %% this constraint set
                            d2gk_all_cols = sparse(nk, var.N);
                            d2gk_all_cols(:, 1:nk) = d2gk;
                            d2gkt_full = sparse(var.N, var.N);
                            d2gkt_full(:, 1:nk) = d2gk_all_cols';
                        end
                    else                    %% selected rows of x
                        jj = var.varsets_idx(vs);   %% indices for var set
                        d2gk_all_cols = sparse(nk, var.N);
                        d2gk_all_cols(:, jj) = d2gk;
                        d2gkt_full = sparse(var.N, var.N);
                        d2gkt_full(:, jj) = d2gk_all_cols';
                    end
                    d2gt = d2gt + d2gkt_full;
                end
            end
            d2g = d2gt';
        end

        function obj = display_soln(obj, var, soln, iseq, varargin)
            % Display solution values for nonlinear constraints.
            % ::
            %
            %   nln.display_soln(var, soln, iseq)
            %   nln.display_soln(var, soln, iseq, name)
            %   nln.display_soln(var, soln, iseq, name, idx_list)
            %   nln.display_soln(var, soln, iseq, fid)
            %   nln.display_soln(var, soln, iseq, fid, name)
            %   nln.display_soln(var, soln, iseq, fid, name, idx_list)
            %
            % Displays the solution values for all nonlinear constraints
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
            %               - ``lower`` - variable lower bounds
            %               - ``upper`` - variable upper bounds
            %   iseq (boolean) : true for equality constraints, false for
            %       inequality constraints
            %   fid (fileID) : fileID of open file to write to (default is
            %       1 for standard output)
            %   name (char array) : *(optional)* name of individual subset
            %   idx_list (cell array) : *(optional)* indices of individual
            %       subset

            [fid, name, idx, idxs, hdr1] = obj.display_soln_std_args(varargin{:});

            if obj.N
                if iseq
                    hdr2 = {'    val    lambda', ...
                            ' -------- --------' };
                    if isempty(soln.lambda)
                        lam = NaN(obj.N, 1);
                    else
                        lam = soln.lambda.eqnonlin;
                    end
                else
                    hdr2 = {'    val      ub      mu_ub', ...
                            ' -------- -------- --------' };
                    if isempty(soln.lambda)
                        mu_u = NaN(obj.N, 1);
                    else
                        mu_u = soln.lambda.ineqnonlin;
                    end
                end
                v = obj.eval(var, soln.x);

                %% print header rows
                obj.display_soln_print_headers(fid, hdr1, hdr2);

                %% print data
                none = '- ';
                for k = 1:length(idxs)
                    obj.display_soln_print_row(fid, idxs(k));
                    if iseq
                        if abs(lam(idxs(k))) < obj.mu_thresh()
                            mu_ub = sprintf(none);
                        else
                            mu_ub = obj.sprintf_num(8, lam(idxs(k)));
                        end
                        fprintf(fid, '%9s%9s\n', obj.sprintf_num(8, v(idxs(k))), ...
                            mu_ub);
                    else
                        if abs(mu_u(idxs(k))) < obj.mu_thresh()
                            mu_ub = sprintf(none);
                        else
                            mu_ub = obj.sprintf_num(8, mu_u(idxs(k)));
                        end
                        fprintf(fid, '%9s%9s%9s\n', obj.sprintf_num(8, v(idxs(k))), ...
                            '0', mu_ub);
                    end
                end

                %% print footer rows
                fprintf(fid, '%s\n', [hdr1{2} hdr2{2}]);
                if iseq
                    fprintf(fid, '%7s %-28s%9s%9s\n', '', 'Min', ...
                        obj.sprintf_num(8, min(v)), ...
                        obj.sprintf_num(8, min(lam)));
                    fprintf(fid, '%7s %-28s%9s%9s\n', '', 'Max', ...
                        obj.sprintf_num(8, max(v)), ...
                        obj.sprintf_num(8, max(lam)));
                else
                    fprintf(fid, '%7s %-28s%9s%9s%9s\n', '', 'Max', ...
                        obj.sprintf_num(8, max(v)), '0', ...
                        obj.sprintf_num(8, max(mu_u)));
                end
                fprintf(fid, '\n');
            end
        end

        function varargout = get_soln(obj, var, soln, iseq, varargin)
            % Fetch solution values for specific named/indexed subsets.
            % ::
            %
            %   vals = nln.get_soln(var, soln, iseq, name)
            %   vals = nln.get_soln(var, soln, iseq, name, idx_list)
            %   vals = nln.get_soln(var, soln, iseq, tags, name)
            %   vals = nln.get_soln(var, soln, iseq, tags, name, idx_list)
            %
            % Returns named/indexed nonlinear constraint results for a solved
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
            %   iseq (boolean) : true for equality constraints, false for
            %       inequality constraints
            %   tags (char array or cell array of char arrays) : names of
            %       desired outputs, default is ``{'g', 'lam', 'dg'}`` with
            %       valid values:
            %
            %           - ``'g'`` or ``'h'`` - constraint value :math:`\g(\x)`
            %           - ``'lam'`` or ``'mu'`` - shadow price :math:`\lam` on
            %             constraint
            %           - ``'dg'`` or ``'dh'`` - constraint Jacobian
            %             :math:`\g_\x = \der{\g}{\x}`
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
            %     [g, lam, dg] = nln.get_soln(var, soln, false, 'flow');
            %     dg_Pmis_5_3 = nln.get_soln(var, soln, true, 'dg', 'Pmis', {5,3});
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
                if any(ismember({'dg', 'dh'}, tags(1:nargout)))
                    [g, dg] = obj.eval(var, soln.x, name, idx);
                elseif ismember('g', tags(1:nargout))
                    g = obj.eval(var, soln.x, name, idx);
                end
                for k = 1:nargout
                    switch tags{k}
                        case {'g', 'h'}
                            varargout{k} = g;
                        case {'dg', 'dh'}
                            varargout{k} = dg;
                        case {'lam', 'mu'}
                            if iseq
                                varargout{k} = soln.lambda.eqnonlin(i1:iN);
                            else
                                varargout{k} = soln.lambda.ineqnonlin(i1:iN);
                            end
                        otherwise
                            error('mp.sm_nln_constraint.get_soln: unknown tag ''%s''', tags{k});
                    end
                end
            end     %% if N
        end

        function ps = parse_soln(obj, soln, iseq, stash)
            % Parse solution for nonlinear constraints.
            % ::
            %
            %   ps = nln.parse_soln(soln, iseq)
            %
            % Parse a full solution struct into parts corresponding to
            % individual nonlinear constraint subsets.
            %
            % Inputs:
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
            %               - ``lower`` - variable lower bounds
            %               - ``upper`` - variable upper bounds
            %   iseq (boolean) : true for equality constraints, false for
            %       inequality constraints
            %   stash (boolean) : if true, store return value in :attr:`soln`
            %       property
            %
            % Output:
            %   ps (struct) : parsed solution, struct where each field listed
            %       below is a struct whos names are the names of the relevant
            %       nonlinear constraint subsets and values are scalars for
            %       named sets, arrays for named/indexed sets:
            %
            %           - ``lam`` - equality constraint shadow prices
            %           - ``mu`` - inequality constraint shadow prices

            ps = []; params = [];
            if iseq
                if obj.get_N() && isfield(soln.lambda, 'eqnonlin')
                    params = struct('src', soln.lambda.eqnonlin, ...
                                    'dst',     'lam'  );
                end
            else
                if obj.get_N && isfield(soln.lambda, 'ineqnonlin')
                    params = struct('src', soln.lambda.ineqnonlin, ...
                                    'dst',      'mu' );
                end
            end
            if ~isempty(params)
                ps = obj.parse_soln_fields(params);
            end

            if nargin > 3 && stash
                obj.soln = ps;
            end
        end
    end     %% methods

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

            default_tags = {'g', 'lam', 'dg'};
        end
    end     %% methods (Access=protected)
end         %% classdef
