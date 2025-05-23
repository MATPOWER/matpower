classdef sm_variable < mp.set_manager_opt_model
% mp.sm_variable -  MP Set Manager class for variables.
% ::
%
%   var = mp.sm_variable()
%   var = mp.sm_variable(label)
%
% MP Set Manager class for variables. Manages variable initial values,
% lower and upper bounds, and variable type, along with indexing.
%
% By convention, ``var`` is the variable name used for mp.sm_variable objects.
%
% mp.sm_variable Properties:
%   * cache - struct for caching aggregated parameters for the set
%
% mp.sm_variable Methods:
%   * sm_variable - constructor
%   * add - add a subset of variables, with initial value, bounds, and var type
%   * params - return initial values, lower bounds, upper bounds, and var type
%   * set_params - modify parameter data for variables
%   * display_soln - display solution values for variables
%   * get_soln - fetch solution values for specific named/indexed subsets
%   * parse_soln - parse solution for variables
%   * varsets_idx - return vector of indices into full :math:`\x` corresponding to ``vs``
%   * varsets_len - return the total number of variables specified by ``vs``
%   * varsets_x - return subset of :math:`\x` specified by ``vs``
%   * varsets_cell2struct - convert ``vs`` from cell array to struct array
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
        % struct for caching aggregated parameters for variables
        cache = [];
    end     %% properties

    methods
        function obj = sm_variable(varargin)
            % Constructor.
            % ::
            %
            %   var = mp.sm_variable(label)

            es = struct();  %% empty struct
            obj@mp.set_manager_opt_model(varargin{:});
            obj.data = struct( ...
                'v0', es, ...
                'vl', es, ...
                'vu', es, ...
                'vt', es );
        end

        function obj = add(obj, name, idx, varargin)
            % Add a subset of variables with initial value, bounds, type.
            % ::
            %
            %   var.add(name, N, v0, vl, vu, vt)
            %   var.add(name, N, v0, vl, vu)
            %   var.add(name, N, v0, vl)
            %   var.add(name, N, v0)
            %   var.add(name, N)
            %
            %   var.add(name, idx_list, N, v0, vl, vu, vt)
            %   var.add(name, idx_list, N, v0, vl, vu)
            %   var.add(name, idx_list, N, v0, vl)
            %   var.add(name, idx_list, N, v0)
            %   var.add(name, idx_list, N)
            %
            % Add a named, and possibly indexed, subset of variables to
            % the set.
            %
            % Inputs:
            %   name (char array) : name of subset/block of variables to add
            %   idx_list (cell array) : *(optional)* index list for subset/block
            %       of variables to add (for an indexed subset)
            %   N (integer) : number of variables in the subset
            %   v0 (double) : *(optional, default = 0)* scalar or
            %       :math:`N \times 1` vector of variable initial values
            %   vl (double) : *(optional, default = -Inf)* scalar or
            %       :math:`N \times 1` vector of variable lower bounds
            %   vu (double) : *(optional, default = Inf)* scalar or
            %       :math:`N \times 1` vector of variable upper bounds
            %   vt (char array) : *(optional, default =* ``'C'`` *)* scalar or
            %       :math:`1 \times N` char array of variable types, where
            %       accepted types are:
            %
            %       - ``'C'`` - continuous
            %       - ``'I'`` - integer
            %       - ``'B'`` - binary
            %
            % Examples::
            %
            %   var.add('V', nb, V0, Vmin, Vmax, 'C');
            %
            %   var.init_indexed_name('x', {2, 3});
            %   for i = 1:2
            %       for j = 1:3
            %           var.add('x', {i, j}, nx(i,j), ...);
            %       end
            %   end
            %
            % See also params, set_params.

            %% call parent to handle standard indexing
            add@mp.set_manager_opt_model(obj, name, idx, varargin{:});

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
            v0 = []; vl = []; vu = []; vt = [];
            if nargs >= 1
                v0 = args{1};
                if nargs >= 2
                    vl = args{2};
                    if nargs >= 3
                        vu = args{3};
                        if nargs >= 4
                            vt = args{4};
                        end
                    end
                end
            end
            if isempty(v0)
                v0 = zeros(N, 1);   %% init to zero by default
            elseif N > 1 && length(v0) == 1     %% expand from scalar as needed
                v0 = v0 * ones(N, 1);
            end
            if isempty(vl)
                vl = -Inf(N, 1);    %% unbounded below by default
            elseif N > 1 && length(vl) == 1     %% expand from scalar as needed
                vl = vl * ones(N, 1);
            end
            if isempty(vu)
                vu = Inf(N, 1);     %% unbounded above by default
            elseif N > 1 && length(vu) == 1     %% expand from scalar as needed
                vu = vu * ones(N, 1);
            end
            if isempty(vt) && N > 0
                vt = 'C';           %% all continuous by default
            end

            %% assign data
            if isempty(idx)
                obj.data.v0.(name) = v0;        %% initial value
                obj.data.vl.(name) = vl;        %% lower bound
                obj.data.vu.(name) = vu;        %% upper bound
                obj.data.vt.(name) = vt;        %% variable type
            else
                %% calls to substruct() are relatively expensive, so we
                %% pre-build the struct for addressing cell array fields
                %% sc = substruct('.', name, '{}', idx);
                sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
                obj.data.v0 = subsasgn(obj.data.v0, sc, v0);    %% initial value
                obj.data.vl = subsasgn(obj.data.vl, sc, vl);    %% lower bound
                obj.data.vu = subsasgn(obj.data.vu, sc, vu);    %% upper bound
                obj.data.vt = subsasgn(obj.data.vt, sc, vt);    %% variable type
            end
            if ~isempty(obj.cache)  %% clear cache of aggregated params
                obj.cache = [];
            end
        end

        function [v0, vl, vu, vt] = params(obj, name, idx)
            % Return initial values, lower bounds, upper bounds and var type.
            % ::
            %
            %   [v0, vl, vu] = var.params()
            %   [v0, vl, vu] = var.params(name)
            %   [v0, vl, vu] = var.params(name, idx_list)
            %   [v0, vl, vu, vt] = var.params(...)
            %
            % Returns the initial value, lower bound, upper bound, and
            % variable type for the full set of variables, if called without
            % input arguments, or for a specific named or named and indexed
            % subset. Values for the full set are cached for subsequent calls.
            %
            % Inputs:
            %   name (char array) : *(optional)* name of subset
            %   idx_list (cell array) : *(optional)* index list for subset
            %
            % Outputs:
            %   v0 (double) : column vector of variable initial values
            %   vl (double) : column vector of variable lower bounds
            %   vu (double) : column vector of variable upper bounds
            %   vt (char array) : char array of variable types:
            %
            %       - ``'C'`` - continuous
            %       - ``'I'`` - integer
            %       - ``'B'`` - binary
            %
            % Examples::
            %
            %   [x0, xmin, xmax] = var.params();
            %   [Pg0, Pmin, Pmax] = var.params('Pg');
            %   [zij0, zijmin, zijmax, ztype] = var.params('z', {i, j});
            %
            % See also add, set_params.

            if nargout > 3
                have_vt = 1;
            else
                have_vt = 0;
            end
            if nargin < 2       %% aggregate
                if isempty(obj.cache)       %% build the aggregate
                    v0 = []; vl = []; vu = []; vt = char([]);
                    %% calls to substruct() are relatively expensive, so we pre-build the
                    %% structs for addressing cell and numeric array fields, updating only
                    %% the subscripts before use
                    sc = struct('type', {'.', '{}'}, 'subs', {'', 1});  %% cell array field
                    for k = 1:obj.NS
                        name = obj.order(k).name;
                        idx = obj.order(k).idx;
                        if isempty(idx)
                            v0 = [ v0; obj.data.v0.(name) ];
                            vl = [ vl; obj.data.vl.(name) ];
                            vu = [ vu; obj.data.vu.(name) ];
                            if have_vt
                                N = obj.idx.N.(name);
                                vt0 = obj.data.vt.(name);
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
                            v0 = [ v0; subsref(obj.data.v0, sc) ];
                            vl = [ vl; subsref(obj.data.vl, sc) ];
                            vu = [ vu; subsref(obj.data.vu, sc) ];
                            if have_vt
                                % (calls to substruct() are relatively expensive ...
                                % sn = substruct('.', name, '()', idx);
                                % ... so replace it with these more efficient lines)
                                sn = sc; sn(2).type = '()';
                                N = subsref(obj.idx.N, sn);
                                vt0 = subsref(obj.data.vt, sc);
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

                    %% cache aggregated parameters
                    obj.cache = struct('v0', v0, 'vl', vl, 'vu', vu, 'vt', vt);
                else                    %% return cached values
                    v0 = obj.cache.v0;
                    vl = obj.cache.vl;
                    vu = obj.cache.vu;
                    vt = obj.cache.vt;
                end
            else                %% individual set
                if isfield(obj.idx.N, name)
                    if nargin < 3 || isempty(idx)
                        v0 = obj.data.v0.(name);
                        vl = obj.data.vl.(name);
                        vu = obj.data.vu.(name);
                        if have_vt
                            N = obj.idx.N.(name);
                            vt0 = obj.data.vt.(name);
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
                        v0 = subsref(obj.data.v0, sc);
                        vl = subsref(obj.data.vl, sc);
                        vu = subsref(obj.data.vu, sc);
                        if have_vt
                            % (calls to substruct() are relatively expensive ...
                            % sn = substruct('.', name, '()', idx);
                            % ... so replace it with these more efficient lines)
                            sn = sc; sn(2).type = '()';
                            N = subsref(obj.idx.N, sn);
                            vt0 = subsref(obj.data.vt, sc);
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

        function obj = set_params(obj, name, idx, params, vals)
            % set_params - Modify parameter data for variables.
            % ::
            %
            %   var.set_params(name, params, vals)
            %   var.set_params(name, idx_list, params, vals)
            %
            % This method can be used to modify parameters for an existing
            % subset of variables.
            %
            % Inputs:
            %   name (char array) : name of subset/block of variables to modify
            %   idx_list (cell array) : *(optional)* index list for subset/block
            %       of variables to modify (for an indexed subset)
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
            % Valid parameter names are ``N``, ``v0``, ``vl``, ``vu``, ``vt``.
            %
            % .. note:: Changing the dimension of a variable subset is not
            %   allowed.
            %
            % Examples::
            %
            %   var.set_params('Pg', 'v0', Pg0);
            %   var.set_params('y', {2,3}, {'vl', 'vu'}, {yl, yu});
            %   var.set_params('Va', 'all', {N, va0, val, vau, vat});
            %
            % See also add, params.

            if nargin < 5
                vals = params;
                params = idx;
                idx = {};
            end

            %% create default list of parameters to update based on set type & inputs
            default_params = {'N', 'v0', 'vl', 'vu', 'vt'};

            %% standardize provided arguments in cell arrays params, vals
            [is_all, np, params, vals] = ...
                obj.set_params_std_args(default_params, params, vals);

            %% get current parameters
            [v0, vl, vu, vt] = obj.params(name, idx);
            N0 = obj.get_N(name, idx);
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
                error('mp.sm_variable.set_params: dimension change for ''%s'' not allowed', obj.nameidxstr(name, idx));
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
                            error('mp.sm_variable.set_params: parameter ''%s'' ''%s'' should have length %d (or 1)', obj.nameidxstr(name, idx), pn{1}, N);
                        end
                    end
                end
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

            %% clear cached parameters
            obj.cache = [];

            %% update dimensions and indexing, if necessary
            dN = N - N0;
            if is_all && dN
                obj.set_params_update_dims(dN, name, idx);
            end
        end

        function obj = display_soln(obj, soln, varargin)
            % Display solution values for variables.
            % ::
            %
            %   var.display_soln(soln)
            %   var.display_soln(soln, name)
            %   var.display_soln(soln, name, idx_list)
            %   var.display_soln(soln, fid)
            %   var.display_soln(soln, fid, name)
            %   var.display_soln(soln, fid, name, idx_list)
            %
            % Displays the solution values for all variables (default) or an
            % individual named or named/indexed subset.
            %
            % Inputs:
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
            %   idx_list (cell array) : *(optional)* indices of individual subset

            [fid, name, idx, idxs, hdr1] = obj.display_soln_std_args(varargin{:});

            if obj.N
                [v0, vl, vu] = obj.params();
                v = soln.x;
                if isempty(soln.lambda)
                    mu_l = NaN(size(v));
                    mu_u = mu_l;
                else
                    mu_l = soln.lambda.lower;
                    mu_u = soln.lambda.upper;
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

        function varargout = get_soln(obj, soln, varargin)
            % Fetch solution values for specific named/indexed subsets.
            % ::
            %
            %   vals = var.get_soln(soln, name)
            %   vals = var.get_soln(soln, name, idx_list)
            %   vals = var.get_soln(soln, tags, name)
            %   vals = var.get_soln(soln, tags, name, idx_list)
            %
            % Returns named/indexed variable results for a solved model,
            % evaluated at the solution found.
            %
            % Inputs:
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
            %       desired outputs, default is ``{'x', 'mu_l', 'mu_u'}``,
            %       with valid values:
            %
            %           - ``'x'`` - value of solution variable
            %           - ``'mu_l'`` - shadow price on variable lower bound
            %           - ``'mu_u'`` - shadow price on variable upper bound
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
            %     [P, muPmin, muPmax] = var.get_soln(soln, 'P');
            %     muRmin_2_3 = var.get_soln(soln, 'mu_l', 'R', {2,3});
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
                for k = 1:nargout
                    switch tags{k}
                        case 'x'
                            varargout{k} = soln.x(i1:iN);
                        case 'mu_l'
                            varargout{k} = soln.lambda.lower(i1:iN);
                        case 'mu_u'
                            varargout{k} = soln.lambda.upper(i1:iN);
                        otherwise
                            error('mp.sm_variable.get_soln: unknown tag ''%s''', tags{k});
                    end
                end
            end     %% if N
        end

        function ps = parse_soln(obj, soln, stash)
            % Parse solution for variables.
            % ::
            %
            %   ps = var.parse_soln(soln)
            %   ps = var.parse_soln(soln, stash)
            %
            % Parse a full solution struct into parts corresponding to
            % individual variable subsets.
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
            %               - ``lower`` - variable lower bounds
            %               - ``upper`` - variable upper bounds
            %   stash (boolean) : if true, store return value in :attr:`soln`
            %       property
            %
            % Output:
            %   ps (struct) : parsed solution, struct where each field listed
            %       below is a struct whos names are the names of the relevant
            %       variable subsets and values are scalars for named sets,
            %       arrays for named/indexed sets:
            %
            %           - ``val`` - variable values
            %           - ``mu_l`` - variable lower bound shadow prices
            %           - ``mu_u`` - variable upper bound shadow prices

            params = struct('src', soln.x, 'dst', 'val');
            if isfield(soln.lambda, 'lower')
                params(end+1).src = soln.lambda.lower;
                params(end  ).dst = 'mu_l';
            end
            if isfield(soln.lambda, 'upper')
                params(end+1).src = soln.lambda.upper;
                params(end  ).dst = 'mu_u';
            end
            ps = obj.parse_soln_fields(params);

            if nargin > 2 && stash
                obj.soln = ps;
            end
        end

        function kk = varsets_idx(obj, vs)
            % varsets_idx - Return vector of indices into full :math:`\x` corresponding to ``vs``.
            % ::
            %
            %   k = var.varsets_idx(vs)
            %
            % Returns a vector of indices into the full variable :math:`\x`
            % corresponding to the variable sub-vector specified by ``vs``.
            %
            % Input:
            %   vs (struct array) : variable set, struct array of ``name``,
            %       ``idx`` pairs of indexed named subsets of variables
            %
            % Output:
            %   k (integer) : vector of indices into full variable :math:`\x`
            %       for sub-vector specified by ``vs``.
            %
            % See also varsets_x.

            persistent sn;
            if isempty(vs)
                kk = (1:obj.N);
            else
                vsN = length(vs);
                k = cell(1, vsN);       %% indices for varsets

                %% calls to substruct() are relatively expensive, so we pre-build the
                %% struct for addressing numeric array fields, updating only
                %% the subscripts before use
                if isempty(sn)
                    sn = struct('type', {'.', '()'}, 'subs', {'', 1});
                end

                ki = 0;
                for v = 1:vsN
                    vname = vs(v).name;
                    vidx = vs(v).idx;
                    if isempty(vidx)
                        i1 = obj.idx.i1.(vname);        %% starting index in full x
                        iN = obj.idx.iN.(vname);        %% ending index in full x
                    else
                        % (calls to substruct() are relatively expensive ...
                        % sn = substruct('.', vname, '()', vidx);
                        % ... so replace it with these more efficient lines)
                        sn(1).subs = vname;
                        sn(2).subs = vidx;
                        i1 = subsref(obj.idx.i1, sn);   %% starting index in full x
                        iN = subsref(obj.idx.iN, sn);   %% ending index in full x
                    end
                    if isscalar(i1)         %% simple named set, or indexed named set
                        ki = ki + 1;
                        k{ki} = (i1:iN);                %% single set of indices for varset
                    else                    %% multi-dim named set w/no index specified
                        ii1 = permute(i1, ndims(i1):-1:1);
                        iiN = permute(iN, ndims(i1):-1:1);
                        for j = 1:length(ii1(:))
                            ki = ki + 1;
                            k{ki} = (ii1(j):iiN(j));    %% multiple sets of indices for varset
                        end
                    end
                end
                kk = [k{:}];
            end
        end

        function nv = varsets_len(obj, vs)
            % varsets_len - Return the total number of variables specified by ``vs``.
            % ::
            %
            %   nv = var.varsets_len(vs)
            %
            % Return the total number of elements in the variable sub-vector
            % specified by ``vs``.
            %
            % Input:
            %   vs (struct array) : variable set, struct array of ``name``,
            %       ``idx`` pairs of indexed named subsets of variables
            %
            % Output:
            %   nv (integer) : total number of elements in the variable
            %       sub-vector specified by ``vs``.
            %
            % See also varsets_cell2struct.

            persistent sn;
            if isempty(vs)
                nv = obj.N;
            else
                nv = 0;

                %% calls to substruct() are relatively expensive, so we pre-build the
                %% struct for addressing numeric array fields, updating only
                %% the subscripts before use
                if isempty(sn)
                    sn = struct('type', {'.', '()'}, 'subs', {'', 1});
                end

                for v = 1:length(vs)
                    idx = vs(v).idx;
                    if isempty(idx)
                        N = obj.idx.N.(vs(v).name);
                    else
                        % (calls to substruct() are relatively expensive ...
                        % sn = substruct('.', vs(v).name, '()', vs(v).idx);
                        % ... so replace it with these more efficient lines)
                        sn(1).subs = vs(v).name;
                        sn(2).subs = idx;
                        N = subsref(obj.idx.N, sn);
                    end
                    nv = nv + sum(N(:));
                end
            end
        end

        function xx = varsets_x(obj, x, vs, return_vector)
            % varsets_x - Return subset of :math:`\x` specified by ``vs``.
            % ::
            %
            %   xx = obj.varsets_x(x, vs)
            %   xx = obj.varsets_x(x, vs, return_vector)
            %
            % Returns a cell array of sub-vectors of :math:`\x` specified by
            % ``vs``, or a stacked variable vector.
            %
            % Input:
            %   x (double) : full variable vector :math:`\x`
            %   vs (struct array) : variable set, struct array of ``name``,
            %       ``idx`` pairs of indexed named subsets of variables
            %   return_vector : if present and true, returns a stacked vector
            %       instead of a cell array
            %
            % Output:
            %   xx (cell array or double) : cell array of sub-vectors of
            %       :math:`\x` specified by ``vs`` or vector of stacked
            %       sub-vectors, depending on presence and value of
            %       ``return_vector``; returns full stacked variable vector
            %       if ``vs`` is missing or empty.
            %
            % See also varsets_len.

            persistent sn;
            if isempty(vs)          %% all rows of x
                xx = x;
            else                    %% selected rows of x
                vsN = length(vs);
                xx = cell(1, vsN);

                %% calls to substruct() are relatively expensive, so we pre-build the
                %% struct for addressing numeric array fields, updating only
                %% the subscripts before use
                if isempty(sn)
                    sn = struct('type', {'.', '()'}, 'subs', {'', 1});
                end

                ki = 0;
                for v = 1:length(vs)
                    vname = vs(v).name;
                    vidx = vs(v).idx;
                    if isempty(vidx)
                        i1 = obj.idx.i1.(vname);    %% starting row in full x
                        iN = obj.idx.iN.(vname);    %% ending row in full x
                    else
                        % (calls to substruct() are relatively expensive ...
                        % sn = substruct('.', vname, '()', vidx);
                        % ... so replace it with these more efficient lines)
                        sn(1).subs = vname;
                        sn(2).subs = vidx;
                        i1 = subsref(obj.idx.i1, sn);   %% starting row in full x
                        iN = subsref(obj.idx.iN, sn);   %% ending row in full x
                    end
                    if isscalar(i1)         %% simple named set, or indexed named set
                        ki = ki + 1;
                        xx{ki} = x(i1:iN);  %% single set of indices for varset
                    else                    %% multi-dim named set w/no index specified
                        ii1 = permute(i1, ndims(i1):-1:1);
                        iiN = permute(iN, ndims(i1):-1:1);
                        for j = 1:length(ii1(:))
                            ki = ki + 1;
                            xx{ki} = x(ii1(j):iiN(j));  %% multiple sets of indices for varset
                        end
                    end
                end
                if nargin > 3 && any(return_vector)
                    xx = vertcat(xx{:});
                end
            end
        end
    end     %% methods

    methods (Access=protected)
        function default_tags = get_soln_default_tags(obj)
            % Return default tags for get_soln().
            % ::
            %
            %   default_tags = var.get_soln_default_tags()
            %
            % Output:
            %   default_tags (cell array) : tags defining the default outputs
            %       of get_soln(), namely ``{'x', 'mu_l', 'mu_u'}``
            %
            % See also get_soln.

            default_tags = {'x', 'mu_l', 'mu_u'};
        end
    end     %% methods (Access=protected)

    methods (Static)
        function vs = varsets_cell2struct(vs)
            % varsets_cell2struct - Convert ``vs`` from cell array to struct array.
            % ::
            %
            %   vs = mp.sm_variable.varsets_cell2struct(vs)
            %
            % Converts ``vs`` from a cell array to a struct array,
            % if necessary.
            %
            % Input:
            %   vs (cell or struct array) : variable set, cell array of names
            %       of variable subsets, or a struct array of ``name``, ``idx``
            %       pairs of indexed named subsets of variables
            %
            % Output:
            %   vs (struct array) : variable set, struct array of ``name``,
            %       ``idx`` pairs of indexed named subsets of variables
            %
            % See also varsets_len.

            %% convert varsets from cell to struct array if necessary
            if ~isempty(vs) && iscell(vs)
                empty_cells = cell(1, length(vs));
                [empty_cells{:}] = deal({});    %% empty cell arrays
                vs = struct('name', vs, 'idx', empty_cells);
            end
        end
    end     %% methods (Static)
end         %% classdef
