classdef set_manager_opt_model < mp.set_manager
% mp.set_manager_opt_model -  MP Set Manager base class for opt_model fields.
% ::
%
%   sm = mp.set_manager_opt_model(label)
%
% Implements functionality to handle parameters and solution data for
% set types used to implement properties of the opt_model class.
%
% mp.set_manager_opt_model Properties:
%   * soln - struct for storing parsed solution values
%
% mp.set_manager_opt_model Methods:
%   * params - *(abstract)* return set-type-specific parameter data
%   * set_params - *(abstract)* modify set-type-specific parameter data
%   * display_soln - display solution values
%   * has_parsed_soln - return true if parsed solution is available
%
% By convention, ``sm`` is the variable name used for mp.set_manager_opt_model
% objects.

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

    properties
        % *(struct)* for storing parsed solution values
        soln
    end     %% properties

    methods
        function obj = set_manager_opt_model(varargin)
            % Constructor.
            % ::
            %
            %   sm = mp.set_manager_opt_model(label)

            obj@mp.set_manager(varargin{:});
        end

        function rv = params(obj, name, idx)
            % Return set-type-specific parameters.
            % ::
            %
            %   [...] = sm.params()
            %   [...] = sm.params(name)
            %   [...] = sm.params(name, idx_list)
            %
            % .. note:: This abstract method must be implemented by a
            %   subclass.
            %
            % Returns set-type-specific parameters for the full set, if called
            % without input arguments, or for a specific named or named and
            % indexed subset.
            %
            % Inputs:
            %   name (char array) : *(optional)* name of subset
            %   idx_list (cell array) : *(optional)* index list for subset
            %
            % Outputs are determined by the implementing subclass.
            %
            % See also mp.set_manager.add, set_params.
        end

        function obj = set_params(obj, name, idx, params, vals)
            % Modify parameter data.
            % ::
            %
            %   sm.set_params(name, params, vals)
            %   sm.set_params(name, idx_list, params, vals)
            %
            % .. note:: This abstract method must be implemented by a
            %   subclass.
            %
            % This method can be used to modify set-type-specific parameters
            % for an existing subset.
            %
            % Inputs:
            %   name (char array) : name of subset/block of entities to modify
            %   idx_list (cell array) : *(optional)* index list for subset/block
            %       of entities modify (for an indexed subset)
            %   params : can be one of three options:
            %
            %       - ``'all'`` - indicates that ``vals`` is a cell array
            %         whose elements correspond to the input parameters of
            %         the :meth:`add() <mp.set_manager.add>` method
            %       - name of a parameter - ``val`` is the value of that
            %         parameter
            %       - cell array of parameter names - ``vals`` is a cell array
            %         of corresponding values
            %   vals : new value or cell array of new values corresponding to
            %       ``params``
            %
            % Valid parameter names are defined by the implementing subclass.
            %
            % See also mp.set_manager.add, params.
        end

        function obj = display_soln(obj, var, soln, varargin)
            % Display solution values for generic set type.
            % ::
            %
            %   sm.display_soln(var, soln)
            %   sm.display_soln(var, soln, name)
            %   sm.display_soln(var, soln, name, idx_list)
            %   sm.display_soln(var, soln, fid)
            %   sm.display_soln(var, soln, fid, name)
            %   sm.display_soln(var, soln, fid, name, idx_list)
            %
            % Displays the solution values for all elements (default)
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
                %% print header rows
                hdr2 = {'', ...
                        '' };
                obj.display_soln_print_headers(fid, hdr1, hdr2);

                %% print data
                none = '- ';
                for k = 1:length(idxs)
                    obj.display_soln_print_row(fid, idxs(k));
                    fprintf(fid, '\n');
                end

                %% print footer rows
                fprintf(fid, '%s\n', [hdr1{2} hdr2{2}]);
                fprintf(fid, '\n');
            end
        end

        function TorF = has_parsed_soln(obj)
            % Return true if parsed solution is available.
            % ::
            %
            %   TorF = sm.has_parsed_soln()
            %
            % Output:
            %   TorF (boolean) : true if parsed solution is available in
            %       :attr:`soln` property; format of :attr:`soln` depends
            %       on implementing subclass

            TorF = ~isempty(obj.soln);
        end

        function ps = parse_soln(obj, soln, stash)
            % Parse solution.
            % ::
            %
            %   ps = sm.parse_soln(soln)
            %
            % Parse a full solution struct into parts corresponding to
            % individual subsets.
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
            %   ps (struct) : parsed solution, struct where each field is a
            %       struct whos names are the names of the relevant subsets
            %       and values are scalars for named sets, arrays for
            %       named/indexed sets; fields depend on implementing class

        end
    end     %% methods

    methods (Access=protected)
        function [is_all, np, params, vals] = set_params_std_args(obj, default_params, params, vals)
            % Standardize input args for use in subclass set_params method.
            % ::
            %
            %   [is_all, np, params, vals] = sm.set_params_std_args(default_params, params, vals)
            %
            % See also set_params.

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
        end

        function obj = set_params_update_dims(obj, dN, name, idx)
            % Update parameter dimensions for use in subclass set_params method.
            % ::
            %
            %   sm.set_params_update_dims(dN, name, idx_list)
            %
            % See also set_params.

            %% calls to substruct() are relatively expensive, so we pre-build the
            %% struct for addressing num array fields
            %% sn = substruct('.', name, '()', idx);
            sn = struct('type', {'.', '()'}, 'subs', {'', 1});  %% num array field
            update = 0;             %% not yet reached set being updated
            update_i1 = 0;          %% flag to indicate whether to update i1
            for k = 1:obj.NS
                o = obj.order(k);
                if ~update && strcmp(o.name, name) && isequal(o.idx, idx)
                    update = 1;     %% arrived at set being updated
                end
                if update
                    if isempty(o.idx)   %% simple named set
                        if update_i1
                            obj.idx.i1.(o.name) = obj.idx.i1.(o.name) + dN;
                        else
                            obj.idx.N.(o.name) = obj.idx.N.(o.name) + dN;
                        end
                        obj.idx.iN.(o.name) = obj.idx.iN.(o.name) + dN;
                    else                %% indexed named set
                        sn(1).subs = o.name;
                        sn(2).subs = o.idx;
                        if update_i1
                            v = subsref(obj.idx.i1, sn);
                            obj.idx.i1 = subsasgn(obj.idx.i1, sn, v + dN);
                        else
                            v = subsref(obj.idx.N, sn);
                            obj.idx.N = subsasgn(obj.idx.N, sn, v + dN);
                        end
                        v = subsref(obj.idx.iN, sn);
                        obj.idx.iN = subsasgn(obj.idx.iN, sn, v + dN);
                    end
                    update_i1 = 1;  %% update i1 from here on out
                end
            end
            obj.N = obj.N + dN;
        end

        function [fid, name, idx, idxs, hdr1] = display_soln_std_args(obj, varargin)
            % Standardize input args for use in subclass display_soln method.
            % ::
            %
            %   [fid, name, idx, idxs, hdr1] = sm.display_soln_std_args()
            %   [...] = sm.display_soln_std_args(name)
            %   [...] = sm.display_soln_std_args(name, idx_list)
            %   [...] = sm.display_soln_std_args(fid)
            %   [...] = sm.display_soln_std_args(fid, name)
            %   [...] = sm.display_soln_std_args(fid, name, idx_list)
            %
            % Inputs:
            %   fid (fileID) : fileID of open file to write to (default is
            %       1 for standard output)
            %   name (char array) : *(optional)* name of individual subset
            %   idx_list (cell array) : *(optional)* indices of individual
            %       subset
            %
            % Outputs:
            %   fid (fileID) : fileID of open file to write to (default is
            %       1 for standard output)
            %   name (char array) : *(optional)* name of individual subset
            %   idx (cell array) : *(optional)* indices of individual subset
            %   idxs (integer) : vector of indices of individual elements
            %   hdr1 (cell array) : 2 element cell array char arrays; first is
            %       column headings, second is column separator row; for
            %       ``idx`` and ``description`` columns
            %
            % See also display_soln.

            if nargin < 2 || ischar(varargin{1})
                fid = 1;
                args = varargin;
            else
                fid = varargin{1};
                args = varargin(2:end);
            end
            nargs = length(args);

            idx = [];
            name = [];
            if nargs >= 1
                name = args{1};
                if nargs >= 2
                    idx = args{2};
                end
            end

            if obj.N
                if isempty(name)            %% all indices for set type
                    idxs = (1:obj.N);
                elseif isempty(idx)         %% all indices for set type & name
                    idxs = [];
                    for k = 1:length(obj.order)
                        if strcmp(obj.order(k).name, name)
                            i1 = obj.idx.i1.(name)(obj.order(k).idx{:});
                            iN = obj.idx.iN.(name)(obj.order(k).idx{:});
                            idxs = [idxs (i1:iN)];
                        end
                    end
                else                        %% indices for name, idx
                    idxs = (obj.idx.i1.(name)(idx{:}):obj.idx.iN.(name)(idx{:}));
                end
                hdr1 = {'  idx    description                ', ...
                        '------- ----------------------------' };
            else
                idxs = [];
                hdr1 = {};
            end
        end

        function obj = display_soln_print_headers(obj, fid, hdr1, hdr2)
            % Print headers for solution display.
            % ::
            %
            %   obj.display_soln_print_headers(fid, hdr2)
            %
            % Inputs:
            %   fid (fileID) : fileID of open file to write to (1 for standard
            %       output)
            %   hdr1 (cell array) : 2 element cell array char arrays; first is
            %       column headings, second is column separator row; for
            %       ``idx`` and ``description`` columns
            %   hdr2 (cell array) : 2 element cell array char arrays; first is
            %       column headings, second is column separator row; for
            %       set-type-specific columns
            %
            % See also display_soln.

            fprintf(fid, '=====  %s  =====\n', obj.label);
            for h = 1:length(hdr1)
                fprintf(fid, '%s\n', [hdr1{h} hdr2{h}]);
            end
        end

        function obj = display_soln_print_row(obj, fid, i)
            % Print idx and description columns for row i of solution display.
            % ::
            %
            %   obj.display_soln_print_row(fid, i)
            %
            % Inputs:
            %   fid (fileID) : fileID of open file to write to (1 for standard
            %       output)
            %   i (integer) : index of element to be printed in this row
            %
            % See also display_soln.

            ii = sprintf('%d', i);
            fmt = sprintf('%%-%ds', length(ii)+ceil((7-length(ii))/2));
            fprintf(fid, '%7s %-28s', sprintf(fmt, ii), obj.describe_idx(i));
        end

        function default_tags = get_soln_default_tags(obj)
            % Return default tags for use in get_soln method.
            % ::
            %
            %   default_tags = sm.get_soln_default_tags()
            %
            % .. note:: This protected abstract method must be implemented by
            %    a subclass.
            %
            % Output:
            %   default_tags (cell array) : tags defining the default outputs
            %       of get_soln()
            %
            % See also get_soln.

            default_tags = {};
        end

        function [tags, name, idx, N, i1, iN] = get_soln_std_args(obj, tags, name, idx)
            % Standardize input args for use in subclass get_soln method.
            % ::
            %
            %   [tags, name, idx_list, N, i1, iN] = sm.get_soln_std_args(tags, name, idx_list)
            %
            % See also get_soln.

            %% input arg handling
            if nargin == 2              %% obj.get_soln(soln, name)
                idx = [];
                name = tags;
                tags = {};
            elseif nargin == 3
                if ischar(name)         %% obj.get_soln(soln, tags, name)
                    idx = [];
                else                    %% obj.get_soln(soln, name, idx)
                    idx = name;
                    name = tags;
                    tags = {};
                end
            end

            %% set up tags for default outputs
            if isempty(tags)
                tags = obj.get_soln_default_tags();
            elseif ~iscell(tags)
                tags = { tags };
            end

            %% set up indexing
            if isempty(idx)         %% name, no idx provided
                N = obj.idx.N.(name);
                i1 = obj.idx.i1.(name);         %% starting row index
                iN = obj.idx.iN.(name);         %% ending row index
                if ~isscalar(N)
                    error('mp.set_manager_opt_model.get_soln_std_args: %s set ''%s'' requires an IDX_LIST arg', class(obj), name);
                end
            else                    %% indexed named set
                %% calls to substruct() are relatively expensive, so we pre-build the
                %% structs for addressing cell and numeric array fields, updating only
                %% the subscripts before use
                sn = struct('type', {'.', '()'}, 'subs', {name, idx});  %% num array field
                N = subsref(obj.idx.N, sn);
                i1 = subsref(obj.idx.i1, sn);   %% starting row index
                iN = subsref(obj.idx.iN, sn);   %% ending row index
            end
        end

        function ps = parse_soln_fields(obj, params)
            % Parse solution fields for subclass parse_soln method.
            % ::
            %
            %   ps = sm.parse_soln_fields(params)
            %
            % Parse solution fields.
            %
            % Input:
            %   params (struct) : struct array with fields:
            %
            %           - ``src`` - values from full solution struct
            %           - ``dst`` - name of destination field in parsed
            %             solution struct
            %
            % Output:
            %   ps (struct) : parsed solution struct, fields names match
            %       those provided in ``params.dst``, values are structs
            %       whose field names correspond to the named subsets in the
            %       set

            %% calls to substruct() are relatively expensive, so we pre-build the
            %% structs for addressing cell and numeric array fields, updating only
            %% the subscripts before use
            persistent sn;
            persistent sc;
            if isempty(sc)
                sc = struct('type', {'.', '{}'}, 'subs', {'', 1});  %% cell array field
            end
            if isempty(sn)
                sn = struct('type', {'.', '()'}, 'subs', {'', 1});  %% num array field
            end

            ps = struct();      %% parsed solution

            np = length(params);
            have_param = zeros(np, 1);
            for j = 1:np
                have_param(j) = ~isempty(params(j).src);
            end
            for k = 1:obj.NS
                name = obj.order(k).name;
                idx = obj.order(k).idx;
                if isempty(idx)
                    N = obj.idx.N.(name);
                else
                    sn(1).subs = name;
                    sn(2).subs = idx;
                    N = subsref(obj.idx.N, sn);
                    need_init = all([idx{:}] == 1);
                end
                if N
                    for j = 1:np
                        if have_param(j)    %% parameter is available
                            dname = params(j).dst;  %% destination field name
                            if isempty(idx)
                                i1 = obj.idx.i1.(name);
                                iN = obj.idx.iN.(name);
                                ps.(dname).(name)  = params(j).src(i1:iN);
                            else
                                if need_init
                                    param_names = fieldnames(obj.data);
                                    ps.(dname).(name) = cell(size(obj.data.(param_names{1}).(name)));
                                end
                                i1 = subsref(obj.idx.i1, sn);    %% starting row index
                                iN = subsref(obj.idx.iN, sn);    %% ending row index
                                sc(1).subs = name;
                                sc(2).subs = idx;
                                ps.(dname) = ...
                                    subsasgn(ps.(dname), sc, params(j).src(i1:iN));
                            end
                        end
                    end     %% for j
                end
            end     %% for k
        end

        function str = sprintf_num(obj, width, val)
            val = full(val);
            if all(isnan(val))
                fmt = sprintf('%%%ds', width);
                str = sprintf(fmt, ' ');
            else
                fmt = sprintf('%%%dg', width);
                str = sprintf(fmt, val);
            end
            if length(str) > width
                fmt = sprintf('%%%d.*g', width);
                for p = width-2:-1:0
                    str = sprintf(fmt, p, val);
                    if length(str) <= width
                        break;
                    end
                end
            end
            assert(length(str) <= width)
        end

        function v = mu_thresh(obj)
            v = 1e-7;
        end

        function v = num_inf(obj)
            v = 1e10;
        end
    end     %% methods (Access=protected)
end         %% classdef
