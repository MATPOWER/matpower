classdef set_manager < handle
% mp.set_manager -  MP Set Manager base class.
% ::
%
%   sm = mp.set_manager(label)
%
% Implements functionality to manage the indexing of various named and
% indexed blocks of elements, such as variables, constraints, etc. This class
% helps keep track of the ordering and indexing of the various blocks as they
% are added to the object. Subclasses may implement additional functionality
% to handle data associated with each block.
%
% mp.set_manager Properties:
%   * label - the label of the set, e.g. ``'VARIABLES'``, ``'LINEAR CONSTRAINTS'``
%   * idx - indexing information
%   * N - total number of entities in the full set
%   * NS - number of named or named/indexed subsets or blocks
%   * order - struct array of names/indices of blocks in order
%   * data - struct of additional set-type-specific data for each block
%
% mp.set_manager Methods:
%   * set_manager - constructor
%   * copy - make a duplicate (shallow copy) of the object
%   * to_struct - convert object data *to* a struct
%   * from_struct - copy object data *from* a struct
%   * add - add a named (and optionally indexed) subset of entities
%   * describe_idx - provide/display name and index label for given indices
%   * display - display summary of indexing of subsets in object
%   * get_N - return the number of elements in the set
%   * init_indexed_name - initialize dimensions for an indexed named set
%   * set_type_idx_map - map index back to named subset & index within set
%
% By convention, ``sm`` is the variable name used for generic mp.set_manager
% objects.
%
% **Simple vs. Indexed Named Subsets**
%
% A subset or block of entities can be added as a simple named subset
% identified by a single name (e.g. ``foo``), or as an indexed named subset
% identified by a name that is indexed by one or more indices (e.g.
% ``bar{4,3}``). For an indexed named subset, the dimensions of the indexed
% subset must be supplied by calling init_indexed_name() before calling add().
%
% **Example**
%
% Suppose the goal is to create and manage a set of variables related to the
% charging of a fleet of :math:`n` electric vehicles over a week, where we have
% variables representing daily charging amounts for each vehicle, daily
% discharge amounts from driving, and the final battery state at the end of
% the week. These could be organized as 7 blocks of ``charge`` variables,
% 1 for each day of the week, 7 blocks of ``discharge`` variables, and
% 1 ``battery_state`` block, where each block consists of an :math:`n \times 1`
% vector of variables.
%
% ::
%
%   n = 20;
%   sm = mp.set_manager('VARIABLES');
%   sm.init_indexed_name('charge', 7);
%   sm.init_indexed_name('discharge', 7);
%   for d = 1:7
%       sm.add('charge', {d}, n, ...);
%       sm.add('discharge', {d}, n, ...);
%   end
%   sm.add('battery_state', n, ...);
%
% The full set of variables will be a :math:`(2 \times 7 + 1) \times n` vector
% assembled in the order the variables are added, that is::
%
%   x = [
%       charge{1};
%       discharge{1};
%       charge{2};
%       discharge{2};
%           .
%           .
%           .
%       charge{7};
%       discharge{7};
%       battery_state ]

% v-----v  TO BE DELETED  v-----v
% The following is the structure of the data in the object. Each field of
% ``idx`` or ``data`` is a struct whose field names are the names of the
% corresponding blocks of elements of that type (e.g. variables, constraints,
% etc.). They are found in order in the corresponding .order field. The
% description next to these fields gives the meaning of the value for each
% named sub-field. E.g. var.data.v0.Pg contains a vector of initial values
% for the 'Pg' block of variables.
%
%   obj
%       .var        - data for 'var' set type, e.g. variable sets that
%                     make up the full optimization variable x
%           .idx
%               .i1 - starting index within x
%               .iN - ending index within x
%               .N  - number of elements in this variable set
%           .N      - total number of elements in x
%           .NS     - number of variable sets or named blocks
%           .data   - additional set-type-specific data for each block
%           .order  - struct array of names/indices for variable
%                     blocks in the order they appear in x
%               .name   - name of the block, e.g. Pg
%               .idx    - indices for name, {2,3} => Pg(2,3)
%       .userdata   - any user defined data
%           .(user defined fields)
% ^-----^  TO BE DELETED  ^-----^

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

    properties
        % *(char array)* label used as header for display
        label

        % *(struct)* indexing information, with the following 3 fields:
        %
        %   - ``i1`` -  starting index of subset within full set
        %   - ``iN`` -  ending index of subset within full set
        %   - ``N`` -  number of elements in this subset
        %
        % Each of the 3 fields above is a struct whose field names are the
        % names of the corresponding blocks of elements. They are found in
        % order in the corresponding field of the :attr:`order` property. The
        % value is a scalar (for simple named blocks) or an array (for
        % indexed named blocks).
        %
        % For example, for a simple block named ``foo``, ``idx.i1.foo``
        % is the starting index of the ``foo`` block in the full set. For
        % a ``bar`` blocked indexed by ``i`` and ``j``, ``idx.iN.bar(i,j)``
        % is the ending index of the ``bar{i,j}`` block in the full set.
        % Similarly, ``idx.N.bar(i,j)`` is the number of elements in that
        % block, i.e. ``idx.iN.bar(i,j) - idx.i1.bar(i,j) + 1``.
        idx = struct('i1', struct(), 'iN', struct(), 'N', struct());

        N = 0;  % *(integer)* total number of individual elements in set
        NS = 0; % *(integer)* total number of subsets in set

        % *(struct)* struct array of names/indices for blocks in the order
        % they appear in the full set, with fields:
        %
        %   - ``name`` - name of the block, e.g. ``bar``
        %   - ``idx`` - cell array of indices for the name, e.g. ``{2,3} ==> bar{2,3}``
        order = struct('name', [], 'idx', []);

        % *(struct)* additional set-type-specific data for each block
        data = struct();
    end     %% properties

    methods
        function obj = set_manager(label)
            % Constructor.
            % ::
            %
            %   sm = mp.set_manager(label)
            %   sm = mp.set_manager(s)
            %
            % Input:
            %   label (char array) : label used in display for this set type
            %   s (struct) : a struct to use as an input to from_struct()

            if isstruct(label)
                obj = obj.from_struct(label);
            elseif ischar(label)
                obj.label = label;
            else
                error('mp.set_manager: constructor input argument must be char array or struct');
            end
        end

        function new_obj = copy(obj, cls)
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
            % Output:
            %   new_sm (mp.set_manager) : duplicate of original object
            %
            % Make a shallow copy of the object by copying each of the
            % top-level properties.

            %% create new object
            if nargin < 2
                cls = class(obj);
            end
            new_obj = eval(sprintf('%s(''%s'')', cls, obj.label));

            %% get the names of the properties
            if have_feature('octave')
                s1 = warning('query', 'Octave:classdef-to-struct');
                warning('off', 'Octave:classdef-to-struct');
            end
            props = fieldnames(obj);
            if have_feature('octave')
                warning(s1.state, 'Octave:classdef-to-struct');
            end

            %% copy them
            for k = 1:length(props)
                new_obj.(props{k}) = obj.(props{k});
            end
        end

        function s =  to_struct(obj)
            % Convert object data *to* a struct.
            % ::
            %
            %   s = sm.to_struct()
            %
            % Converts the object data *to* a struct that can later be
            % converted back to an identical object using mp.struct2object.
            % Useful for saving the object data to a MAT-file in Octave.

            s = nested_struct_copy(struct(), obj);
            s.class_ = class(obj);
            s.constructor_args_ = { obj.label };
        end

        function obj = from_struct(obj, s)
            % Copy object data *from* a struct.
            % ::
            %
            %   sm.from_struct(s)
            %
            % Called by function mp.struct2object, after creating the object
            % to copy the object data *from* a struct. Useful for recreating
            % the object after loading struct data from a MAT-file in Octave.

            obj = nested_struct_copy(obj, s);
        end

        function obj = add(obj, name, idx, varargin)
            % Add a named (and optionally indexed) subset of entities.
            % ::
            %
            %   sm.add(name, N, ...)
            %   sm.add(name, idx_list, N, ...)
            %
            % This base class method handles the indexing part. Subclasses are
            % expected to override it to handle any data that goes with each
            % subset added for the given set type.
            %
            % For example::
            %
            %   % Variable Set
            %   sm_var.add(name, idx_list, N, v0, vl, vu, vt);
            %
            %   % Linear Constraint Set
            %   sm_lin.add(name, idx_list, N, A, l, u, varsets);
            %
            %   % Nonlinear Equality Constraint Set
            %   sm_nle.add(name, idx_list, N, fcn, hess, computed_by, varsets);
            %
            %   % Nonlinear Inequality Constraint Set
            %   sm_nli.add(name, idx_list, N, fcn, hess, computed_by, varsets);
            %
            %   % Quadratic Cost Set
            %   sm_qdc.add(name, idx_list, N, cp, varsets);
            %
            %   % General Nonlinear Cost Set
            %   sm_nlc.add(name, idx_list, N, fcn, varsets);

            %% handle input args
            if iscell(idx)
                N = varargin{1};
                args = varargin(2:end);
            else
                N = idx;
                idx = {};
                args = varargin;
            end

            %% add general indexing info about this named set
            if isempty(idx)     %% simple named set
                %% prevent duplicate name in set of specified type
                if isfield(obj.idx.N, name)
                    error('mp.set_manager.add: %s set named ''%s'' already exists', obj.label, name);
                end

                %% add indexing info about this set
                obj.idx.i1.(name)  = obj.N + 1;     %% starting index
                obj.idx.iN.(name)  = obj.N + N;     %% ending index
                obj.idx.N.(name)   = N;             %% number in set
                obj.N  = obj.idx.iN.(name);         %% number of elements of this type
                obj.NS = obj.NS + 1;                %% number of sets of this type
                obj.order(obj.NS).name = name;      %% add name to ordered list of sets
                obj.order(obj.NS).idx  = {};
            else                %% indexed named set
                %% calls to substruct() are relatively expensive, so we pre-build the
                %% struct for addressing numeric array fields
                %% sn = substruct('.', name, '()', idx);
                sn = struct('type', {'.', '()'}, 'subs', {name, idx});  %% num array field

                %% prevent duplicate name in set of specified type
                if subsref(obj.idx.i1, sn) ~= 0
                    str = '%d'; for m = 2:length(idx), str = [str ',%d']; end
                    nname = sprintf(['%s(' str, ')'], name, idx{:});
                    error('mp.set_manager.add: %s set named ''%s'' already exists', obj.label, nname);
                end

                %% add indexing info about this set
                obj.idx.i1  = subsasgn(obj.idx.i1, sn, obj.N + 1);  %% starting index
                obj.idx.iN  = subsasgn(obj.idx.iN, sn, obj.N + N);  %% ending index
                obj.idx.N   = subsasgn(obj.idx.N,  sn, N);          %% number in set
                obj.N  = subsref(obj.idx.iN, sn);   %% number of elements of this type
                obj.NS = obj.NS + 1;                %% number of sets of this type
                obj.order(obj.NS).name = name;      %% add name to ordered list of sets
                obj.order(obj.NS).idx  = idx;       %% add indices to ordered list of sets
            end
        end

        function label = describe_idx(obj, idxs)
            % Provide/display name and index label for given indices.
            % ::
            %
            %   label = sm.describe_idx(idxs)
            %
            % Returns char arrays describing (name and optionally index) the
            % element of the full set that corresponds to the indices in
            % ``idxs``. The return value is a char array if ``idxs`` is a
            % scalar, otherwise it is a cell array of char arrays of the same
            % dimension as ``idxs``.
            %
            % For example, this method can tell you that element 361
            % corresponds to element 5 of the ``foo`` block (``'foo(5)'``) or
            % to element 17 of the ``bar{8,5}`` block (``'bar{8,5}(17)'``)
            %
            % Input:
            %   idxs (integer) : scalar or array of indices into full set
            %
            % Output:
            %   label (char array or cell array) : label or cell array of
            %       same dimensions as ``idxs`` of labels containing the
            %       corresponding name (possibly indexed) and index within
            %       that name.
            %
            % Examples::
            %
            %   label = sm.describe_idx(87));
            %   labels = sm.describe_idx([38; 49; 93]));
            %
            % See also set_type_idx_map.

            label = cell(size(idxs));       %% pre-allocate return cell array

            if ~isempty(idxs)
                s = obj.set_type_idx_map(idxs);
            end
            for i = 1:length(idxs(:))
                idx = s(i).idx;
                if isempty(idx)
                    label{i} = sprintf('%s(%d)', s(i).name, s(i).i);
                else
                    if length(idx) <= 1
                        idxstr = sprintf('%d', idx{1});
                    else
                        idxstr = [sprintf('%d', idx{1}) sprintf(',%d', idx{2:end})];
                    end
                    label{i} = sprintf('%s{%s}(%d)', s(i).name, idxstr, s(i).i);
                end
            end
            if isscalar(idxs)               %% return scalar
                label = label{1};
            end
        end

        function display(obj, tag)
            % Display summary of indexing of subsets in object.
            %
            % This method is called automatically when omitting a semicolon
            % on a line that retuns an object of this class.
            %
            % Displays the details of the indexing for each subset or
            % block, as well as the total number of elements and subsets.

            if obj.NS
                fmt = '%-26s %6s %8s %8s %8s\n';
                fprintf(fmt, obj.label, 'name', 'i1', 'iN', 'N');
                fprintf(fmt, repmat('=', 1, length(obj.label)), '------', '-----', '-----', '------');
                idx = obj.idx;
                fmt = '%10d:%22s %8d %8d %8d\n';
                for k = 1:obj.NS
                    name = obj.order(k).name;
                    if isempty(obj.order(k).idx)
                        fprintf(fmt, k, name, idx.i1.(name), idx.iN.(name), idx.N.(name));
                    else
                        vsidx = obj.order(k).idx;
                        str = '%d'; for m = 2:length(vsidx), str = [str ',%d']; end
                        s = substruct('.', name, '()', vsidx);
                        nname = sprintf(['%s{' str, '}'], name, vsidx{:});
                        fprintf(fmt, k, nname, ...
                                subsref(idx.i1, s), subsref(idx.iN, s), subsref(idx.N, s));
                    end
                end
                if nargin < 2
                    tag = '';
                else
                    tag = [tag '.'];
                end
                fmt = sprintf('%%10d = %%sNS%%%dd = %%sN\\n\\n', 36-length(tag));
                fprintf(fmt, obj.NS, tag, obj.N, tag);
            else
                fprintf('%-26s  :  <none>\n', obj.label);
            end
        end

        function N = get_N(obj, name, idx)
            % Return the number of elements in the set.
            % ::
            %
            %   N = sm.get_N()
            %   N = sm.get_N(name)
            %   N = sm.get_N(name, idx_list)
            %
            % Returns either the total number of elements in the set or the
            % number corresponding to a specified named block, or indexed named
            % subset.
            %
            % Examples::
            %
            %   N = sm.get_N()             % total number of elements in set
            %   N = sm.get_N(name)         % # of elements in named set
            %   N = sm.get_N(name, idx)    % # of elements in indexed named set

            if nargin < 2
                N = obj.N;
            else
                if isfield(obj.idx.N, name)
                    if nargin < 3 || isempty(idx)
                        N = obj.idx.N.(name);
                    else
                        % s1 = substruct('.', name, '()', idx);
                        sn = struct('type', {'.', '()'}, 'subs', {name, idx});  %% num array field
                        N = subsref(obj.idx.N, sn);
                    end
                else
                    N = 0;
                end
            end
        end

        function obj = init_indexed_name(obj, name, dim_list)
            % Initialize dimensions for an indexed named set.
            % ::
            %
            %   sm.init_indexed_name(name, dim_list)
            %
            % A subset or block can be identified by a single ``name``, such
            % as ``foo``, or by a ``name`` that is indexed by one or more
            % indices, such as ``'bar{3,4}'``. For an indexed named set, before
            % adding the indexed subsets themselves, the dimensions of the
            % indexed set of names must be initialized by calling this method.
            %
            % Inputs:
            %   name (char array) : name of the block or subset
            %   dim_list (cell array) : dimensions of the indexing of the name
            %
            % Example::
            %
            %   % elements with indexed named set 'R{i,j}'
            %   sm.init_indexed_name('R', {2, 3});
            %   for i = 1:2
            %       for j = 1:3
            %           sm.add('R', {i, j}, ...);
            %       end
            %   end

            %% prevent duplicate name in set of specified type
            if isfield(obj.idx.N, name)
                error('mp.set_manager.init_indexed_name: %s set named ''%s'' already exists', ...
                    st_label, name);
            end

            %% use column vector if single dimension
            if length(dim_list) == 1
                dim_list = {dim_list{:}, 1};
            end

            %% add general info about this named set
            zero_vector = zeros(dim_list{:});
            obj.idx.i1.(name) = zero_vector;    %% starting index
            obj.idx.iN.(name) = zero_vector;    %% ending index
            obj.idx.N.(name)  = zero_vector;    %% number of entities

            %% initialize set-type-specific data
            fn = fieldnames(obj.data);
            empty_cell  = cell(dim_list{:});
            for k = 1:length(fn)
                obj.data.(fn{k}).(name) = empty_cell;
            end
        end

        function s = set_type_idx_map(obj, idxs, group_by_name)
            % Map index back to named subset & index within set.
            % ::
            %
            %   s = sm.set_type_idx_map()
            %   s = sm.set_type_idx_map(idxs)
            %   s = sm.set_type_idx_map(idxs, group_by_name)
            %
            % Returns a struct of same dimensions as ``idxs`` specifying, for
            % each index, the corresponding named set and element within the
            % named set.
            %
            % Inputs:
            %   idxs (integer) : *(optional)* scalar or array of indices into
            %       full set *(default, if empty or not provided, is*
            %       ``[1:ns]'`` *where* ``ns`` *is the full dimension of the
            %       set corresponding to the all elements)*
            %   group_by_name (logical) : *(default = false)* if true, then the
            %       results are consolidated, with a single entry in ``s`` for
            %       each unique name/idx pair, where the ``i`` and ``j`` fields
            %       are vectors. In this case ``s`` is 1 dimensional.
            %
            % Output:
            %   s (struct) : return struct with the following fields:
            %
            %       - ``name`` : name of corresponding set
            %       - ``idx`` : cell array of indices for the name, if named
            %         set is indexed
            %       - ``i`` : index of element within the set
            %       - ``j`` : *(only if* ``group_by_name == 1`` *)*,
            %         corresponding index of full set, equal to a particular
            %         element of ``idxs``
            %
            % Examples::
            %
            %     s = var.set_type_idx_map(87));
            %     s = lin.set_type_idx_map([38; 49; 93]));
            %     s = var.set_type_idx_map());
            %     s = lin.set_type_idx_map([], 1));
            %
            % See also describe_idx.

            %% default args
            if nargin < 3
                group_by_name = 0;
                if nargin < 2
                    idxs = [];
                end
            end

            NS = obj.NS;

            %% special case : everything and grouped by name
            if group_by_name && isempty(idxs)
                %% pre-allocate return struct
                c = cell(NS, 1);    %% pre-allocate return cell array
                s = struct('name', c, 'idx', c, 'i', c);
                for k = 1:NS;
                    name = obj.order(k).name;   %% name of individual set
                    idx = obj.order(k).idx;     %% idx of individual set
                    s(k).name = name;
                    if isempty(idx)
                        i1 = obj.idx.i1.(name);
                        iN = obj.idx.iN.(name);
                        N = obj.idx.N.(name);
                        s(k).i = [1:N]';
                        s(k).j = [i1:iN]';
                    else
                        s(k).idx = idx;
                        ss = substruct('.', name, '()', idx);
                        i1 = subsref(obj.idx.i1, ss);
                        iN = subsref(obj.idx.iN, ss);
                        N = subsref(obj.idx.N, ss);
                        s(k).i = [1:N]';
                        s(k).j = [i1:iN]';
                    end
                end
            else    %% general case
                if isempty(idxs)
                    idxs = [1:obj.N]';      %% all indices
                end

                %% check for invalid idxs
                if any(idxs > obj.N)
                    error('mp.set_manager.set_type_idx_map: IDXS must not exceed maximum index (%d)', obj.N);
                end
                if any(idxs < 1)
                    error('mp.set_manager.set_type_idx_map: IDXS must be positive');
                end

                %% pre-allocate return struct
                c = cell(size(idxs));       %% pre-allocate return cell array
                s = struct('name', c, 'idx', c, 'i', c);

                %% sort idxs so we can go through loop only once
                [sorted_idxs, jj] = sort(idxs(:), 'descend');
                k = NS;                         %% index into set type (decrementing)
                name = obj.order(k).name;       %% name of individual set
                idx = obj.order(k).idx;         %% idx of individual set
                for i = 1:length(sorted_idxs)   %% index into sorted_idxs
                    ii = sorted_idxs(i);        %% idx of interest
                    j = jj(i);                  %% index into s and original idxs
                    while k > 0
                        if isempty(idx)
                            i1 = obj.idx.i1.(name);
                            if ii >= i1
                                s(j).name = name;
                                s(j).i = ii - i1 + 1;
                                break;
                            else
                                k = k - 1;
                                name = obj.order(k).name;
                                idx = obj.order(k).idx;
                            end
                        else
                            ss = substruct('.', name, '()', idx);
                            i1 = subsref(obj.idx.i1, ss);
                            if ii >= i1
                                s(j).name = name;
                                s(j).i = ii - i1 + 1;
                                s(j).idx = idx;
                                break;
                            else
                                k = k - 1;
                                name = obj.order(k).name;
                                idx = obj.order(k).idx;
                            end
                        end
                    end
                end

                %% consolidate indices into vectors for each unique
                %% name/idx pair, if requested
                if group_by_name
                    %% extract fields
                    [name, idx, i] = deal(cell(size(idxs)));
                    [name{:}] = deal(s.name);
                    [idx{:}]  = deal(s.idx);
                    [i{:}]    = deal(s.i);      i = cell2mat(i);

                    %% find unique name/idx
                    name_idx = cellfun(@join_name_idx, name, idx, ...
                        'UniformOutput', 0);
                    [c, ia, ic] = unique(name_idx);

                    %% recreate struct, grouped by name/idx
                    c0 = cell(size(c));
                    s = struct('name', name(ia), 'idx', idx(ia), 'i', c0, 'j', c0);
                    for k = 1:length(ia)
                        s(k).i  = i(ic == k);
                        s(k).j = idxs(ic == k);
                    end
                end
            end
        end
    end     %% methods

    methods (Access=protected)
        function str = nameidxstr(obj, name, idx)
            str = sprintf('%s%s', name, idxstr(idx));
        end
    end     %% methods (Access=protected)
end         %% classdef


function name_idx = join_name_idx(name, idx)
    if isempty(idx)
        name_idx = name;
    else
        name_idx = [name sprintf('_%d', idx{:})];
    end
end

function str = idxstr(idx)
    if isempty(idx)
        str = '';
    elseif length(idx) == 1
        str = sprintf('(%d)', idx{1});
    else
        str = ['(' sprintf('%d', idx{1}) sprintf(',%d', idx{2:end}) ')'];
    end
end
