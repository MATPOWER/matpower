classdef mp_idx_manager < handle
% mp_idx_manager - |MATPOWER| Index Manager abstract class
%
% .. note::
%    .. deprecated:: 5.0 Please use mp.set_manager instead. The functionality
%       previously implemented in mp_idx_manager, a parent container class
%       for managing various set types, has been moved to mp.set_manager,
%       a base class for the set type objects themselves.
%
% ::
%
%   A MATPOWER Index Manager object can be used to manage the indexing of
%   various named and indexed blocks of various set types, such as variables,
%   constraints, etc. This class helps keep track of the ordering and
%   indexing of the various blocks as they are added to the object.
%
%   The types of named sets to be managed by the class are defined by the
%   DEF_SET_TYPES method, which assigns a struct to the 'set_types' field.
%
%   Properties
%       set_types   - a struct defined by DEF_SET_TYPES method
%       userdata    - a struct containing arbitrary data added by the user
%
%   Private Methods
%       add_named_set - Adds named subset of a particular type to the
%           object.
%
%       def_set_types - (must be implemented in the subclass) Returns a
%           struct defining the various set types, where the key is the
%           set type name, which must also be declared as a property in
%           the object's class, and the value is a string name used for
%           display purposes or a mp.set_manager object.
%
%           E.g.
%               function obj = def_set_types(obj)
%                   obj.set_types = struct(...
%                           'var', 'variable', ...
%                           'lin', 'linear constraint' ...
%                       );
%               end
%
%       init_set_types - Initializes the structures needed to track the
%           ordering and indexing of each set type and can be overridden
%           to initialize any additional data to be stored with each block
%           of each set type. Ideally, this would be called at the end
%           of the MP_IDX_MANAGER constructor, but this is not possible
%           in Octave 5.2 and earlier due to a bug related to altering
%           fields of an object not yet fully constructed.
%           Fixed in Octave v6.x: https://savannah.gnu.org/bugs/?52614
%           The workaround is to be sure that this method is called by
%           some subclass method(s) after object construction, but before
%           any other use (e.g. in display() and add_*() methods).
%
%       valid_named_set_type - Returns a label for the given named set
%           type if valid, empty otherwise.
%
%   Public Methods
%       add_<SET-TYPE> - (must be implemented in the subclass) Used to
%           add a block of the given set type to the object, using the
%           private ADD_NAMED_SET method internally followed by the
%           handling of any data specific to that set type.
%
%       copy - Makes a shallow copy of the object.
%
%       describe_idx - Identifies indices of a given set type.
%           E.g. variable 361 corresponds to Pg(68)
%
%       display - (must be implemented in the subclass) Displays the
%           object (called automatically when you omit the semicolon at
%           the command-line).
%
%       display_set - Prints to screen the indexing details for the
%           specified set type. Intended to be called by DISPLAY method.
%
%       get_idx - Returns index structure(s) for specified set type(s),
%           with starting/ending indices and number of elements for
%           each named (and optionally indexed) block.
%
%       get_userdata - Retreives values of user data stored in the object.
%
%       get - Return the value of any individual field.
%
%       getN - Returns the number of elements of any given set type,
%           optionally for a single named block.
%
%       init_indexed_name - Initializes the dimensions for an indexed
%           named set.
%
%       params_<SET-TYPE> - (must be implemented in the subclass)
%           Returns set-type-specific data for a given type.
%
%   The following is the structure of the data in the object, using a set
%   type named 'var' for illustration. Each field of .idx or .data is a
%   struct whose field names are the names of the corresponding blocks of
%   elements of that type (e.g. variables, constraints, etc.). They are
%   found in order in the corresponding .order field. The description next
%   to these fields gives the meaning of the value for each named sub-field.
%   E.g. obj.var.data.v0.Pg contains a vector of initial values for the 'Pg'
%   block of variables.
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

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%    es = struct();

    properties
        userdata = struct();
        set_types
    end     %% properties

    methods
        function obj = mp_idx_manager(s)
            % Constructor.
            % ::
            %
            %   obj = mp_idx_manager()
            %   obj = mp_idx_manager(a_struct)
            %   obj = mp_idx_manager(an_obj)

            if nargin > 0
                if isa(s, 'mp_idx_manager')
                    %% this copy constructor will not be inheritable under
                    %% Octave until the fix has been included for:
                    %%      https://savannah.gnu.org/bugs/?52614
                    if have_feature('octave')
                        s1 = warning('query', 'Octave:classdef-to-struct');
                        warning('off', 'Octave:classdef-to-struct');
                    end
                    props = fieldnames(s);
                    if have_feature('octave')
                        warning(s1.state, 'Octave:classdef-to-struct');
                    end
                    [~, k] = ismember('set_types', props);
                    props(k) = [];  %% remove 'set_types'
                    for k = 1:length(props)
                        obj = copy_prop(s, obj, props{k});
                    end
                elseif isstruct(s)
                    props = fieldnames(obj);
                    for k = 1:length(props)
                        if isfield(s, props{k})
                            obj = copy_prop(s, obj, props{k});
                        end
                    end
                else
                    error('mp_idx_manager.mp_idx_manager: input must be an ''mp_idx_manager'' object or a struct');
                end
            end

            obj.def_set_types();

            %% The INIT_SET_TYPES() method should ideally be (1) called here
            %% in the MP_IDX_MANAGER constructor and (2) skipped if constructed
            %% from an existing object, neither of which work in Octave 5.2
            %% and earlier due to inheritance issues in constructors:
            %%   https://savannah.gnu.org/bugs/?52614
            %%
            %% WORKAROUND: In some subclass method(s), call INIT_SET_TYPES()
            %%             automatically as needed after construction, but
            %%             before use, checking first to ensure that it
            %%             hasn't already been called.
%             if isempty(obj.????)        %% skip if already initialized (e.g.
%                 obj.init_set_types();   %% constructed from existing object)
%             end
        end

        function obj = init_set_types(obj)
            % Initialize indexing structures for each set type.

            %% Can allow def_set_types() to return a struct whose values
            %% are either char arrays or mp.set_manager objects.

            %% initialize each (set_type) field with base data structure
            for f = fieldnames(obj.set_types)'
                nis = obj.set_types.(f{1});
                if ischar(nis)
                    nis = mp.set_manager(nis);
                else
                    nis = nis.copy();
                end
                obj.(f{1}) = nis;
            end
        end

        function new_obj = copy(obj)
            % Duplicate the object.

            %% initialize copy
            new_obj = eval(class(obj));  %% create new object
            new_obj.init_set_types();

            %% copy properties/fields
            if have_feature('octave')
                s1 = warning('query', 'Octave:classdef-to-struct');
                warning('off', 'Octave:classdef-to-struct');
            end
            props = fieldnames(obj);
            if have_feature('octave')
                warning(s1.state, 'Octave:classdef-to-struct');
            end
            [~, k] = ismember('set_types', props);
            props(k) = [];  %% remove 'set_types'
            for k = 1:length(props)
                new_obj = copy_prop(obj, new_obj, props{k});
            end
        end

        function display_set(obj, stype, sname)
            % Display indexing information for a given set type.

            obj.(stype).display(stype);
        end

        obj = add_named_set(obj, set_type, name, idx, N, varargin)

        label = describe_idx(obj, set_type, idxs)

        varargout = get_idx(obj, varargin)

        rv = get_userdata(obj, name)

        val = get(obj, varargin)

        N = getN(obj, set_type, name, idx)

        obj = init_indexed_name(obj, set_type, name, dim_list)

        s = set_type_idx_map(obj, set_type, idxs, group_by_name)

        str = valid_named_set_type(obj, set_type)
    end     %% methods
end         %% classdef

function d = copy_prop(s, d, prop)
    if isa(s.(prop), 'mp.sm_quad_cost')
        d.(prop) = s.(prop).copy('mp.sm_quad_cost_legacy');
    elseif isa(s.(prop), 'mp.set_manager')
        d.(prop) = s.(prop).copy();
    elseif isa(d.(prop), 'mp.set_manager')
        d.(prop) = nested_struct_copy( ...
            d.(prop), s.(prop));
    else
        d.(prop) = s.(prop);
    end
end
