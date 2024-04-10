classdef (Abstract) dm_element < handle
% mp.dm_element - Abstract base class for |MATPOWER| **data model element** objects.
%
% A data model element object encapsulates all of the input and output data
% for a particular element type. All data model element classes inherit from
% mp.dm_element and each element type typically implements its own subclass.
% A given data model element object contains the data for all instances of
% that element type, stored in one or more table data structures.
%
% Defines the following columns in the main data table, which are inherited
% by all subclasses:
%
%   ==============  ============  ==========================================
%   Name            Type          Description
%   ==============  ============  ==========================================
%   ``uid``         *integer*     unique ID
%   ``name``        *char array*  element name
%   ``status``      *boolean*     true = online, false = offline
%   ``source_uid``  *undefined*   intended for any info required to link
%                                 back to element instance in source data
%   ==============  ============  ==========================================
%
% By convention, data model element variables are named ``dme`` and data model
% element class names begin with ``mp.dme``.
%
% In addition to being containers for the data itself, data model elements are
% responsible for handling the on/off status of each element, preparation of
% parameters needed by network and mathematical models, definition of
% connections with other elements, defining solution data to be updated when
% exporting, and pretty-printing of data to the console or file.
%
% Elements that create nodes (e.g. buses) are called **junction** elements.
% Elements that define ports (e.g. generators, branches, loads) can connect
% the ports of a particular instance to the nodes of a particular instance of
% a junction element by specifying two pieces of information for each port:
%
%   - the **type** of junction element it connects to
%   - the **index** of the specific junction element
%
% mp.dm_element Properties:
%   * tab - main data table
%   * nr - total number of rows in table
%   * n - number of online elements
%   * ID2i - max(ID) x 1 vector, maps IDs to row indices
%   * on - ``n`` x 1 vector of row indices of online elements
%   * off - (``nr``-``n``) x 1 vector of row indices of offline elements
%   * i2on - ``nr`` x 1 vector mapping row index to index in :attr:`on`/:attr:`off` respectively
%
% mp.dm_element Methods:
%   * name - get name of element type, e.g. ``'bus'``, ``'gen'``
%   * label - get singular label for element type, e.g. ``'Bus'``, ``'Generator'``
%   * labels - get plural label for element type, e.g. ``'Buses'``, ``'Generators'``
%   * cxn_type - type(s) of junction element(s) to which this element connects
%   * cxn_idx_prop - name(s) of property(ies) containing indices of junction elements
%   * cxn_type_prop - name(s) of property(ies) containing types of junction elements
%   * table_exists - check for existence of data in main data table
%   * main_table_var_names - names of variables (columns) in main data table
%   * export_vars - names of variables to be exported by DMCE to data source
%   * export_vars_offline_val - values of export variables for offline elements
%   * dm_converter_element - get corresponding data model converter element
%   * copy - create a duplicate of the data model element object
%   * count - determine number of instances of this element in the data
%   * initialize - initialize (online/offline) status of each element
%   * ID - return unique ID's for all or indexed rows
%   * init_status - initialize ``status`` column
%   * update_status - update (online/offline) status based on connectivity, etc
%   * build_params - extract/convert/calculate parameters for online elements
%   * rebuild - rebuild object, calling count(), initialize(), build_params()
%   * display - display the data model element object
%   * pretty_print - pretty-print data model element to console or file
%   * pp_have_section - true if pretty-printing for element has specified section
%   * pp_rows - indices of rows to include in pretty-printed output
%   * pp_get_headers - get pretty-printed headers for this element/section
%   * pp_get_footers - get pretty-printed footers for this element/section
%   * pp_data - pretty-print the data for this element/section
%   * pp_have_section_cnt - true if pretty-printing for element has **counts** section
%   * pp_data_cnt - pretty-print the **counts** data for this element
%   * pp_have_section_sum - true if pretty-printing for element has **summary** section
%   * pp_data_sum - pretty-print the **summary** data for this element
%   * pp_have_section_ext - true if pretty-printing for element has **extremes** section
%   * pp_data_ext - pretty-print the **extremes** data for this element
%   * pp_have_section_det - true if pretty-printing for element has **details** section
%   * pp_get_title_det - get title of **details** section for this element
%   * pp_get_headers_det - get pretty-printed **details** headers for this element
%   * pp_get_footers_det - get pretty-printed **details** footers for this element
%   * pp_data_det - pretty-print the **details** data for this element
%   * pp_data_row_det - get pretty-printed row of **details** data for this element
%
% See the :ref:`sec_dm_element` section in the |MATPOWER-Dev-Manual| for more
% information.
%
% See also mp.data_model.

%   MATPOWER
%   Copyright (c) 2020-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        tab     % *(table)* main data table
        nr      % *(integer)* total number of rows in table
        n       % *(integer)* number of online elements
        ID2i    % *(integer)* max(ID) x 1 vector, maps IDs to row indices
        on      % *(integer)* ``n`` x 1 vector of row indices of online elements
        off     % *(integer)* (``nr``-``n``) x 1 vector of row indices of offline elements
        i2on    % *(integer)* ``nr`` x 1 vector mapping row index to index in :attr:`on`/:attr:`off` respectively
    end     %% properties

    methods
        function name = name(obj)
            % Get name of element type, e.g. ``'bus'``, ``'gen'``.
            % ::
            %
            %   name = dme.name()
            %
            % Output:
            %   name (char array) : name of element type, must be a valid
            %       struct field name
            %
            % Implementation provided by an element type specific subclass.

            name = '';      %% e.g. 'bus', 'gen'
        end

        function label = label(obj)
            % Get singular label for element type, e.g. ``'Bus'``, ``'Generator'``.
            % ::
            %
            %   label = dme.label()
            %
            % Output:
            %   label (char array) : user-visible label for element type,
            %       when singular
            %
            % Implementation provided by an element type specific subclass.

            label = '';     %% e.g. 'Bus', 'Generator'
        end

        function label = labels(obj)
            % Get plural label for element type, e.g. ``'Buses'``, ``'Generators'``.
            % ::
            %
            %   label = dme.labels()
            %
            % Output:
            %   label (char array) : user-visible label for element type,
            %       when plural
            %
            % Implementation provided by an element type specific subclass.

            label = '';     %% e.g. 'Buses', 'Generators'
        end

        function name = cxn_type(obj)
            % Type(s) of junction element(s) to which this element connects.
            % ::
            %
            %   name = dme.cxn_type()
            %
            % Output:
            %   name (char array or cell array of char arrays) : name(s) of
            %       type(s) of junction elements, i.e. node-creating elements
            %       (e.g. ``'bus'``), to which this element connects
            %
            % Assuming an element with *nc* connections, there are three
            % options for the return value:
            %
            %   1. Single char array with one type that applies to all
            %      connections, cxn_idx_prop() returns *empty*.
            %   2. Cell array with *nc* elements, one for each connection,
            %      cxn_idx_prop() returns *empty*.
            %   3. Cell array of valid junction element types, cxn_idx_prop()
            %      return value *not empty*.
            %
            % See the :ref:`sec_dm_element_cxn` section in the
            % |MATPOWER-Dev-Manual| for more information.
            %
            % Implementation provided by an element type specific subclass.
            %
            % See also cxn_idx_prop, cxn_type_prop.

            name = '';
        end

        function name = cxn_idx_prop(obj)
            % Name(s) of property(ies) containing indices of junction elements.
            % ::
            %
            %   name = dme.cxn_idx_prop()
            %
            % Output:
            %   name (char array or cell array of char arrays) : name(s) of
            %       property(ies) containing indices of junction elements that
            %       define connections (e.g. ``{'fbus', 'tbus'}``)
            %
            % See the :ref:`sec_dm_element_cxn` section in the
            % |MATPOWER-Dev-Manual| for more information.
            %
            % Implementation provided by an element type specific subclass.
            %
            % See also cxn_type, cxn_type_prop.

            name = '';
        end

        function name = cxn_type_prop(obj)
            % Name(s) of property(ies) containing types of junction elements.
            % ::
            %
            %   name = dme.cxn_type_prop()
            %
            % Output:
            %   name (char array or cell array of char arrays) : name(s) of
            %       properties containing type of junction elements for each
            %       connection
            %
            %       *Note:* If not empty, dimension must match cxn_idx_prop()
            %
            % This is only used if the junction element type can vary by
            % individual element, e.g. some elements of this type connect to
            % one kind of bus, some to another kind. Otherwise, it returns
            % an empty string and the junction element types for the
            % connections are determined solely by cxn_type().
            %
            % See the :ref:`sec_dm_element_cxn` section in the
            % |MATPOWER-Dev-Manual| for more information.
            %
            % Implementation provided by an element type specific subclass.
            %
            % See also cxn_type, cxn_idx_prop.

            name = '';
        end

        function TorF = table_exists(obj)
            % Check for existence of data in main data table.
            % ::
            %
            %   TorF = dme.table_exists()
            %
            % Output:
            %   TorF (boolean) : true if main data table is not empty

            TorF = ~isempty(obj.tab);
        end

        function names = main_table_var_names(obj)
            % Names of variables (columns) in main data table.
            % ::
            %
            %   names = dme.main_table_var_names()
            %
            % Output:
            %   names (cell array of char arrays) : names of variables
            %       (columns) in main table
            %
            % This base class includes the following variables
            % ``{'uid', 'name', 'status', 'source_uid'}`` which are common to
            % all element types and should therefore be included in all
            % subclasses. That is, subclass methods should append their
            % additional fields to those returned by this parent method. For
            % example, a subclass method would like something like the
            % following::
            %
            %   function names = main_table_var_names(obj)
            %       names = horzcat( main_table_var_names@mp.dm_element(obj), ...
            %       {'subclass_var1', 'subclass_var2'} );
            %   end

            names = {'uid', 'name', 'status', 'source_uid'};
        end

        function vars = export_vars(obj)
            % Names of variables to be exported by DMCE to data source.
            % ::
            %
            %   vars = dme.export_vars()
            %
            % Output:
            %   vars (cell array of char arrays) : names of variables to export
            %
            % Return the names of the variables the data model converter
            % element needs to export to the data source. This is typically
            % the list of variables updated by the solution process, e.g. bus
            % voltages, line flows, etc.

            vars = {};
        end

        function s = export_vars_offline_val(obj)
            % Values of export variables for offline elements.
            % ::
            %
            %   s = dme.export_vars_offline_val()
            %
            % Output:
            %   s (struct) : keys are export variable names, values are
            %       the corresponding values to assign to these variables
            %       for offline elements.
            %
            % Returns a struct defining the values of export variables for
            % offline elements. Called by mp.mm_element.data_model_update
            % to define how to set export variables for offline elements.
            %
            % Export variables not found in the struct are not modified.
            %
            % For example, ``s = struct('va', 0, 'vm', 1)`` would assign
            % the value 0 to the ``va`` variable and 1 to the ``vm`` variable
            % for any offline elements.
            %
            % See also export_vars.

            s = struct();
        end

        function dmce = dm_converter_element(obj, dmc, name)
            % Get corresponding data model converter element.
            % ::
            %
            %   dmce = dme.dm_converter_element(dmc)
            %   dmce = dme.dm_converter_element(dmc, name)
            %
            % Inputs:
            %   dmc (mp.dm_converter) : data model converter object
            %   name (char array) : *(optional)* name of element type
            %       *(default is name of this object)*
            %
            % Output:
            %   dmce (mp.dmc_element) : data model converter element object

            if nargin < 3
                name = obj.name;
            end
            dmce = dmc.elements.(name);
        end

        function new_obj = copy(obj)
            % Create a duplicate of the data model element object.
            % ::
            %
            %   new_dme = dme.copy()
            %
            % Output:
            %   new_dme (mp.dm_element) : copy of data model element object

            %% make shallow copy of object
            new_obj = eval(class(obj));  %% create new object
            if have_feature('octave')
                s1 = warning('query', 'Octave:classdef-to-struct');
                warning('off', 'Octave:classdef-to-struct');
            end
            props = fieldnames(obj);
            if have_feature('octave')
                warning(s1.state, 'Octave:classdef-to-struct');
            end
            for k = 1:length(props)
                new_obj.(props{k}) = obj.(props{k});
            end
        end

        function nr = count(obj, dm)
            % Determine number of instances of this element in the data.
            % 
            % Store the count in the :attr:`nr` property.
            % ::
            %
            %   nr = dme.count(dm);
            %
            % Input:
            %   dm (mp.data_model) : data model
            %
            % Output:
            %   nr (integer) : number of instances (rows of data)
            %
            % Called for each element by the
            % :meth:`count() <mp.data_model.count>` method of mp.data_model
            % during the **count** stage of a data model build.
            %
            % See the :ref:`sec_building_data_model` section in the
            % |MATPOWER-Dev-Manual| for more information.

            if isempty(obj.tab)
                nr = 0;
            else
                nr = size(obj.tab, 1);
            end
            obj.nr = nr;
        end

        function obj = initialize(obj, dm)
            % Initialize a newly created data model element object.
            % ::
            %
            %   dme.initialize(dm)
            %
            % Input:
            %   dm (mp.data_model) : data model
            %
            % Initialize the (online/offline) status of each element and
            % create a mapping of ID to row index in the :attr:`ID2i` element
            % property, then call init_status().
            %
            % Called for each element by the
            % :meth:`initialize() <mp.data_model.initialize>` method of
            % mp.data_model during the **initialize** stage of a data model
            % build.
            %
            % See the :ref:`sec_building_data_model` section in the
            % |MATPOWER-Dev-Manual| for more information.

            %% create ID --> index mapping
            nr = obj.nr;
            ID = obj.ID;
            maxID = max(ID);
            if ~isempty(ID)
                if nr > 5000 && maxID > 10 * nr %% use sparse map (save memory)
                    ID2i = sparse(ID, ones(nr,1), 1:nr, maxID, 1);
                else                            %% use dense map (faster)
                    ID2i = accumarray([ID ones(nr,1)], 1:nr, [maxID 1]);
                end
            else
                ID2i = [];
            end
            obj.ID2i = ID2i;

            %% initial on/off status
            obj.init_status(dm);
        end

        function uid = ID(obj, idx)
            % Return unique ID's for all or indexed rows.
            % ::
            %
            %   uid = dme.ID()
            %   uid = dme.ID(idx)
            %
            % Input:
            %   idx (integer) : *(optional)* row index vector
            %
            % Return an *nr* x 1 vector of unique IDs for all rows, i.e. a map
            % of row index to unique ID or, if a row index vector is provided
            % just the ID's of the indexed rows.

            if nargin > 1 && ~isempty(idx)
                uid = obj.tab.uid(idx);
            else
                uid = obj.tab.uid;
            end
        end

        function obj = init_status(obj, dm)
            % Initialize ``status`` column.
            % ::
            %
            %   dme.init_status(dm)
            %
            % Input:
            %   dm (mp.data_model) : data model
            %
            % Called by initialize(). Does nothing in the base class.
        end

        function obj = update_status(obj, dm)
            % Update (online/offline) status based on connectivity, etc.
            % ::
            %
            %   dme.update_status(dm)
            %
            % Input:
            %   dm (mp.data_model) : data model
            %
            % Update status of each element based on connectivity or other
            % criteria and define element properties containing number and
            % row indices of online elements (:attr:`n` and :attr:`on`),
            % indices of offline elements (:attr:`off`), and mapping
            % (:attr:`i2on`) of row indices to corresponding entries in
            % :attr:`on` or :attr:`off`.
            %
            % Called for each element by the
            % :meth:`update_status() <mp.data_model.update_status>` method of
            % mp.data_model during the **update status** stage of a data model
            % build.
            %
            % See the :ref:`sec_building_data_model` section in the
            % |MATPOWER-Dev-Manual| for more information.

            status = obj.tab.status;
            if ~isempty(status)
                obj.on  = find(  status );
                obj.off = find( ~status );
                obj.n   = length(obj.on);
                obj.i2on = zeros(obj.nr, 1);
                obj.i2on(obj.on ) = (1:obj.n);
                obj.i2on(obj.off) = (1:length(obj.off));
            end
        end

        function obj = build_params(obj, dm)
            % Extract/convert/calculate parameters for online elements.
            % ::
            %
            %   dme.build_params(dm)
            %
            % Input:
            %   dm (mp.data_model) : data model
            %
            % Extract/convert/calculate parameters as necessary for online
            % elements from the original data tables (e.g. p.u. conversion,
            % initial state, etc.) and store them in element-specific
            % properties.
            %
            % Called for each element by the
            % :meth:`build_params() <mp.data_model.build_params>` method of
            % mp.data_model during the **build parameters** stage of a data
            % model build.
            %
            % See the :ref:`sec_building_data_model` section in the
            % |MATPOWER-Dev-Manual| for more information.
            %
            % Does nothing in the base class.
        end

        function obj = rebuild(obj, dm)
            % Rebuild object, calling count(), initialize(), build_params().
            % ::
            %
            %   dme.rebuild(dm)
            %
            % Input:
            %   dm (mp.data_model) : data model
            %
            % Typically used after modifying data in the main table.

            obj.count(dm);
            obj.initialize(dm);
            obj.build_params(dm);
        end

        function display(obj)
            % Display the data model element object.
            %
            % This method is called automatically when omitting a semicolon
            % on a line that retuns an object of this class.
            %
            % Displays the details of the elements, including total number
            % of rows, number of online elements, and the main data table.

%             if have_feature('octave')
%                 struct(obj)
%             else
%                 display@handle(obj)
%             end
            fprintf('DATA MODEL ELEMENT NAME  : %s\n', obj.name);
            fprintf('DATA MODEL ELEMENT CLASS : %s\n', class(obj));
            fprintf('    # OF ROWS            : %d\n', obj.nr);
            fprintf('    # OF ONLINE ELEMENTS : %d\n', obj.n);
            disp(obj.tab);
        end

        function obj = pretty_print(obj, dm, section, out_e, mpopt, fd, pp_args)
            % Pretty print data model element to console or file.
            % ::
            %
            %   dme.pretty_print(dm, section, out_e, mpopt, fd, pp_args)
            %
            % Inputs:
            %   dm (mp.data_model) : data model
            %   section (char array) : section identifier, e.g. ``'cnt'``,
            %       ``'sum'``, ``'ext'``, or ``'det'``, for **counts**,
            %       **summary**, **extremes**, or **details** sections,
            %       respectively
            %   out_e (boolean) : output control flag for this element/section
            %   mpopt (struct) : |MATPOWER| options struct
            %   fd (integer) : *(optional, default = 1)* file identifier to
            %       use for printing, (1 for standard output, 2 for standard
            %       error)
            %   pp_args (struct) : arbitrary struct of additional pretty
            %       printing arguments passed to all sub-methods, allowing
            %       a single sub-method to be used for multiple output
            %       portions (e.g. for active and reactive power) by passing
            %       in a different argument; by convention, arguments for a
            %       ``branch`` element, for example, are passed in
            %       ``pp_args.branch``, etc.

            if out_e && obj.pp_have_section(section, mpopt, pp_args);
                %% get indices of relevant rows
                rows = obj.pp_rows(dm, section, out_e, mpopt, pp_args);

                %% if there are rows to show
                if ~isempty(rows) && rows(1) ~= 0
                    %% title & headers
                    h = obj.pp_get_headers(dm, section, out_e, mpopt, pp_args);
                    for k = 1:length(h)
                        fprintf(fd, '%s\n', h{k});
                    end

                    %% data
                    obj.pp_data(dm, section, rows, out_e, mpopt, fd, pp_args);

                    %% footers
                    f = obj.pp_get_footers(dm, section, out_e, mpopt, pp_args);
                    for k = 1:length(f)
                        fprintf(fd, '%s\n', f{k});
                    end
                end
            end
        end

        function TorF = pp_have_section(obj, section, mpopt, pp_args)
            % True if pretty-printing for element has specified section.
            % ::
            %
            %   TorF = dme.pp_have_section(section, mpopt, pp_args)
            %
            % Inputs:
            %    : see pretty_print() for details
            %
            % Output:
            %   TorF (boolean) : true if output includes specified section
            %
            % Implementation handled by section-specific *pp_have_section*
            % methods or :meth:`pp_have_section_other() <mp.dme_shared_opf.pp_have_section_other>`.
            %
            % See also pp_have_section_cnt, pp_have_section_sum,
            % pp_have_section_ext, pp_have_section_det.

            switch section
                case 'cnt'
                    TorF = obj.pp_have_section_cnt(mpopt, pp_args);
                case 'sum'
                    TorF = obj.pp_have_section_sum(mpopt, pp_args);
                case 'ext'
                    TorF = obj.pp_have_section_ext(mpopt, pp_args);
                case 'det'
                    TorF = obj.pp_have_section_det(mpopt, pp_args);
                otherwise
                    TorF = obj.pp_have_section_other(section, mpopt, pp_args);
            end
        end

        function rows = pp_rows(obj, dm, section, out_e, mpopt, pp_args)
            % Indices of rows to include in pretty-printed output.
            % ::
            %
            %   rows = dme.pp_rows(dm, section, out_e, mpopt, pp_args)
            %
            % Inputs:
            %    : see pretty_print() for details
            %
            % Output:
            %   rows (integer) : index vector of rows to be included in output
            %
            %       - 0 = no rows
            %       - -1 = all rows
            %
            % Includes all rows by default.

            switch section
                case {'cnt', 'sum', 'ext', 'det'}
                    rows = -1;  %% all rows
                otherwise
                    rows = obj.pp_rows_other(dm, section, out_e, mpopt, pp_args);
            end
        end

        function h = pp_get_headers(obj, dm, section, out_e, mpopt, pp_args)
            % Get pretty-printed headers for this element/section.
            % ::
            %
            %   h = dme.pp_get_headers(dm, section, out_e, mpopt, pp_args)
            %
            % Inputs:
            %    : see pretty_print() for details
            %
            % Output:
            %   h (cell array of char arrays) : lines of pretty printed header
            %       output for this element/section
            %
            % Empty by default for counts, summary and extremes sections,
            % and handled by pp_get_headers_det() for details section.

            switch section
                case {'cnt', 'sum', 'ext'}
                    h = {};
                case 'det'
                    h = obj.pp_get_headers_det(dm, out_e, mpopt, pp_args);
                otherwise
                    h = obj.pp_get_headers_other(dm, section, out_e, mpopt, pp_args);
            end
        end

        function f = pp_get_footers(obj, dm, section, out_e, mpopt, pp_args)
            % Get pretty-printed footers for this element/section.
            % ::
            %
            %   f = dme.pp_get_footers(dm, section, out_e, mpopt, pp_args)
            %
            % Inputs:
            %    : see pretty_print() for details
            %
            % Output:
            %   f (cell array of char arrays) : lines of pretty printed footer
            %       output for this element/section
            %
            % Empty by default for counts, summary and extremes sections,
            % and handled by pp_get_headers_det() for details section.

            switch section
                case {'cnt', 'sum', 'ext'}
                    f = {};
                case 'det'
                    f = obj.pp_get_footers_det(dm, out_e, mpopt, pp_args);
                otherwise
                    f = obj.pp_get_footers_other(dm, section, out_e, mpopt, pp_args);
            end
        end

        function obj = pp_data(obj, dm, section, rows, out_e, mpopt, fd, pp_args)
            % Pretty-print the data for this element/section.
            % ::
            %
            %   dme.pp_data(dm, section, rows, out_e, mpopt, fd, pp_args)
            %
            % Inputs:
            %   rows (integer) : indices of rows to include, from pp_rows()
            %   ... : see pretty_print() for details of other inputs
            %
            % Implementation handled by section-specific *pp_data*
            % methods or :meth:`pp_data_other() <mp.dme_shared_opf.pp_data_other>`.
            %
            % See also pp_data_cnt, pp_data_sum, pp_data_ext, pp_data_det.

            switch section
                case 'cnt'
                    obj.pp_data_cnt(dm, rows, out_e, mpopt, fd, pp_args);
                case 'sum'
                    obj.pp_data_sum(dm, rows, out_e, mpopt, fd, pp_args);
                case 'ext'
                    obj.pp_data_ext(dm, rows, out_e, mpopt, fd, pp_args);
                case 'det'
                    obj.pp_data_det(dm, rows, out_e, mpopt, fd, pp_args);
                otherwise
                    obj.pp_data_other(dm, section, rows, out_e, mpopt, fd, pp_args);
            end
        end

        function TorF = pp_have_section_cnt(obj, mpopt, pp_args)
            % True if pretty-printing for element has **counts** section.
            % ::
            %
            %   TorF = dme.pp_have_section_cnt(mpopt, pp_args)
            %
            % Default is **true**.
            %
            % See also pp_have_section.

            TorF = true;    %% all elements have count section by default
        end

        function obj = pp_data_cnt(obj, dm, rows, out_e, mpopt, fd, pp_args)
            % Pretty-print the **counts** data for this element.
            % ::
            %
            %   dme.pp_data_cnt(dm, rows, out_e, mpopt, fd, pp_args)
            %
            % See also pp_data.

            if obj.n
                on = sprintf('%d', obj.n);
            else
                on = '-';
            end
            if obj.nr - obj.n
                off = sprintf('%d', obj.nr - obj.n);
            else
                off = '-';
            end
            fprintf(fd, '  %-20s%7s %7s %7d\n', obj.labels, on, off, obj.nr);
        end

        function TorF = pp_have_section_sum(obj, mpopt, pp_args)
            % True if pretty-printing for element has **summary** section.
            % ::
            %
            %   TorF = dme.pp_have_section_sum(mpopt, pp_args)
            %
            % Default is **false**.
            %
            % See also pp_have_section.

            TorF = false;   %% no summary section for elements by default
        end

        function obj = pp_data_sum(obj, dm, rows, out_e, mpopt, fd, pp_args)
            % Pretty-print the **summary** data for this element.
            % ::
            %
            %   dme.pp_data_sum(dm, rows, out_e, mpopt, fd, pp_args)
            %
            % Does nothing by default.
            %
            % See also pp_data.
        end

        function TorF = pp_have_section_ext(obj, mpopt, pp_args)
            % True if pretty-printing for element has **extremes** section.
            % ::
            %
            %   TorF = dme.pp_have_section_ext(mpopt, pp_args)
            %
            % Default is **false**.
            %
            % See also pp_have_section.

            TorF = false;   %% no ext section for elements by default
        end

        function obj = pp_data_ext(obj, dm, rows, out_e, mpopt, fd, pp_args)
            % Pretty-print the **extremes** data for this element.
            % ::
            %
            %   dme.pp_data_ext(dm, rows, out_e, mpopt, fd, pp_args)
            %
            % Does nothing by default.
            %
            % See also pp_data.
        end

        function TorF = pp_have_section_det(obj, mpopt, pp_args)
            % True if pretty-printing for element has **details** section.
            % ::
            %
            %   TorF = dme.pp_have_section_det(mpopt, pp_args)
            %
            % Default is **false**.
            %
            % See also pp_have_section.

            TorF = false;   %% no det section for elements by default
        end

        function str = pp_get_title_det(obj, mpopt, pp_args)
            % Get title of **details** section for this element.
            % ::
            %
            %   str = dme.pp_get_title_det(mpopt, pp_args)
            %
            % Inputs:
            %    : see pretty_print() for details
            %
            % Output:
            %   str (char array) : title of details section, e.g.
            %       ``'Bus Data'``, ``'Generator Data'``, etc.
            %
            % Called by pp_get_headers_det() to insert title into detail
            % section header.

            if obj.pp_have_section_det(mpopt, pp_args)
                str = sprintf('%s Data', obj.label);
            else
                str = '';
            end
        end

        function h = pp_get_headers_det(obj, dm, out_e, mpopt, pp_args)
            % Get pretty-printed **details** headers for this element.
            % ::
            %
            %   h = dme.pp_get_headers_det(dm, out_e, mpopt, pp_args)
            %
            % See also pp_get_headers.

            str = obj.pp_get_title_det(mpopt, pp_args);
            if isempty(str)
                h = {};
            else
                h = dm.pp_section_label(str);
            end
        end

        function f = pp_get_footers_det(obj, dm, out_e, mpopt, pp_args)
            % Get pretty-printed **details** footers for this element.
            % ::
            %
            %   f = dme.pp_get_footers_det(dm, out_e, mpopt, pp_args)
            %
            % Empty by default.
            %
            % See also pp_get_footers.

            f = {};
        end

        function obj = pp_data_det(obj, dm, rows, out_e, mpopt, fd, pp_args)
            % Pretty-print the **details** data for this element.
            % ::
            %
            %   dme.pp_data_det(dm, rows, out_e, mpopt, fd, pp_args)
            %
            % Calls pp_data_row_det() for each row.
            %
            % See also pp_data, pp_data_row_det.

            if obj.pp_have_section_det(mpopt, pp_args)
                assert(rows == -1, 'dm_elemnt/pp_data_det: ''rows'' is expected to be -1, indicating all rows');
                for k = 1:obj.nr
                    fprintf(fd, '%s\n', ...
                        obj.pp_data_row_det(dm, k, out_e, mpopt, fd, pp_args));
                end
            end
        end

        function str = pp_data_row_det(obj, dm, k, out_e, mpopt, fd, pp_args)
            % Get pretty-printed row of **details** data for this element.
            % ::
            %
            %   str = dme.pp_data_row_det(dm, k, out_e, mpopt, fd, pp_args)
            %
            % Inputs:
            %   k (integer) : index of row to print
            %   ... : see pretty_print() for details of other inputs
            %
            % Output:
            %   str (char array) : row of data *(without newline)*
            %
            % Called by pp_data_det() for each row.

            str = '';
        end
    end     %% methods
end         %% classdef
