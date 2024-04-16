classdef (Abstract) dmc_element < handle
% mp.dmc_element- Abstract base class for **data model converter element** objects.
%
% A data model converter element object implements the functionality needed
% to import and export a particular element type from and to a given
% data format. All data model converter element classes inherit from
% mp.dmc_element and each element type typically implements its own subclass.
%
% By convention, data model converter element variables are named ``dmce``
% and data model converter element class names begin with ``mp.dmce``.
%
% Typically, much of the import/export functionality for a particular concrete
% subclass can be defined simply by implementing the table_var_map() method.
%
% mp.dmc_element Methods:
%   * name - get name of element type, e.g. ``'bus'``, ``'gen'``
%   * data_model_element - get corresponding data model element
%   * data_field - get name of field in data source corresponding to default data table
%   * data_subs - get subscript reference struct for accessing data source
%   * data_exists - check if default field exists in data source
%   * get_import_spec - get import specification
%   * get_export_spec - get export specification
%   * get_import_size - get dimensions of data to be imported
%   * get_export_size - get dimensions of data to be exported
%   * table_var_map - get variable map for import/export
%   * import - import data from data source into data model element
%   * import_table_values - import table values for given import specification
%   * get_input_table_values - get values to insert in data model element table
%   * import_col - extract and optionally modify values from data source column
%   * export - export data from data model element to data source
%   * export_table_values - export table values for given import specification
%   * init_export_data - initialize data source for export from data model element
%   * default_export_data_table - create default (empty) data table for data source
%   * default_export_data_nrows - get number of rows default_export_data_table()
%   * export_col - export a variable (table column) to the data source
%
% See the :ref:`sec_dmc_element` section in the |MATPOWER-Dev-Manual| for more
% information.
%
% See also mp.dm_converter.

%   MATPOWER
%   Copyright (c) 2021-2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties


    methods
        function name = name(obj)
            % Get name of element type, e.g. ``'bus'``, ``'gen'``.
            % ::
            %
            %   name = dmce.name()
            %
            % Output:
            %   name (char array) : name of element type, must be a valid
            %       struct field name
            %
            % Implementation provided by an element type specific subclass.

            name = '';      %% e.g. 'bus', 'gen'
        end

        function dme = data_model_element(obj, dm, name)
            % Get the corresponding data model element.
            % ::
            %
            %   dme = dmce.data_model_element(dm)
            %   dme = dmce.data_model_element(dm, name)
            %
            % Inputs:
            %   dm (mp.data_model) : data model object
            %   name (char array) : *(optional)* name of element type
            %       *(default is name of this object)*
            %
            % Output:
            %   dme (mp.dm_element) : data model element object

            if nargin < 3
                name = obj.name;
            end
            dme = dm.elements.(name);
        end

        function df = data_field(obj)
            % Get name of field in data source corresponding to default data table.
            % ::
            %
            %   df = dmce.data_field()
            %
            % Output:
            %   df (char array) : field name

            df = '';        %% name of field in d for default data table
        end

        function s = data_subs(obj)
            % Get subscript reference struct for accessing data source.
            % ::
            %
            %   s = dmce.data_subs()
            %
            % Output:
            %   s (struct) : same as the ``s`` input argument to the built-in
            %       :func:`subsref`, to access this element's data in data
            %       source, with fields:
            %
            %           - ``type`` -- character vector or string containing
            %             ``'()'``, ``'{}'``, or ``'.'`` specifying the
            %             subscript type
            %           - ``subs`` -- cell array, character vector, or string
            %             containing the actual subscripts
            %
            % The default implementation in this base class uses the
            % return value of the data_field() method to access a field of
            % the data source struct. That is::
            %
            %   s = struct('type', '.', 'subs', dmce.data_field());

            s = struct('type', '.', 'subs', obj.data_field());
        end

        function TorF = data_exists(obj, d)
            % Check if default field exists in data source.
            % ::
            %
            %   TorF = dmce.data_exists(d)
            %
            % Input:
            %   d : data source, type depends on the implementing subclass
            %       (e.g. |MATPOWER| case struct for mp.dm_converter_mpc2)
            %
            % Output:
            %   TorF (boolean) : true if field exists
            %
            % Check if value returned by data_field() exists as a field in ``d``.

            TorF = isfield(d, obj.data_field());
        end

        function spec = get_import_spec(obj, dme, d)
            % Get import specification.
            % ::
            %
            %   spec = dmce.get_import_spec(dme, d)
            %
            % Inputs:
            %   dme (mp.dm_element) : data model element object
            %   d : data source, type depends on the implementing subclass
            %       (e.g. |MATPOWER| case struct for mp.dm_converter_mpc2)
            %
            % Output:
            %   spec (struct) : import specification, with
            %       keys:
            %
            %           - ``'subs'`` - subscript reference struct for accessing
            %             data source, as returned by data_subs()
            %           - ``'nr'``, ``'nc'``, ``'r'`` - number of rows, number
            %             of columns, row index vector, as returned by
            %             get_import_size()
            %           - ``'vmap'`` - variable map, as returned by
            %             table_var_map()
            %
            % See also get_export_spec.


            spec = obj.get_spec(dme, d, 1);
        end

        function spec = get_export_spec(obj, dme, d)
            % Get export specification.
            % ::
            %
            %   spec = dmce.get_export_spec(dme, d)
            %
            % Inputs:
            %   dme (mp.dm_element) : data model element object
            %   d : data source, type depends on the implementing subclass
            %       (e.g. |MATPOWER| case struct for mp.dm_converter_mpc2)
            %
            % Output:
            %   spec (struct) : export specification, see get_import_spec()
            %
            % See also get_import_spec.

            spec = obj.get_spec(dme, d, 0);
        end

        function [nr, nc, r] = get_import_size(obj, d)
            % Get dimensions of data to be imported.
            % ::
            %
            %   [nr, nc, r] = dmce.get_import_size(d)
            %
            % Input:
            %   d : data source, type depends on the implementing subclass
            %       (e.g. |MATPOWER| case struct for mp.dm_converter_mpc2)
            %
            % Outputs:
            %   nr (integer) : number of rows of data
            %   nc (integer) : number of columns of data
            %   r (integer) : optional index vector *(empty by default)* of
            %       rows in data source field that correspond to data to be
            %       imported

            if obj.data_exists(d)
                %% use size of default table
                [nr, nc] = size(subsref(d, obj.data_subs()));
            else
                [nr, nc] = deal(0);
            end
            r = [];                         %% all rows
        end

        function [nr, nc, r] = get_export_size(obj, dme)
            % Get dimensions of data to be exported.
            % ::
            %
            %   [nr, nc, r] = dmce.get_export_size(dme)
            %
            % Input:
            %   dme (mp.dm_element) : data model element object
            %
            % Outputs:
            %   nr (integer) : number of rows of data
            %   nc (integer) : number of columns of data
            %   r (integer) : optional index vector *(empty by default)* of
            %       rows in main table of ``dme`` that correspond to data to
            %       be exported

            [nr, nc] = size(dme.tab);   %% use size of default table
            r = [];                     %% all rows
        end

        function vmap = table_var_map(obj, dme, d)
            % Get variable map for import/export.
            % ::
            %
            %   vmap = dmce.table_var_map(dme, d)
            %
            % Inputs:
            %   dme (mp.dm_element) : data model element object
            %   d : data source, type depends on the implementing subclass
            %       (e.g. |MATPOWER| case struct for mp.dm_converter_mpc2)
            %
            % Output:
            %   vmap (struct) : variable map, see :numref:`tab_var_map` in the
            %       |MATPOWER-Dev-Manual| for details
            %
            % This method initializes each entry to ``{'col', []}`` by default,
            % so subclasses only need to assign ``vmap.(vn){2}`` for columns
            % that map directly from a column of the data source.

            %% initialize with vmap.(<name>) = {'col', []}, for all <name>
            names = dme.main_table_var_names();
            vals = cell(size(names));
            [vals{:}] = deal({'col', []});
            vmap = cell2struct(vals, names, 2);
        end

        function dme = import(obj, dme, d, var_names, ridx)
            % Import data from data source into data model element.
            % ::
            %
            %   dme = dmce.import(dme, d, var_names, ridx)
            %
            % Inputs:
            %   dme (mp.dm_element) : data model element object
            %   d : data source, type depends on the implementing subclass
            %       (e.g. |MATPOWER| case struct for mp.dm_converter_mpc2)
            %   var_names (cell array) : *(optional)* list of names of
            %       variables (columns of main table) to import *(default is
            %       all variables)*
            %   ridx (integer) : *(optional)* vector of row indices of data
            %       to import *(default is all rows)*
            %
            % Output:
            %   dme (mp.dm_element) : updated data model element object
            %
            % See also export.

            if nargin < 5
                ridx = [];
                if nargin < 4
                    var_names = {};
                end
            end

            %% get conversion spec
            spec = obj.get_import_spec(dme, d);

            %% import main table
            dme = obj.import_table_values(dme, d, spec, var_names, ridx);
        end

        function dme = import_table_values(obj, dme, d, spec, var_names, ridx)
            % Import table values for given import specification.
            % ::
            %
            %   dme = dmce.import_table_values(dme, d, spec, var_names, ridx)
            %
            % Inputs:
            %   dme (mp.dm_element) : data model element object
            %   d : data source, type depends on the implementing subclass
            %       (e.g. |MATPOWER| case struct for mp.dm_converter_mpc2)
            %   spec (struct) : import specification, see get_import_spec()
            %   var_names (cell array) : *(optional)* list of names of
            %       variables (columns of main table) to import *(default is
            %       all variables)*
            %   ridx (integer) : *(optional)* vector of row indices of data
            %       to import *(default is all rows)*
            %
            % Output:
            %   dme (mp.dm_element) : updated data model element object
            %
            % Called by import().

            if nargin < 6
                ridx = [];
                if nargin < 5
                    var_names = {};
                end
            end

            if dme.table_exists()
                update = 1;     %% update values in existing dme table
            else
                update = 0;     %% use values to create dme table
                var_names = dme.main_table_var_names();
            end

            %% get values from input data structure
            vals = obj.get_input_table_values(d, spec, var_names, ridx);

            if ~isempty(vals)
                %% remove empty values
                %% Octave 5.2 has issues with this ...
                % k = find(cellfun(@isempty, vals));
                %% ... so we use this instead ...
                k = find(cellfun(@(x)prod(size(x)) == 0, vals));
                if ~isempty(k)
                    var_names(k) = [];
                    vals(k) = [];
                end

                if update       %% update existing dme table w/values
                    if isempty(ridx)    %% all rows
                        for k = 1:length(var_names)
                            dme.tab.(var_names{k}) = vals{k};
                        end
                    else                %% indexed rows
                        for k = 1:length(var_names)
                            dme.tab.(var_names{k})(ridx, :) = vals{k};
                        end
                    end
                else            %% create dme table from values
                    if ~isempty(vals)
                        table_class = mp_table_class();
                        dme.tab = table_class(vals{:}, 'VariableNames', var_names);
                    end

                    %% check for unique uid's if not generated
                    if strcmp(var_names{1}, 'uid') && ~strcmp(spec.vmap.uid{1}, 'IDs')
                        if length(unique(vals{1})) ~= spec.nr
                            error('mp.dmc_element.import_table_values: ''uid'' values must be unique\ndata contains only %d unique ''uid'' value(s) for %d ''%s'' elements\n', ...
                                length(unique(vals{1})), spec.nr, obj.name);
                        end
                    end
                end
            end
        end

        function vals = get_input_table_values(obj, d, spec, var_names, ridx)
            % Get values to insert in data model element table.
            % ::
            %
            %   vals = dmce.get_input_table_values(d, spec, var_names, ridx)
            %
            % Inputs:
            %   d : data source, type depends on the implementing subclass
            %       (e.g. |MATPOWER| case struct for mp.dm_converter_mpc2)
            %   spec (struct) : import specification, see get_import_spec()
            %   var_names (cell array) : *(optional)* list of names of
            %       variables (columns of main table) to import *(default is
            %       all variables)*
            %   ridx (integer) : *(optional)* vector of row indices of data
            %       to import *(default is all rows)*
            %
            % Output:
            %   vals (cell array) : values to assign to table columns in data
            %       model element
            %
            % Called by import_table_values().

            if nargin < 5
                ridx = [];
            end
            nr = spec.nr;

            if nr
                %% initialize variable values
                vals = cell(size(var_names));

                %% assign variable values
                for k = 1:length(var_names)
                    vn = var_names{k};
                    vm = spec.vmap.(vn);
                    switch vm{1}    %% switch on type of mapping
                        case {'col'}     %% column of default table
                            val = obj.import_col(d, spec, vn, vm{2:end});
                        case 'num'      %% scalar numerical value
                            val = vm{2} * ones(nr, 1);
                        case {'cell'}   %% cell array of values
                            val = cell(nr, 1);
                            [val{:}] = deal(vm{2});
                        case 'IDs'      %% [1:nr]', consecutive IDs
                            val = [1:nr]';
                        case 'r'        %% r
                            val = spec.r;
                        case 'fcn'      %% general function
                            if length(vm) > 1 && isa(vm{2}, 'function_handle')
                                import_fcn = vm{2};
                                val = import_fcn(obj, d, spec, vn);
                            end
                        otherwise
                            error('mp.dmc_element.get_input_table_values: %d is an unknown var map type', vm{1});
                    end
                    if ~isempty(ridx)
                        vals{k} = val(ridx, :);
                    else
                        vals{k} = val;
                    end
                end
            else            %% table does not exist or is empty
                vals = [];
            end
        end

        function vals = import_col(obj, d, spec, vn, c, sf)
            % Extract and optionally modify values from data source column.
            % ::
            %
            %   vals = dmce.import_col(d, spec, vn, c, sf)
            %
            % Inputs:
            %   d : data source, type depends on the implementing subclass
            %       (e.g. |MATPOWER| case struct for mp.dm_converter_mpc2)
            %   spec (struct) : import specification, see get_import_spec()
            %   vn (char array) : variable name
            %   c (integer) : column index for data in data source
            %   sf (double or function handle) : *(optional)* scale factor,
            %       function is called as `sf(dmce, vn)`
            %
            % Output:
            %   vals (cell array) : values to assign to table columns in data
            %       model element
            %
            %
            % Called by get_input_table_values().

            if nargin > 5
                if isa(sf, 'function_handle')
                    sf = sf(obj, vn);
                end
            else
                sf = 1;
            end

            %% default to zeros if mpc table does not have this many columns
            if c > spec.nc
                vals = zeros(spec.nr, 1);
            else
                data = subsref(d, spec.subs);
                if spec.nr && isempty(spec.r)
                    vals = sf * data(:, c);
                else
                    vals = sf * data(spec.r, c);
                end
            end
        end

        function d = export(obj, dme, d, var_names, ridx)
            % Export data from data model element to data source.
            % ::
            %
            %   d = dmce.export(dme, d, var_names, ridx)
            %
            % Inputs:
            %   dme (mp.dm_element) : data model element object
            %   d : data source, type depends on the implementing subclass
            %       (e.g. |MATPOWER| case struct for mp.dm_converter_mpc2)
            %   var_names (cell array) : *(optional)* list of names of
            %       variables (columns of main table) to export *(default is
            %       all variables)*
            %   ridx (integer) : *(optional)* vector of row indices of data
            %       to export *(default is all rows)*
            %
            % Output:
            %   d : updated data source
            %
            % See also import.

            if nargin < 5
                ridx = [];
                if nargin < 4
                    var_names = {};
                end
            end

            %% get conversion spec
            spec = obj.get_export_spec(dme, d);

            %% export main table
            d = obj.export_table_values(dme, d, spec, var_names, ridx);
        end

        function d = export_table_values(obj, dme, d, spec, var_names, ridx)
            % Export table values for given import specification.
            % ::
            %
            %   d = dmce.export_table_values(dme, d, spec, var_names, ridx)
            %
            % Inputs:
            %   dme (mp.dm_element) : data model element object
            %   d : data source, type depends on the implementing subclass
            %       (e.g. |MATPOWER| case struct for mp.dm_converter_mpc2)
            %   spec (struct) : export specification, see get_export_spec()
            %   var_names (cell array) : *(optional)* list of names of
            %       variables (columns of main table) to export *(default is
            %       all variables)*
            %   ridx (integer) : *(optional)* vector of row indices of data
            %       to export *(default is all rows)*
            %
            % Output:
            %   d : updated data source
            %
            % Called by export().

            if nargin < 6
                ridx = [];
                if nargin < 5
                    var_names = {};
                end
            end

            if isempty(var_names)       %% default to all variables
                var_names = dme.main_table_var_names();
            end

            if ~obj.data_exists(d)      %% initialize data field in d
                d = obj.init_export_data(dme, d, spec);
            end

            for k = 1:length(var_names)
                vn = var_names{k};
                vm = spec.vmap.(vn);
                switch vm{1}    %% switch on type of mapping
                    case {'col'}     %% column of default table
                        d = obj.export_col(dme, d, spec, vn, ridx, vm{2:end});
                    case {'num', 'cell', 'IDs', 'r'}
                        %% do nothing
                    case 'fcn'      %% general function
                        if length(vm) > 2 && isa(vm{3}, 'function_handle')
                            export_fcn = vm{3};
                            d = export_fcn(obj, dme, d, spec, vn, ridx);
                        end
                    otherwise
                        error('mp.dmc_element.export_table_values: %d is an unknown var map type', vm{1});
                end
            end
        end

        function d = init_export_data(obj, dme, d, spec)
            % Initialize data source for export from data model element.
            % ::
            %
            %   d = dmce.init_export_data(dme, d, spec)
            %
            % Inputs:
            %   dme (mp.dm_element) : data model element object
            %   d : data source, type depends on the implementing subclass
            %       (e.g. |MATPOWER| case struct for mp.dm_converter_mpc2)
            %   spec (struct) : export specification, see get_export_spec()
            %
            % Output:
            %   d : updated data source
            %
            % Called by export_table_values().

            d = subsasgn(d, spec.subs, obj.default_export_data_table(spec));
        end

        function dt = default_export_data_table(obj, spec)
            % Create default (empty) data table for data source.
            % ::
            %
            %   dt = dmce.default_export_data_table(spec)
            %
            % Input:
            %   spec (struct) : export specification, see get_export_spec()
            %
            % Output:
            %   dt : data table for data source, type depends on implementing
            %       subclass
            %
            % Called by init_export_data().

            nr = obj.default_export_data_nrows(spec);
            dt = zeros(nr, spec.nc);
        end

        function nr = default_export_data_nrows(obj, spec)
            % Get number of rows for default_export_data_table().
            % ::
            %
            %   nr = default_export_data_nrows(spec)
            %
            % Input:
            %   spec (struct) : export specification, see get_export_spec()
            %
            % Output:
            %   nr (integer) : number of rows
            %
            % Called by default_export_data_table().

            if ~isempty(spec.r)
                nr = max(spec.nr, max(spec.r));
            else
                nr = spec.nr;
            end
        end

        function d = export_col(obj, dme, d, spec, vn, ridx, c, sf)
            % Export a variable (table column) to the data source.
            % ::
            %
            %   d = dmce.export_col(dme, d, spec, vn, ridx, c, sf)
            %
            % Inputs:
            %   dme (mp.dm_element) : data model element object
            %   d : data source, type depends on the implementing subclass
            %       (e.g. |MATPOWER| case struct for mp.dm_converter_mpc2)
            %   spec (struct) : export specification, see get_export_spec()
            %   vn (char array) : variable name
            %   ridx (integer) : *(optional)* vector of row indices of data
            %       to export *(default is all rows)*
            %   c (integer) : column index for data in data source
            %   sf (double or function handle) : *(optional)* scale factor,
            %       function is called as `sf(dmce, vn)`
            %
            % Output:
            %   d : updated data source
            %
            % Called by export_table_values().

            if nargin > 7
                if isa(sf, 'function_handle')
                    sf = sf(obj, vn);
                end
            else
                sf = 1;
            end

            %% default to zeros if mpc table does not have this many columns
            ss = spec.subs;
            if sf
                ss(end+1).type = '()';
                if spec.nr && isempty(spec.r)
                    if isempty(ridx)
                        ss(end).subs = {':', c};
                    else
                        ss(end).subs = {ridx, c};
                    end
                else
                    if isempty(ridx)
                        ss(end).subs = {spec.r, c};
                    else
                        ss(end).subs = {spec.r(ridx), c};
                    end
                end
                if isempty(ridx)
                    val = dme.tab.(vn) / sf;
                else
                    val = dme.tab.(vn)(ridx, :) / sf;
                end
                d = subsasgn(d, ss, val);
            end
        end
    end     %% methods

    methods (Access=protected)
        function spec = get_spec(obj, dme, d, import)
            % Implementation of get_import_spec() and get_export_spec().
            % ::
            %
            %   spec = dmce.get_spec(dme, d, import)

            if import == 1
                [nr, nc, r] = obj.get_import_size(d);
            else
                [nr, nc, r] = obj.get_export_size(dme);
            end
            spec = struct( ...
                'subs', obj.data_subs(), ...
                'nr', nr, ...
                'nc', nc, ...
                'r',  r, ...
                'vmap', obj.table_var_map(dme, d) ...
            );
        end
    end     %% methods (Access=protected)
end         %% classdef
