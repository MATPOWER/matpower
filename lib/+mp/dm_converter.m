classdef (Abstract) dm_converter < mp.element_container
% mp.dm_converter - Abstract base class for |MATPOWER| **data model converter** objects.
%
% A data model converter provides the ability to convert data between a
% data model and a specific data source or format, such as the PSS/E RAW
% format or version 2 of the |MATPOWER| case format. It is used, for example,
% during the import stage of the data model build process.
%
% A data model converter object is primarily a container for data model
% converter element (mp.dmc_element) objects. Concrete data model converter
% classes are specific to the type or format of the data source.
%
% By convention, data model converter variables are named ``dmc`` and data
% model converter class names begin with ``mp.dm_converter``.
%
% mp.dm_converter Methods:
%   * format_tag - return char array identifier for data source/format
%   * copy - make duplicate of object
%   * build - create and add element objects
%   * import - import data from a data source into a data model
%   * export - export data from a data model to a data source
%   * init_export - initialize a data source for export
%   * save - save data source to a file
%   * display - display the data model converter object
%
% See the :ref:`sec_dm_converter` section in the |MATPOWER-Dev-Manual| for more
% information.
%
% See also mp.data_model, mp.task.

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
        function tag = format_tag(obj)
            % Return a short char array identifier for data source/format.
            % ::
            %
            %   tag = dmc.format_tag()
            %
            % E.g. the subclass for the |MATPOWER| case format
            % returns ``'mpc2'``.
            %
            % *Note: This is an abstract method that must be implemented
            % by a subclass.*

            error('mp.dm_converter.format_tag: must be implemented by sub-class');
        end

        function new_obj = copy(obj)
            % Create a duplicate of the data model converter object, calling
            % the :meth:`copy() <mp.dmc_element.copy>` method on each element.
            % ::
            %
            %   new_dmc = dmc.copy()

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

            %% make copies of each individual element
            new_obj.elements = new_obj.elements.copy();
        end

        function obj = build(obj)
            % Create and add data model converter element objects.
            % ::
            %
            %   dmc.build()
            %
            % Create the data model converter element objects by instantiating
            % each class in the :attr:`element_classes <mp.element_container.element_classes>`
            % property and adding the resulting object to the
            % :attr:`elements <mp.element_container.elements>` property.

            obj.elements = mp.mapped_array();
            for c = obj.element_classes
                dmce = c{1}();      %% element constructor
                obj.elements.add_elements(dmce, dmce.name);
            end
        end

        function dm = import(obj, dm, d)
            % Import data from a data source into a data model.
            % ::
            %
            %   dm = dmc.import(dm, d)
            %
            % Inputs:
            %   dm (mp.data_model) : data model
            %   d : data source, type depends on the implementing subclass
            %       (e.g. |MATPOWER| case struct for mp.dm_converter_mpc2)
            %
            % Output:
            %   dm (mp.data_model) : updated data model
            %
            % Calls the import() method for each data model converter element
            % and its corresponding data model element.

            for k = 1:length(obj.elements)
                dmce = obj.elements{k};
                dme = dmce.data_model_element(dm);
                dme = dmce.import(dme, d);
            end
        end

        function d = export(obj, dm, d)
            % Export data from a data model to a data source.
            % ::
            %
            %   d = dmc.export(dm, d)
            %
            % Inputs:
            %   dm (mp.data_model) : data model
            %   d : data source, type depends on the implementing subclass
            %       (e.g. |MATPOWER| case struct for mp.dm_converter_mpc2)
            %
            % Output:
            %   d : updated data source
            %
            % Calls the export() method for each data model converter element
            % and its corresponding data model element.

            if nargin < 3 || isempty(d)
                %% initialize empty data struct
                d = obj.init_export(dm);
                all_vars = 1;
            else
                all_vars = 0;
            end

            for k = 1:length(dm.elements)
                dme = dm.elements{k};
                if obj.elements.has_name(dme.name)
                    dmce = dme.dm_converter_element(obj);
                    if all_vars
                        vars = dme.main_table_var_names();
                    else
                        vars = dme.export_vars();
                    end
                    d = dmce.export(dme, d, vars);
                end
            end
        end

        function d = init_export(obj, dm)
            % Initialize a data source for export.
            % ::
            %
            %   d = dmc.export(dm)
            %
            % Input:
            %   dm (mp.data_model) : data model
            %
            % Output:
            %   d : new empty data source, type depends on the implementing
            %       subclass (e.g. |MATPOWER| case struct for
            %       mp.dm_converter_mpc2)
            %
            % Creates a new data source of the appropriate type in
            % preparation for calling export().

            d = struct();
        end

        function fname_out = save(obj, fname, d)
            % Save data source to a file.
            % ::
            %
            %   fname_out = dmc.save(fname, d)
            %
            % Inputs:
            %   fname (char array) - file name to use for saving
            %   d : data source, type depends on the implementing subclass
            %       (e.g. |MATPOWER| case struct for mp.dm_converter_mpc2)
            %
            % Output:
            %   fname_out (char array) : final file name after saving,
            %       possibly modified from input (e.g. extension added)
            %
            % *Note: This is an abstract method that must be implemented
            % by a subclass.*

            error('mp.dm_converter.save: must be implemented by subclass');
        end

        function display(obj)
            % Display the data model converter object.
            %
            % This method is called automatically when omitting a semicolon
            % on a line that retuns an object of this class.
            %
            % Displays the details of the data model converter elements.

%             if have_feature('octave')
%                 struct(obj)
%             else
%                 display@handle(obj)
%             end
            fprintf('DATA CONVERTER CLASS : %s\n', class(obj));

            %% elements
            fprintf('\nELEMENTS\n')
            fprintf('========\n')
            fprintf(' i  name          class\n');
            fprintf('-- -----------   --------------------\n');
            for k = 1:length(obj.elements)
                dmce = obj.elements{k};
                fprintf('%2d  %-13s %s\n', k, dmce.name, class(dmce));
%                 dmce
            end
        end
    end     %% methods
end         %% classdef
