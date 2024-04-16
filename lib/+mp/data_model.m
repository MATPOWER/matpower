classdef data_model < mp.element_container
% mp.data_model - Base class for |MATPOWER| **data model** objects.
%
% The data model object encapsulates the input data provided by the user
% for the problem of interest and the output data presented back to the user
% upon completion. It corresponds roughly to the ``mpc`` (|MATPOWER| case)
% and ``results`` structs used throughout the legacy |MATPOWER| implementation,
% but encapsulated in an object with additional functionality. It includes
% tables of data for each type of element in the system.
%
% A data model object is primarily a container for data model element
% (mp.dm_element) objects. Concrete data model classes may be specific to
% the task.
%
% By convention, data model variables are named ``dm`` and data model class
% names begin with ``mp.data_model``.
%
% mp.data_model Properties:
%   * base_mva - system per unit MVA base
%   * base_kva - system per unit kVA base
%   * source - source of data, e.g. ``mpc`` (|MATPOWER| case struct)
%   * userdata - arbitrary user data
%
% mp.data_model Methods:
%   * data_model - constructor, assign  default data model element classes
%   * copy - make duplicate of object
%   * build - create, add, and build element objects
%   * count - count instances of each element and remove if count is zero
%   * initialize - initialize (online/offline) status of each element
%   * update_status - update (online/offline) status based on connectivity, etc
%   * build_params - extract/convert/calculate parameters for online elements
%   * online - get number of online elements of named type
%   * display - display the data model object
%   * pretty_print - pretty print data model to console or file
%   * pp_flags - from options, build flags to control pretty printed output
%   * pp_section_label - construct section header lines for output
%   * pp_section_list - return list of section tags
%   * pp_have_section - return true if section exists for object
%   * pp_section - pretty print the given section
%   * pp_get_headers - construct pretty printed lines for section headers
%   * pp_get_headers_cnt - construct pretty printed lines for **cnt** section headers
%   * pp_get_headers_ext - construct pretty printed lines for **ext** section headers
%   * pp_data - pretty print the data for the given section
%   * set_bus_v_lims_via_vg - set gen bus voltage limits based on gen voltage setpoints
%
% See the :ref:`sec_data_model` section in the |MATPOWER-Dev-Manual| for more
% information.
%
% See also mp.task, mp.net_model, mp.math_model, mp.dm_converter.

%   MATPOWER
%   Copyright (c) 1996-2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        % *(double)* system per unit MVA base, for balanced single-phase
        % systems/sections, must be provided if system includes any
        % ``'bus'`` elements
        base_mva

        % *(double)* system per unit kVA base, for unbalanced 3-phase
        % systems/sections, must be provided if system includes any
        % ``'bus3p'`` elements
        base_kva
        source      % source of data, e.g. ``mpc`` (|MATPOWER| case struct)
        userdata = struct();    % *(struct)* arbitrary user data
    end     %% properties

    methods
        %% constructor
        function obj = data_model()
            % Constructor, assign default data model element classes.
            % ::
            %
            %   dm = mp.data_model()

            %% call parent constructor
            obj@mp.element_container();
            obj.element_classes = ...
                { @mp.dme_bus, @mp.dme_gen, @mp.dme_load, ...
                    @mp.dme_branch, @mp.dme_shunt };
        end

        function new_obj = copy(obj)
            % Create a duplicate of the data model object, calling the
            % :meth:`copy() <mp.dm_element.copy>` method on each element.
            % ::
            %
            %   new_dm = dm.copy()
        
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

            %% make copies of each individual element
            new_obj.elements = new_obj.elements.copy();
        end

        function obj = build(obj, d, dmc)
            % Create and add data model element objects.
            % ::
            %
            %   dm.build(d, dmc)
            %
            % Inputs:
            %   d : data source, type depends on the implementing subclass
            %       (e.g. |MATPOWER| case struct for mp.dm_converter_mpc2)
            %   dmc (mp.dm_converter) : data model converter
            %
            % Create the data model element objects by instantiating each
            % class in the :attr:`element_classes <mp.element_container.element_classes>`
            % property and adding the resulting object to the
            % :attr:`elements <mp.element_container.elements>` property.
            % Then proceed through the following additional build stages for
            % each element.
            %
            %   - Import
            %   - Count
            %   - Initialize
            %   - Update status
            %   - Build parameters
            %
            % See the :ref:`sec_building_data_model` section in the
            % |MATPOWER-Dev-Manual| for more information.

            %% create empty element objects for each class
            obj.elements = mp.mapped_array();
            for k = 1:length(obj.element_classes)
                dme_class = obj.element_classes{k};
                dme = dme_class();      %% element constructor
                obj.elements.add_elements(dme, dme.name);
            end

            dmc.import(obj, d);     %% import data from external format
            obj.count();            %% count element objects for each class
            obj.initialize();       %% get IDs, self status
            obj.update_status();    %% update status wrt connectivity, define on/off
            obj.build_params();     %% get parameters
        end

        function obj = count(obj)
            % Count instances of each element and remove if count is zero.
            % ::
            %
            %   dm.count()
            %
            % Call each element's :meth:`count() <mp.dm_element.count>` method
            % to determine the number of instances of that element in the data,
            % and remove the element type from
            % :attr:`elements <mp.element_container.elements>` if the count
            % is 0.
            %
            % Called by build() to perform its **count** stage. See the
            % :ref:`sec_building_data_model` section in the
            % |MATPOWER-Dev-Manual| for more information.

            for k = length(obj.elements):-1:1
                obj.elements{k}.count(obj);     %% set nr
                if obj.elements{k}.nr == 0
                    obj.elements.delete_elements(k);
                end
            end
        end

        function obj = initialize(obj)
            % Initialize (online/offline) status of each element.
            % ::
            %
            %   dm.initialize()
            %
            % Call each element's :meth:`initialize() <mp.dm_element.initialize>`
            % method to initialize statuses and create ID to row index
            % mappings.
            %
            % Called by build() to perform its **initialize** stage. See the
            % :ref:`sec_building_data_model` section in the
            % |MATPOWER-Dev-Manual| for more information.

            for k = 1:length(obj.elements)
                obj.elements{k}.initialize(obj);
            end
        end

        function obj = update_status(obj)
            % Update (online/offline) status based on connectivity, etc.
            % ::
            %
            %   dm.update_status()
            %
            % Call each element's :meth:`update_status() <mp.dm_element.update_status>`
            % method to update statuses based on connectivity or other criteria
            % and define element properties containing number and row indices
            % of online elements, indices of offline elements, and mapping of
            % row indices to indices in online and offline element lists.
            %
            % Called by build() to perform its **update status** stage. See the
            % :ref:`sec_building_data_model` section in the
            % |MATPOWER-Dev-Manual| for more information.

            for k = 1:length(obj.elements)
                obj.elements{k}.update_status(obj);
            end
        end

        function obj = build_params(obj)
            % Extract/convert/calculate parameters for online elements.
            % ::
            %
            %   dm.build_params()
            %
            % Call each element's :meth:`build_params() <mp.dm_element.build_params>`
            % method to build parameters as necessary for online elements from
            % the original data tables (e.g. p.u. conversion, initial
            % state, etc.) and store them in element-specific properties.
            %
            % Called by build() to perform its **build parameters** stage. See
            % the :ref:`sec_building_data_model` section in the
            % |MATPOWER-Dev-Manual| more information.

            for k = 1:length(obj.elements)
                obj.elements{k}.build_params(obj);
            end
        end

        function n = online(obj, name)
            % Get number of online elements of named type.
            % ::
            %
            %   n = dm.online(name)
            %
            % Input:
            %   name (char array) : name of element type (e.g. ``'bus'``,
            %       ``'gen'``) as returned by the element's
            %       :meth:`name() <mp.dm_element.name>` method
            %
            % Output:
            %   n (integer) : number of online elements
            
            if obj.elements.has_name(name)
                n = obj.elements.(name).n;
            else
                n = 0;
            end
        end

        function display(obj)
            % Display the data model object.
            %
            % This method is called automatically when omitting a semicolon
            % on a line that retuns an object of this class.
            %
            % Displays the details of the data model elements.

%             if have_feature('octave')
%                 struct(obj)
%             else
%                 display@handle(obj)
%             end
            fprintf('DATA MODEL CLASS : %s\n', class(obj));
            fprintf('        base_mva : %g\n', obj.base_mva);

            %% elements
            fprintf('\nELEMENTS\n')
            fprintf('========\n')
            fprintf(' i  name            nr         n      class\n');
            fprintf('-- -----------   --------  --------  --------------------\n');
            for k = 1:length(obj.elements)
                dme = obj.elements{k};
                fprintf('%2d  %-13s %6d %9d    %s\n', k, dme.name, dme.nr, dme.n, class(dme));
%                 dme
            end

            %% user data
            fields = fieldnames(obj.userdata);
            if ~isempty(fields)
                fprintf('\nUSER DATA\n')
                fprintf('=========\n')
                fprintf('  name                         size       class\n');
                fprintf(' ------------------------   -----------  --------------------\n');
                for k = 1:length(fields)
                    f = obj.userdata.(fields{k});
                    [m, n] = size(f);
                    fprintf('  %-24s  %5dx%-5d   %s\n', fields{k}, m, n, class(f));
                end
            end
        end

        function [obj, out_] = pretty_print(obj, mpopt, fd)
            % Pretty print data model to console or file.
            % ::
            %
            %   dm.pretty_print(mpopt)
            %   dm.pretty_print(mpopt, fd)
            %   [dm, out] = dm.pretty_print(mpopt, fd)
            %
            % Inputs:
            %   mpopt (struct) : |MATPOWER| options struct
            %   fd (integer) : *(optional, default = 1)* file identifier to
            %       use for printing, (1 for standard output, 2 for standard
            %       error)
            %
            % Outputs:
            %   dm (mp.data_model) : the data model object
            %   out (struct) : struct of output control flags
            %
            % Displays the model parameters to a pretty-printed text format.
            % The result can be output either to the console or to a file.
            %
            % The output is organized into sections and each element type
            % controls its own output for each section. The default sections
            % are:
            %
            %   - **cnt** - counts, number of online, offline, and total
            %     elements of this type
            %   - **sum** - summary, e.g. total amount of capacity, load,
            %     line loss, etc.
            %   - **ext** - extremes, e.g. min and max voltages, nodal prices, etc.
            %   - **det** - details, table of detailed data, e.g. voltages, prices
            %     for buses, dispatch, limits for generators, etc.

            if nargin < 3
                fd = 1;     %% print to stdio by default
                if nargin < 2
                    mpopt = mpoption();
                end
            end

            %% get output flags
            out = obj.pp_flags(mpopt);

            if out.any
                sections = obj.pp_section_list(out);  %% e.g. cnt, sum, ext, det, etc.
                for s = 1:length(sections)
                    out_s = out.sec.(sections{s});
                    if out_s.any && obj.pp_have_section(sections{s}, mpopt)
                        obj.pp_section(sections{s}, out_s, mpopt, fd);
                    end
                end
            end

            %% return output control flags
            if nargout > 1
                out_ = out;
            end
        end

        function [out, add] = pp_flags(obj, mpopt)
            % From options, build flags to control pretty printed output.
            % ::
            %
            %   [out, add] = dm.pp_flags(mpopt)
            %
            % Input:
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Outputs:
            %   out (struct) : struct of output control flags
            %       ::
            %
            %           out
            %             .all    (-1, 0 or 1)
            %             .any    (0 or 1)
            %             .sec
            %               .cnt
            %                 .all    (-1, 0 or 1)
            %                 .any    (0 or 1)
            %               .sum    (same as cnt)
            %               .ext    (same as cnt)
            %               .det
            %                 .all    (-1, 0 or 1)
            %                 .any    (0 or 1)
            %                 .elm
            %                   .<name>    (0 or 1)
            %
            %       where <name> is the name of the corresponding element
            %       type.
            %
            %   add (struct) : additional data for subclasses to use
            %       ::
            %
            %           add
            %             .s0
            %               .<name> = 0
            %             .s1
            %               .<name> = 1
            %             .suppress    (-1, 0 or 1)
            %             .names    (cell array of element names)
            %             .ne       (number of element names)
            %
            % See also pretty_print.

            %% output control flags
            ne = length(obj.elements);
            names = cellfun(@(k)obj.elements{k}.name, num2cell(1:ne), ...
                            'UniformOutput', false );

            %% suppress output detail (default for large systems)
            suppress = mpopt.out.suppress_detail;
            if suppress == -1
                if obj.elements.has_name('bus') && ...
                        obj.elements.bus.nr > 500
                    suppress = 1;
                else
                    suppress = 0;
                end
            end

            %% element struct with all fields == 0, names are dm element names
            s0 = cell2struct(num2cell(zeros(1,ne)), names, 2);
            s1 = cell2struct(num2cell(ones( 1,ne)), names, 2);
            s2 = struct('all', 1, 'any', 1);

            out = struct( ...
                'all', mpopt.out.all, ...
                'any', 0, ...
                'sec', struct( ...
                    'cnt', s2, ...
                    'sum', s2, ...
                    'ext', s2, ...
                    'det', struct( ...
                        'all', -1, ...
                        'any', 0, ...
                        'elm', s1 ) ) );
            out.sec.cnt.all = out.all == 1 || (out.all == -1 && mpopt.out.sys_sum);
            out.sec.sum.all = out.sec.cnt.all;
            out.sec.ext.all = out.sec.cnt.all;
            out.sec.cnt.any = out.sec.cnt.all;
            out.sec.sum.any = out.sec.cnt.all;
            out.sec.ext.any = out.sec.cnt.all;
%             out.area_sum = out.all == 1 || (out.all == -1 && ~suppress && mpopt.out.area_sum);

            %% update detail options
            for k = 1:ne
                out.sec.det.elm.(names{k}) = out.all == 1 || ...
                    ( out.all == -1 && ~suppress && ...
                        (~isfield(mpopt.out, names{k}) || mpopt.out.(names{k})) );
            end
            out.sec.det.any = any(cell2mat( ...
                    cellfun(@double, struct2cell(out.sec.det.elm), ...
                            'UniformOutput', false) ));

            %% update any field
            out.any = out.all == 1 || ...
                        (out.all == -1 && (out.sec.cnt.all || out.sec.det.any ));

            %% return additional data
            if nargout > 1
                add = struct('s1', s1, 's0', s0, 'suppress', suppress, ...
                    'names', {names}, 'ne', ne);
            end
        end

        function h = pp_section_label(obj, label, blank_line)
            % Construct pretty printed lines for section label.
            % ::
            %
            %   h = dm.pp_section_label(label, blank_line)
            %
            % Inputs:
            %   label (char array) : label for the section header
            %   blank_line (boolean) : include a blank line before the section
            %       label if true
            %
            % Output:
            %   h (cell array of char arrays) : individual lines of section label
            %
            % See also pretty_print.
            
            if nargin < 3
                blank_line = 1; %% include a blank line before the section label
            end

            width = 80;
            line1 = repmat('=', 1, width);
            line2 = ['|     ' label repmat(' ', 1, width-length(label)-7) '|'];
            if blank_line
                h = {'', line1, line2, line1};
            else
                h = {line1, line2, line1};
            end
        end

        function sections = pp_section_list(obj, out)
            % Return list of section tags.
            % ::
            %
            %   sections = dm.pp_section_list(out)
            %
            % Input:
            %   out (struct) : struct of output control flags (see pp_flags()
            %       for details)
            %
            % Output:
            %   sections (cell array of char arrays) : e.g. ``{'cnt', 'sum',
            %       'ext', 'det'}``
            %
            % See also pretty_print.

            sections = {'cnt', 'sum', 'ext', 'det'};
        end

        function TorF = pp_have_section(obj, section, mpopt)
            % Return true if section exists for object with given options.
            % ::
            %
            %   TorF = dm.pp_have_section(section, mpopt)
            %
            % Inputs:
            %   section (char array) : e.g. ``'cnt'``, ``'sum'``,
            %       ``'ext'``, or ``'det'``
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Output:
            %   TorF (boolean) : true if section exists
            %
            % See also pretty_print.

            %% section per-element info
            for k = 1:length(obj.elements)
                dme = obj.elements{k};
                TorF = dme.pp_have_section(section, mpopt, struct());
                if TorF
                    break;
                end
            end
        end

        function obj = pp_section(obj, section, out_s, mpopt, fd)
            % Pretty print the given section.
            % ::
            %
            %   dm.pp_section(section, out_s, mpopt, fd)
            %
            % Inputs:
            %   section (char array) : e.g. ``'cnt'``, ``'sum'``,
            %       ``'ext'``, or ``'det'``
            %   out_s (struct) : output control flags for the section,
            %       ``out_s = out.sec.(section)``
            %   mpopt (struct) : |MATPOWER| options struct
            %   fd (integer) : *(optional, default = 1)* file identifier to
            %       use for printing, (1 for standard output, 2 for standard
            %       error)
            %
            % See also pretty_print.

            %% section title & headers
            h = obj.pp_get_headers(section, out_s, mpopt);
            for k = 1:length(h)
                fprintf(fd, '%s\n', h{k});
            end

            %% section system info
            obj.pp_data(section, out_s, mpopt, fd);

            %% section per-element info
            for k = 1:length(obj.elements)
                dme = obj.elements{k};
                if out_s.all == -1
                    out_e = out_s.elm.(dme.name);
                else
                    out_e = out_s.all;
                end
                pp_args = struct('model', mpopt.model);
                dme.pretty_print(obj, section, out_e, mpopt, fd, pp_args);
            end
        end

        function h = pp_get_headers(obj, section, out_s, mpopt)
            % Construct pretty printed lines for section headers.
            % ::
            %
            %   h = dm.pp_get_headers(section, out_s, mpopt)
            %
            % Inputs:
            %   section (char array) : e.g. ``'cnt'``, ``'sum'``,
            %       ``'ext'``, or ``'det'``
            %   out_s (struct) : output control flags for the section,
            %       ``out_s = out.sec.(section)``
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Output:
            %   h (cell array of char arrays) : individual lines of section
            %       headers
            %
            % See also pretty_print.

            switch section
                case 'cnt'
                    h = obj.pp_get_headers_cnt(out_s, mpopt);
                case 'sum'
                    h = {''};
                case 'ext'
                    h = obj.pp_get_headers_ext(out_s, mpopt);
                case 'det'
                    h = {};
                otherwise
                    h = obj.pp_get_headers_other(section, out_s, mpopt);
            end
        end

        function h = pp_get_headers_cnt(obj, out_s, mpopt)
            % Construct pretty printed lines for **cnt** section headers.
            % ::
            %
            %   h = dm.pp_get_headers_cnt(out_s, mpopt)
            %
            % Inputs:
            %   out_s (struct) : output control flags for the section,
            %       ``out_s = out.sec.(section)``
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Output:
            %   h (cell array of char arrays) : individual lines of **cnt**
            %       section headers
            %
            % See also pretty_print, pp_get_headers.

            h = [ obj.pp_section_label('System Summary', 0) ...
                 {  '  elements                on     off    total', ...
                    ' --------------------- ------- ------- -------' } ];
        end

        function h = pp_get_headers_ext(obj, out_s, mpopt)
            % Construct pretty printed lines for **ext** section headers.
            % ::
            %
            %   h = dm.pp_get_headers_cnt(out_s, mpopt)
            %
            % Inputs:
            %   out_s (struct) : output control flags for the section,
            %       ``out_s = out.sec.(section)``
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Output:
            %   h (cell array of char arrays) : individual lines of **ext**
            %       section headers
            %
            % See also pretty_print, pp_get_headers.

            h = {   '', ...
                    '                                           minimum                        maximum', ...
                    '                               -----------------------------  -----------------------------' };
        end

        function h = pp_get_headers_other(obj, section, out_s, mpopt)
            % Construct pretty printed lines for other section headers.
            %
            % Returns nothing in base class, but subclasses can implement
            % other section types (e.g. ``'lim'`` for OPF).
            % ::
            %
            %   h = dm.pp_get_headers_other(section, out_s, mpopt)
            %
            % Inputs:
            %   section (char array) : e.g. ``'cnt'``, ``'sum'``,
            %       ``'ext'``, or ``'det'``
            %   out_s (struct) : output control flags for the section,
            %       ``out_s = out.sec.(section)``
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Output:
            %   h (cell array of char arrays) : individual lines of **ext**
            %       section headers
            %
            % See also pretty_print, pp_get_headers.

            h = {};
        end

        function obj = pp_data(obj, section, out_s, mpopt, fd)
            % Pretty print the data for the given section.
            % ::
            %
            %   dm.pp_data(section, out_s, mpopt, fd)
            %
            % Inputs:
            %   section (char array) : e.g. ``'cnt'``, ``'sum'``,
            %       ``'ext'``, or ``'det'``
            %   out_s (struct) : output control flags for the section,
            %       ``out_s = out.sec.(section)``
            %   mpopt (struct) : |MATPOWER| options struct
            %   fd (integer) : *(optional, default = 1)* file identifier to
            %       use for printing, (1 for standard output, 2 for standard
            %       error)
            %
            % See also pretty_print, pp_section.
        end

        function obj = set_bus_v_lims_via_vg(obj, use_vg)
            % Set gen bus voltage limits based on gen voltage setpoints.
            % ::
            %
            %   dm.set_bus_v_lims_via_vg(use_vg)
            %
            % Input:
            %   use_vg (double) : 1 if voltage setpoint should be used,
            %       0 for original bus voltage bounds, or fractional value
            %       between 0 and 1 for bounds interpolated between the two.

            bus_dme = obj.elements.bus;
            gen_dme = obj.elements.gen;
            gbus = bus_dme.i2on(gen_dme.bus(gen_dme.on));   %% buses of online gens
            nb = bus_dme.n;
            ng = gen_dme.n;

            %% gen connection matrix, element i, j is 1 if, generator j at bus i is ON
            Cg = sparse(gbus, (1:ng)', 1, nb, ng);
            Vbg = Cg * sparse(1:ng, 1:ng, gen_dme.vm_setpoint, ng, ng);
            vm_ub = max(Vbg, [], 2);    %% zero for non-gen buses, else max vm_setpoint of gens @ bus
            ib = find(vm_ub);               %% buses with online gens
            vm_lb = max(2*Cg - Vbg, [], 2); %% same as vm_ub, except min vm_setpoint of gens @ bus
            vm_lb(ib) = 2 - vm_lb(ib);

            if use_vg == 1      %% use vm_setpoint directly
                bus_dme.vm_ub(ib) = vm_ub(ib);  %% ub set by max vm_setpoint @ bus
                bus_dme.vm_lb(ib) = vm_lb(ib);  %% lb set by min vm_setpoint @ bus
                bus_dme.vm_start(ib) = vm_ub(ib);
            elseif use_vg > 0 && use_vg < 1     %% fractional value
                %% use weighted avg between original vm_lb/vm_ub limits and vm_setpoint
                bus_dme.vm_ub(ib) = (1-use_vg) * bus_dme.vm_ub(ib) + use_vg * vm_ub(ib);
                bus_dme.vm_lb(ib) = (1-use_vg) * bus_dme.vm_lb(ib) + use_vg * vm_lb(ib);
            else
                error('mp.data_model.set_bus_v_lims_via_vg: option ''opf.use_vg'' (= %g) cannot be negative or greater than 1', use_vg);
            end

            %% update bus table as well (for output)
            bus_dme.tab.vm_ub(bus_dme.on(ib)) = bus_dme.vm_ub(ib);
            bus_dme.tab.vm_lb(bus_dme.on(ib)) = bus_dme.vm_lb(ib);
        end
    end     %% methods
end         %% classdef
