classdef (Abstract) dm_element < handle
%MP.DM_ELEMENT  Abstract base class for MATPOWER data model elements

%   MATPOWER
%   Copyright (c) 2020-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        tab             %% main data table
        nr              %% total number of rows in table
        n               %% number of online elements
        ID2i            %% max(ID) x 1 vector, maps IDs to row indices
        on              %% n x 1 vector of row indices of online elements
        off             %% (nr-n) x 1 vector of row indices of offline elements
        i2on            %% nr x 1 vector mapping row index to index in on/off respectively
    end     %% properties

    methods
        function name = name(obj)
            name = '';      %% e.g. 'bus', 'gen'
        end

        function label = label(obj)
            label = '';     %% e.g. 'Bus', 'Generator'
        end

        function label = labels(obj)
            label = '';     %% e.g. 'Buses', 'Generators'
        end

        function name = cxn_type(obj)
            %% char array or cell array of char arrays with names of types of
            %% junction elements, i.e. node-creating elements (e.g. 'bus'),
            %% this element can connect to
            name = '';
        end

        function name = cxn_idx_prop(obj)
            %% char array or cell array of char arrays with names of properties
            %% containing indices of junction elements that define connections
            %% (e.g. {'bus_fr', 'bus_to'})
            name = '';
        end

        function name = cxn_type_prop(obj)
            %% char array or cell array of char arrays (dim must match
            %% cxn_idx_prop) with names of properties containing type of
            %% junction elements for each connection, only used if the
            %% junction element type can vary by element (e.g. some lines
            %% connect to one kind of bus, some to another kind)
            name = '';
        end

        function TorF = table_exists(obj)
            TorF = ~isempty(obj.tab);
        end

        function names = table_var_names(obj)
            names = obj.main_table_var_names();
        end

        function names = main_table_var_names(obj)
            names = {'uid', 'name', 'status', 'source_uid'};
        end

        function vars = export_vars(obj)
            vars = {};
        end

        function dmce = dm_converter_element(obj, dmc, name)
            if nargin < 3
                name = obj.name;
            end
            dmce = dmc.elements.(name);
        end

        function new_obj = copy(obj)
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
            if isempty(obj.tab)
                nr = 0;
            else
                nr = size(obj.tab, 1);
            end
            obj.nr = nr;
        end

        function obj = initialize(obj, dm)
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
            %% returns nr x 1 vector of unique IDs, maps row index to ID
            %% (or indexed version)
            if nargin > 1 && ~isempty(idx)
                uid = obj.tab.uid(idx);
            else
                uid = obj.tab.uid;
            end
        end

        function obj = init_status(obj, dm)
        end

        function obj = update_status(obj, dm)
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
        end

        function obj = rebuild(obj, dm)
            obj.count(dm);
            obj.initialize(dm);
            obj.init_status(dm);
            obj.build_params(dm);
        end

        function display(obj)
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
            switch section
                case {'cnt', 'sum', 'ext', 'det'}
                    rows = -1;  %% all rows
                otherwise
                    rows = obj.pp_rows_other(dm, section, out_e, mpopt, pp_args);
            end
        end

        function h = pp_get_headers(obj, dm, section, out_e, mpopt, pp_args)
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
            TorF = true;    %% all elements have count section by default
        end

        function obj = pp_data_cnt(obj, dm, rows, out_e, mpopt, fd, pp_args)
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
            TorF = false;   %% no summary section for elements by default
        end

        function obj = pp_data_sum(obj, dm, rows, out_e, mpopt, fd, pp_args)
        end

        function TorF = pp_have_section_ext(obj, mpopt, pp_args)
            TorF = false;   %% no ext section for elements by default
        end

        function obj = pp_data_ext(obj, dm, rows, out_e, mpopt, fd, pp_args)
        end

        function TorF = pp_have_section_det(obj, mpopt, pp_args)
            TorF = false;   %% no ext section for elements by default
        end

        function str = pp_get_title_det(obj, mpopt, pp_args)
            if obj.pp_have_section_det(mpopt, pp_args)
                str = sprintf('%s Data', obj.label);
            else
                str = '';
            end
        end

        function h = pp_get_headers_det(obj, dm, out_e, mpopt, pp_args)
            str = obj.pp_get_title_det(mpopt, pp_args);
            if isempty(str)
                h = {};
            else
                h = dm.pp_section_label(str);
            end
        end

        function f = pp_get_footers_det(obj, dm, out_e, mpopt, pp_args)
            f = {};
        end

        function obj = pp_data_det(obj, dm, rows, out_e, mpopt, fd, pp_args)
            if obj.pp_have_section_det(mpopt, pp_args)
                assert(rows == -1, 'dm_elemnt/pp_data_det: ''rows'' is expected to be -1, indicating all rows');
                for k = 1:obj.nr
                    fprintf(fd, '%s\n', ...
                        obj.pp_data_row_det(dm, k, out_e, mpopt, fd, pp_args));
                end
            end
        end

        function str = pp_data_row_det(obj, dm, k, out_e, mpopt, fd, pp_args)
            str = '';
        end
    end     %% methods
end         %% classdef
