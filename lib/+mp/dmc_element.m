classdef (Abstract) dmc_element < handle
%MP.DMC_ELEMENT  Abstract base class for data model converter for indv elements

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties


    methods
        function name = name(obj)
            name = '';      %% e.g. 'bus', 'gen'
        end

        function dme = data_model_element(obj, dm, name)
            if nargin < 3
                name = obj.name;
            end
            dme = dm.elements.(name);
        end

        function df = data_field(obj)
            df = '';        %% name of field in d for default data table
        end

        function s = data_subs(obj)
            s = struct('type', '.', 'subs', obj.data_field());
        end

        function TorF = data_exists(obj, d)
            TorF = isfield(d, obj.data_field());
        end

        function spec = get_import_spec(obj, dme, d)
            spec = obj.get_spec(dme, d, 1);
        end

        function spec = get_export_spec(obj, dme, d)
            spec = obj.get_spec(dme, d, 0);
        end

        function spec = get_spec(obj, dme, d, import)
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
                'vmap', obj.table_var_map(dme, d), ...
                'tab_name', 'tab' ...
            );
        end

        function [nr, nc, r] = get_import_size(obj, d)
            if obj.data_exists(d)
                %% use size of default table
                [nr, nc] = size(subsref(d, obj.data_subs()));
            else
                [nr, nc] = deal(0);
            end
            r = [];                         %% all rows
        end

        function [nr, nc, r] = get_export_size(obj, dme)
            [nr, nc] = size(dme.tab);   %% use size of default table
            r = [];                     %% all rows
        end

        function vmap = table_var_map(obj, dme, d)
            %% initialize with vmap.(<name>) = {'col', []}, for all <name>
            names = dme.table_var_names();
            vals = cell(size(names));
            [vals{:}] = deal({'col', []});
            vmap = cell2struct(vals, names, 2);
        end

        function dme = import(obj, dme, d, var_names, ridx)
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
                            error('mp.dmc_element/import_table_values: ''uid'' values must be unique\ndata contains only %d unique ''uid'' value(s) for %d ''%s'' elements\n', ...
                                length(unique(vals{1})), spec.nr, obj.name);
                        end
                    end
                end
            end
        end

        function vals = get_input_table_values(obj, d, spec, var_names, ridx)
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
                            error('mp.dmc_element/get_input_table_values: %d is an unknown var map type', vm{1});
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
                        error('mp.dmc_element/export: %d is an unknown var map type', vm{1});
                end
            end
        end

        function d = init_export_data(obj, dme, d, spec)
            d = subsasgn(d, spec.subs, obj.default_export_data_table(spec));
        end

        function dt = default_export_data_table(obj, spec)
            nr = obj.default_export_data_nrows(spec);
            dt = zeros(nr, spec.nc);
        end

        function nr = default_export_data_nrows(obj, spec)
            if ~isempty(spec.r)
                nr = max(spec.nr, max(spec.r));
            else
                nr = spec.nr;
            end
        end

        function d = export_col(obj, dme, d, spec, vn, ridx, c, sf)
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
end         %% classdef
