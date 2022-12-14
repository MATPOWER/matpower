classdef dmce_test_widget < mp.dmc_element
%MP.DMCE_TEST_WIDGET  Data model converter for test widget element

%   MATPOWER
%   Copyright (c) 2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        use_r = 0;
    end     %% properties

    methods
        function name = name(obj)
            name = 'test_widget';
        end

        function df = data_field(obj)
            df = 'foo';
        end

        function [nr, nc, r] = get_import_size(obj, d)
            [nr, nc, r] = get_import_size@mp.dmc_element(obj, d);
            if obj.use_r
                r = find(d.foo(:, 6) ~= 0);
                nr = length(r);
            end
        end

        function [nr, nc, r] = get_export_size(obj, dme)
            [nr, nc, r] = get_export_size@mp.dmc_element(obj, dme);
            if obj.use_r
                r = dme.tab.source_uid;
                nr = length(r);
            end
        end

        function vmap = table_var_map(obj, dme, d)
            vmap = table_var_map@mp.dmc_element(obj, dme, d);

            if obj.use_r
                src = {'r'};
            else
                src = {'cell', ''}; %% empty char
            end
            sf = @(ob, vn)scale_fcn(ob, vn);
            import_fcn = @(ob, d, spec, vn)widget_name_import(ob, d, spec, vn);
            export_fcn = @(ob, dme, d, spec, vn, ridx)widget_name_export(ob, dme, d, spec, vn, ridx);

            %% mapping for each name, default is {'col', []}
            vmap.uid{2}         = 3;
            vmap.name           = {'fcn', import_fcn, export_fcn};  %% fcns to import/export from/to d.bar.baz
            vmap.status         = {'num', 1};   %% fcn w/logic for d.bus_types
            vmap.source_uid     = src;
            vmap.ids            = {'IDs'};
            vmap.twos           = {'num', 2};
            vmap.alpha{2}       = 1;
            vmap.beta           = {'col', 2, 10};
            vmap.gamma          = {'col', 4, sf};
            vmap.delta          = {'col', 5, sf};
            vmap.epsilon{2}     = 6;
        end

        function d = init_export_data(obj, dme, d, spec)
            d = init_export_data@mp.dmc_element(obj, dme, d, spec);
            if ~isfield(d, 'bar') || ~isfield(d.bar, 'baz')
                nr = obj.default_export_data_nrows(spec);
                d.bar.baz = cell(nr, 2);
            end
        end

        function dt = default_export_data_table(obj, spec)
            nr = obj.default_export_data_nrows(spec);
            dt = zeros(nr, 6);
        end

        function sf = scale_fcn(obj, vn)
            switch vn
                case 'gamma'
                    sf = 2;
                case 'delta'
                    sf = 5;
            end
        end

        function vals = widget_name_import(obj, d, spec, vn)
            if spec.nr && isempty(spec.r)
                vals = d.bar.baz(:, 2);
            else
                vals = d.bar.baz(spec.r, 2);
            end
        end

        function d = widget_name_export(obj, dme, d, spec, vn, ridx)
            if spec.nr && isempty(spec.r)
                if isempty(ridx)
                    d.bar.baz(:, 1) = dme.tab.name;
                else
                    d.bar.baz(ridx, 1) = dme.tab.name(ridx);
                end
            else
                if isempty(ridx)
                    d.bar.baz(spec.r, 1) = dme.tab.name;
                else
                    d.bar.baz(spec.r(ridx), 1) = dme.tab.name(ridx);
                end
            end
        end

        function vals = bus_status_import(obj, d, spec, vn, c)
            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE] = idx_bus;

            if spec.nr && isempty(spec.r)
                vals = d.bus(:, BUS_TYPE) ~= NONE;
            else
                vals = d.bus(spec.r, BUS_TYPE) ~= NONE;
            end
        end
    end     %% methods
end         %% classdef
