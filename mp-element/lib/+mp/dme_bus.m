classdef dme_bus < mp.dm_element
%MP.DME_BUS  MATPOWER data model class for bus data

%   MATPOWER
%   Copyright (c) 2020-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        type        %% node type vector for buses that are on
        vm_start    %% initial voltage magnitudes (p.u.) for buses that are on
        va_start    %% initial voltage angles (radians) for buses that are on
        vm_lb       %% voltage magnitude lower bounds for buses that are on
        vm_ub       %% voltage magnitude upper bounds for buses that are on
    end     %% properties

    methods
        function name = name(obj)
            name = 'bus';
        end

        function label = label(obj)
            label = 'Bus';
        end

        function label = labels(obj)
            label = 'Buses';
        end

        function names = main_table_var_names(obj)
            names = horzcat( main_table_var_names@mp.dm_element(obj), ...
                {'base_kv', 'type', 'area', 'zone', 'vm_lb', 'vm_ub', ...
                 'va', 'vm'});
        end

        function vars = export_vars(obj)
            vars = {'type', 'vm', 'va'};
        end

        function obj = init_status(obj, dm)
            %% overrides mp.dm_element/init_status()

            %% check that all buses have a valid type
            bt = obj.tab.type;
            err = find(~mp.NODE_TYPE.is_valid(bt));
            if ~isempty(err)
                error('mp.dme_bus/init_status: bus %d has an invalid type', err);
            end

            %% temporarily set bus type property with dimensions for all buses
            %% (reduced for online buses only in update_status())
            obj.type = bt;
            status = (bt ~= mp.NODE_TYPE.NONE);     %% bus status
            obj.tab.status = status;
        end

        function obj = update_status(obj, dm)
            %% call parent to fill in on/off
            update_status@mp.dm_element(obj, dm);

            %% update bus type property to correspond to online buses only
            obj.type = obj.type(obj.on);
        end

        function [gbus, ig] = gbus_vector(obj, gen_dme)
            %% buses of online gens
            gbus = obj.i2on(gen_dme.bus(gen_dme.on));
            ig = [];
        end

        function obj = set_vm_start(obj, gen_dme, gbus, ig)
            gen = gen_dme.tab;

            %% pull PV bus voltage magnitudes from gen.vm_setpoint
            vcb = ones(obj.n, 1);   %% create mask of voltage-controlled buses
            vcb(obj.type == mp.NODE_TYPE.PQ) = 0;   %% exclude PQ buses
            %% find indices of online gens at online v-c buses
            k = find(vcb(gbus));
            if isempty(ig)
                obj.vm_start(gbus(k)) = gen.vm_setpoint(gen_dme.on(k));
            else
                obj.vm_start(gbus(k)) = gen.vm_setpoint(gen_dme.on(ig(k)));
            end
        end

        function obj = build_params(obj, dm)
            %% initialize voltage from bus table
            bus = obj.tab;
            obj.va_start = bus.va(obj.on) * pi/180;
            obj.vm_start = obj.tab.vm(obj.on);
            obj.vm_lb = bus.vm_lb(obj.on);
            obj.vm_ub = bus.vm_ub(obj.on);

            %% update bus type and starting vm based on connected gen status
            if dm.elements.is_index_name('gen')
                gen_dme = dm.elements.gen;
                [gbus, ig] = obj.gbus_vector(gen_dme);
                nb = obj.n;
                ng = length(gbus);

                %% update bus types based on connected generator status
                %% gen connection matrix, element i, j is 1 if gen j @ bus i is ON
                Cg = sparse(gbus, (1:ng)', 1, nb, ng);
                bus_gen_status = Cg * ones(ng, 1);  %% num of gens ON at each bus
%                obj.type(obj.type == mp.NODE_TYPE.REF & ~bus_gen_status) = ...
%                    mp.NODE_TYPE.PQ;
                  % above line would affect OPF (not just PF, CPF) where REF is
                  % used only as angle reference and does not require an online gen
                obj.type(obj.type == mp.NODE_TYPE.PV & ~bus_gen_status) = ...
                    mp.NODE_TYPE.PQ;
%                obj.ensure_ref_bus();   %% pick a new ref bus if one does not exist

                obj.set_vm_start(gen_dme, gbus, ig);
            end
        end

        function obj = pp_data_cnt(obj, dm, rows, out_e, mpopt, fd, pp_args)
            %% call parent
            pp_data_cnt@mp.dm_element(obj, dm, rows, out_e, mpopt, fd, pp_args);

            %% print area, zone counts
            fprintf(fd, '  %-36s%7d\n', '  Areas', ...
                length(unique(obj.tab.area)));
            fprintf(fd, '  %-36s%7d\n', '  Zones', ...
                length(unique(obj.tab.zone)));
        end

        function TorF = pp_have_section_ext(obj, mpopt, pp_args)
            TorF = true;
        end

        function obj = pp_data_ext(obj, dm, rows, out_e, mpopt, fd, pp_args)
            %% call parent
            pp_data_ext@mp.dm_element(obj, dm, rows, out_e, mpopt, fd, pp_args);

            %% print bus extremes
            [min_vm, min_vm_i] = min(obj.tab.vm);
            [max_vm, max_vm_i] = max(obj.tab.vm);
            [min_va, min_va_i] = min(obj.tab.va);
            [max_va, max_va_i] = max(obj.tab.va);

            fprintf(fd, '  %-29s %15s @ %11s %16s @ %s\n', ...
                'Bus voltage magnitude', ...
                sprintf('%7.3f p.u.', min_vm), ...
                sprintf('bus %-7d', obj.tab.uid(min_vm_i)), ...
                sprintf('%7.3f p.u.', max_vm), ...
                sprintf('bus %d', obj.tab.uid(max_vm_i)) );
            fprintf(fd, '  %-29s %15s @ %11s %16s @ %s\n', ...
                'Bus voltage angle', ...
                sprintf('%8.2f deg', min_va), ...
                sprintf('bus %-7d', obj.tab.uid(min_va_i)), ...
                sprintf('%8.2f deg', max_va), ...
                sprintf('bus %d', obj.tab.uid(max_va_i)) );
        end

        function TorF = pp_have_section_det(obj, mpopt, pp_args)
            TorF = true;
        end

        function h = pp_get_headers_det(obj, dm, out_e, mpopt, pp_args)
            h = [ pp_get_headers_det@mp.dm_element(obj, dm, out_e, mpopt, pp_args) ...
                {   '                      Voltage', ...
                    ' Bus ID   Status  Mag(pu)  Ang(deg)', ...
                    '--------  ------  -------  --------' } ];
            %%       1234567 -----1 12345.789 12345.789
        end

        function str = pp_data_row_det(obj, dm, k, out_e, mpopt, fd, pp_args)
            str = sprintf('%7d %6d %9.3f %9.3f', ...
                    obj.tab.uid(k), obj.tab.status(k), ...
                    obj.tab.vm(k), obj.tab.va(k));
        end

        function obj = set_bus_type_ref(obj, dm, idx)
            obj.tab.type(obj.on(idx)) = mp.NODE_TYPE.REF;
            obj.type(idx) = mp.NODE_TYPE.REF;
        end

        function obj = set_bus_type_pv(obj, dm, idx)
            obj.tab.type(obj.on(idx)) = mp.NODE_TYPE.PV;
            obj.type(idx) = mp.NODE_TYPE.PV;
        end

        function obj = set_bus_type_pq(obj, dm, idx)
            obj.tab.type(obj.on(idx)) = mp.NODE_TYPE.PQ;
            obj.type(idx) = mp.NODE_TYPE.PQ;
        end
    end     %% methods
end         %% classdef
