classdef dme_bus3p < mp.dm_element
% mp.dme_bus3p - Data model element for 3-phase bus.
%
% Implements the data element model for 3-phase bus elements.
%
% Adds the following columns in the main data table, found in the
% :attr:`tab` property:
%
%   ===========  =========  ========================================
%   Name         Type       Description
%   ===========  =========  ========================================
%   ``type``     *integer*  bus type (1 = PQ, 2 = PV, 3 = ref, 4 = isolated)
%   ``base_kv``  *double*   base voltage *(kV)*
%   ``vm1``      *double*   phase 1 voltage magnitude *(p.u.)*
%   ``vm2``      *double*   phase 2 voltage magnitude *(p.u.)*
%   ``vm3``      *double*   phase 3 voltage magnitude *(p.u.)*
%   ``va1``      *double*   phase 1 voltage angle *(degrees)*
%   ``va2``      *double*   phase 2 voltage angle *(degrees)*
%   ``va3``      *double*   phase 3 voltage angle *(degrees)*
%   ===========  =========  ========================================

%   MATPOWER
%   Copyright (c) 2021-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        type        % node type vector for buses that are on
        vm1_start   % initial phase 1 voltage magnitudes (p.u.) for buses that are on
        vm2_start   % initial phase 2 voltage magnitudes (p.u.) for buses that are on
        vm3_start   % initial phase 3 voltage magnitudes (p.u.) for buses that are on
        va1_start   % initial phase 1 voltage angles (radians) for buses that are on
        va2_start   % initial phase 2 voltage angles (radians) for buses that are on
        va3_start   % initial phase 3 voltage angles (radians) for buses that are on
        vm_control  % true if voltage is controlled, for buses that are on
    end     %% properties

    methods
        function name = name(obj)
            %
            name = 'bus3p';
        end

        function label = label(obj)
            %
            label = '3-ph Bus';
        end

        function label = labels(obj)
            %
            label = '3-ph Buses';
        end

        function names = main_table_var_names(obj)
            %
            names = horzcat( main_table_var_names@mp.dm_element(obj), ...
                {'type', 'base_kv', 'vm1', 'vm2', 'vm3', 'va1', 'va2', 'va3'});
        end

%         function vars = export_vars(obj)
%             vars = {'vm1', 'vm2', 'vm3', 'va1', 'va2', 'va3'};
%         end

%         function s = export_vars_offline_val(obj)
%             s = export_vars_offline_val@mp.dm_element(obj);     %% call parent
%         end

        function obj = init_status(obj, dm)
            %

            %% overrides mp.dm_element/init_status()

            %% check that all buses have a valid type
            bt = obj.tab.type;
            err = find(~mp.NODE_TYPE.is_valid(bt));
            if ~isempty(err)
                error('mp.dme_bus.init_status: bus %d has an invalid type', err);
            end

            %% temporarily set bus type property with dimensions for all buses
            %% (reduced for online buses only in update_status())
            obj.type = bt;
            status = (bt ~= mp.NODE_TYPE.NONE);     %% bus status
            obj.tab.status = status;
        end

        function obj = update_status(obj, dm)
            %

            %% call parent to fill in on/off
            update_status@mp.dm_element(obj, dm);

            %% update bus type property to correspond to online buses only
            obj.type = obj.type(obj.on);

            %% initialize vm_control property to all zeros, to be set
            %% by other elements, online gens, etc.
            obj.vm_control = zeros(size(obj.on));
        end

        function obj = build_params(obj, dm)
            %

            %% initialize voltage from bus table
            bus = obj.tab;
            obj.va1_start = bus.va1(obj.on) * pi/180;
            obj.va2_start = bus.va2(obj.on) * pi/180;
            obj.va3_start = bus.va3(obj.on) * pi/180;
            obj.vm1_start = bus.vm1(obj.on);
            obj.vm2_start = bus.vm2(obj.on);
            obj.vm3_start = bus.vm3(obj.on);

            %% set PV buses without online voltage controls to PQ
            i = find(obj.type == mp.NODE_TYPE.PV & ~obj.vm_control);
            obj.type(i) =  mp.NODE_TYPE.PQ; %% direct assignment to type
                                            %% property (as opposed to use of
                                            %% set_bus_type_pv() method)
                                            %% keeps it from propagating to
                                            %% output tables
        end

        function TorF = pp_have_section_det(obj, mpopt, pp_args)
            %
            TorF = true;
        end

        function h = pp_get_headers_det(obj, dm, out_e, mpopt, pp_args)
            %
            h = [ pp_get_headers_det@mp.dm_element(obj, dm, out_e, mpopt, pp_args) ...
                {   '  3-ph            Phase A Voltage    Phase B Voltage    Phase C Voltage', ...
                    ' Bus ID   Status   (kV)     (deg)     (kV)     (deg)     (kV)     (deg)', ...
                    '--------  ------  -------  -------   -------  -------   -------  -------' } ];
            %%       1234567 -----1 12345.7890 12345.78 1234.6789 12345.78 1234.6789 12345.78
        end

        function str = pp_data_row_det(obj, dm, k, out_e, mpopt, fd, pp_args)
            %
            base_kv = obj.tab.base_kv(k) / sqrt(3);
            str = sprintf('%7d %6d %10.4f %8.2f %9.4f %8.2f %9.4f %8.2f', ...
                    obj.tab.uid(k), obj.tab.status(k), ...
                    obj.tab.vm1(k) * base_kv, obj.tab.va1(k), ...
                    obj.tab.vm2(k) * base_kv, obj.tab.va2(k), ...
                    obj.tab.vm3(k) * base_kv, obj.tab.va3(k) );
        end

%         function obj = set_bus_type_ref(obj, dm, idx)
%             obj.tab.type(obj.on(idx)) = mp.NODE_TYPE.REF;
%             obj.type(idx) = mp.NODE_TYPE.REF;
%         end
% 
%         function obj = set_bus_type_pv(obj, dm, idx)
%             obj.tab.type(obj.on(idx)) = mp.NODE_TYPE.PV;
%             obj.type(idx) = mp.NODE_TYPE.PV;
%         end
% 
%         function obj = set_bus_type_pq(obj, dm, idx)
%             obj.tab.type(obj.on(idx)) = mp.NODE_TYPE.PQ;
%             obj.type(idx) = mp.NODE_TYPE.PQ;
%         end
    end     %% methods
end         %% classdef
