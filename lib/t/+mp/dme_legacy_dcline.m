classdef dme_legacy_dcline < mp.dm_element
% mp.dme_legacy_dcline - Data model element for legacy DC line.
%
% Implements the data element model for legacy DC line elements, with linear
% line losses.
%
%   :math:`p_\mathrm{loss} = \param{l}_0 + \param{l}_1 p_\mathrm{fr}`
%
% Adds the following columns in the main data table, found in the
% :attr:`tab` property:
%
%   ==================  =========  ========================================
%   Name                Type       Description
%   ==================  =========  ========================================
%   ``bus_fr``          *integer*  bus ID (``uid``) of "from" bus
%   ``bus_to``          *integer*  bus ID (``uid``) of "to" bus
%   ``loss0``           *double*   :math:`\param{l}_0`, constant term of loss function (MW)
%   ``loss1``           *double*   :math:`\param{l}_1`, linear coefficient of loss function (MW/MW)
%   ``vm_setpoint_fr``  *double*   per unit "from" bus voltage magnitude setpoint
%   ``vm_setpoint_to``  *double*   per unit "to" bus voltage magnitude setpoint
%   ``p_fr_lb``         *double*   lower bound on MW flow at "from" port
%   ``p_fr_ub``         *double*   upper bound on MW flow at "from" port
%   ``q_fr_lb``         *double*   lower bound on MVAr injection into "from" bus
%   ``q_fr_ub``         *double*   upper bound on MVAr injection into "from" bus
%   ``q_to_lb``         *double*   lower bound on MVAr injection into "to" bus
%   ``q_to_ub``         *double*   upper bound on MVAr injection into "to" bus
%   ``p_fr``            *double*   MW flow at "from" end ("from" --> "to")
%   ``q_fr``            *double*   MVAr injection into "from" bus
%   ``p_to``            *double*   MW flow at "to" end ("from" --> "to")
%   ``q_to``            *double*   MVAr injection into "to" bus
%   ==================  =========  ========================================

%   MATPOWER
%   Copyright (c) 2020-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        fbus    % bus index vector for "from" port (port 1) (all DC lines)
        tbus    % bus index vector for "to" port (port 2) (all DC lines)
        fbus_on % vector of "from" bus indices into online buses (in-service DC lines)
        tbus_on % vector of "to" bus indices into online buses (in-service DC lines)
        loss0   % constant term of loss function (p.u.) (in-service DC lines)
        loss1   % linear coefficient of loss function (in-service DC lines)
        p_fr_start  % initial active power (p.u.) at "from" port (in-service DC lines)
        p_to_start  % initial active power (p.u.) at "to" port (in-service DC lines)
        q_fr_start  % initial reactive power (p.u.) at "from" port (in-service DC lines)
        q_to_start  % initial reactive power (p.u.) at "to" port (in-service DC lines)
        vm_setpoint_fr  % from bus voltage magnitude setpoint (p.u.) (in-service DC lines)
        vm_setpoint_to  % to bus voltage magnitude setpoint (p.u.) (in-service DC lines)
        p_fr_lb % p.u. lower bound on active power flow at "from" port (in-service DC lines)
        p_fr_ub % p.u. upper bound on active power flow at "from" port (in-service DC lines)
        q_fr_lb % p.u. lower bound on reactive power flow at "from" port (in-service DC lines)
        q_fr_ub % p.u. upper bound on reactive power flow at "from" port (in-service DC lines)
        q_to_lb % p.u. lower bound on reactive power flow at "to" port (in-service DC lines)
        q_to_ub % p.u. upper bound on reactive power flow at "to" port (in-service DC lines)
    end     %% properties

    methods
        function name = name(obj)
            %
            name = 'legacy_dcline';
        end

        function label = label(obj)
            %
            label = 'DC Line';
        end

        function label = labels(obj)
            %
            label = 'DC Lines';
        end

        function name = cxn_type(obj)
            %
            name = 'bus';
        end

        function name = cxn_idx_prop(obj)
            %
            name = {'fbus', 'tbus'};
        end

        function names = main_table_var_names(obj)
            %
            names = horzcat( main_table_var_names@mp.dm_element(obj), ...
                {'bus_fr', 'bus_to', 'loss0', 'loss1', ...
                'vm_setpoint_fr', 'vm_setpoint_to', ...
                'p_fr_lb', 'p_fr_ub', ...
                'q_fr_lb', 'q_fr_ub', 'q_to_lb', 'q_to_ub', ...
                'p_fr', 'q_fr', 'p_to', 'q_to'} );
        end

        function vars = export_vars(obj)
            %
            vars = {'p_fr', 'q_fr', 'p_to', 'q_to'};
        end

        function s = export_vars_offline_val(obj)
            %

            s = export_vars_offline_val@mp.dm_element(obj);     %% call parent
            s.p_fr = 0;
            s.q_fr = 0;
            s.p_to = 0;
            s.q_to = 0;
        end

        function TorF = have_cost(obj)
            %
            TorF = 0;
        end

        function obj = initialize(obj, dm)
            %
            initialize@mp.dm_element(obj, dm);  %% call parent

            %% get bus mapping info
            b2i = dm.elements.bus.ID2i;     %% bus num to idx mapping

            %% set bus index vectors for port connectivity
            obj.fbus = b2i(obj.tab.bus_fr);
            obj.tbus = b2i(obj.tab.bus_to);
        end

        function obj = update_status(obj, dm)
            %

            %% get bus status info
            bus_dme = dm.elements.bus;
            bs = bus_dme.tab.status;    %% bus status

            %% update status of branches connected to isolated/offline buses
            obj.tab.status = obj.tab.status & bs(obj.fbus) & ...
                                              bs(obj.tbus);

            %% call parent to fill in on/off
            update_status@mp.dm_element(obj, dm);

            %% for all online DC lines ...
            %% ... set terminal buses (except ref) to PV type
            idx = [obj.fbus(obj.on); obj.tbus(obj.on)]; %% all terminal buses
            idx(bus_dme.type(idx) == mp.NODE_TYPE.REF) = [];    %% except ref
            bus_dme.set_bus_type_pv(dm, idx);

            %% set bus_dme.vm_control
            obj.fbus_on = bus_dme.i2on(obj.fbus(obj.on));
            obj.tbus_on = bus_dme.i2on(obj.tbus(obj.on));
            bus_dme.vm_control(obj.fbus_on) = 1;
            bus_dme.vm_control(obj.tbus_on) = 1;
        end

        function obj = apply_vm_setpoints(obj, dm)
            %

            % set starting bus voltage, if bus is voltage-controlled
            bus_dme = dm.elements.bus;
            i_fr = find(bus_dme.vm_control(obj.fbus_on));
            i_to = find(bus_dme.vm_control(obj.tbus_on));
            bus_dme.vm_start(obj.fbus_on(i_fr)) = obj.vm_setpoint_fr(i_fr);
            bus_dme.vm_start(obj.tbus_on(i_to)) = obj.vm_setpoint_to(i_to);
        end

        function obj = build_params(obj, dm)
            %
            obj.loss0 = obj.tab.loss0(obj.on) / dm.base_mva;
            obj.loss1 = obj.tab.loss1(obj.on);
            obj.p_fr_start = obj.tab.p_fr(obj.on) / dm.base_mva;
            obj.p_to_start = (obj.loss1 - 1) .* obj.p_fr_start + obj.loss0;
            obj.q_fr_start = -obj.tab.q_fr(obj.on) / dm.base_mva;
            obj.q_to_start = -obj.tab.q_to(obj.on) / dm.base_mva;
            obj.vm_setpoint_fr = obj.tab.vm_setpoint_fr(obj.on);
            obj.vm_setpoint_to = obj.tab.vm_setpoint_to(obj.on);
            obj.p_fr_lb = obj.tab.p_fr_lb(obj.on) / dm.base_mva;
            obj.p_fr_ub = obj.tab.p_fr_ub(obj.on) / dm.base_mva;
            obj.q_fr_lb = obj.tab.q_fr_lb(obj.on) / dm.base_mva;
            obj.q_fr_ub = obj.tab.q_fr_ub(obj.on) / dm.base_mva;
            obj.q_to_lb = obj.tab.q_to_lb(obj.on) / dm.base_mva;
            obj.q_to_ub = obj.tab.q_to_ub(obj.on) / dm.base_mva;

            obj.apply_vm_setpoints(dm);
        end

        function TorF = pp_have_section_sum(obj, mpopt, pp_args)
            %
            TorF = true;
        end

        function obj = pp_data_sum(obj, dm, rows, out_e, mpopt, fd, pp_args)
            %

            %% call parent
            pp_data_sum@mp.dm_element(obj, dm, rows, out_e, mpopt, fd, pp_args);

            %% print DC line summary
            fprintf(fd, '  %-29s  %12.2f MW', 'Total DC line losses', ...
                sum(obj.tab.p_fr(obj.on)) - sum(obj.tab.p_to(obj.on)) );
            if mpopt.model(1) ~= 'D'    %% AC model
                fprintf(fd, ' %12.2f MVAr', ...
                    sum(obj.tab.q_fr(obj.on)) + sum(obj.tab.q_to(obj.on)) );
            end
            fprintf(fd, '\n');
        end

        function h = pp_get_headers_det(obj, dm, out_e, mpopt, pp_args)
            %
            h = [ pp_get_headers_det@mp.dm_element(obj, dm, out_e, mpopt, pp_args) ...
                {   ' DC Line    From       To              Power Flow (MW)      Loss    Reactive Inj (MVAr)', ...
                    '   ID      Bus ID    Bus ID   Status    From       To       (MW)      From       To', ...
                    '--------  --------  --------  ------  --------  --------  --------  --------  --------' } ];
            %%       1234567 123456789 123456789 -----1 1234567.90 123456.89 123456.89 123456.89 123456.89
        end

        function TorF = pp_have_section_det(obj, mpopt, pp_args)
            %
            TorF = true;
        end

        function str = pp_data_row_det(obj, dm, k, out_e, mpopt, fd, pp_args)
            %
            str = sprintf('%7d %9d %9d %6d %10.2f %9.2f %9.2f %9.2f %9.2f', ...
                obj.tab.uid(k), obj.tab.bus_fr(k), obj.tab.bus_to(k), ...
                obj.tab.status(k), ...
                obj.tab.p_fr(k), obj.tab.p_to(k), ...
                obj.tab.p_fr(k) - obj.tab.p_to(k), ...
                obj.tab.q_fr(k), obj.tab.q_to(k) );
        end
    end     %% methods
end         %% classdef
