classdef dme_gen < mp.dm_element
% mp.dme_gen - Data model element for generator.
%
% Implements the data element model for generator elements.
%
% Adds the following columns in the main data table, found in the
% :attr:`tab` property:
%
%   =====================  =========  =====================================
%   Name                   Type       Description
%   =====================  =========  =====================================
%   ``bus``                *integer*  bus ID (``uid``)
%   ``vm_setpoint``        *double*   voltage magnitude setpoint *(p.u.)*
%   ``pg_lb``              *double*   active power output lower bound *(MW)*
%   ``pg_ub``              *double*   active power output upper bound *(MW)*
%   ``qg_lb``              *double*   reactive power output lower bound *(MVAr)*
%   ``qg_ub``              *double*   reactive power output upper bound *(MVAr)*
%   ``pg``                 *double*   active power output *(MW)*
%   ``qg``                 *double*   reactive power output *(MVAr)*
%   ``startup_cost_cold``  *double*   cold startup cost *(USD)*
%   ``pc1``                *double*   lower active power output of PQ
%                                     capability curve *(MW)*
%   ``pc2``                *double*   upper active power output of PQ
%                                     capability curve *(MW)*
%   ``qc1_lb``             *double*   lower bound on reactive power output
%                                     at ``pc1`` *(MVAr)*
%   ``qc1_ub``             *double*   upper bound on reactive power output
%                                     at ``pc1`` *(MVAr)*
%   ``qc2_lb``             *double*   lower bound on reactive power output
%                                     at ``pc2`` *(MVAr)*
%   ``qc2_ub``             *double*   upper bound on reactive power output
%                                     at ``pc2`` *(MVAr)*
%   =====================  =========  =====================================

%   MATPOWER
%   Copyright (c) 1996-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus         % bus index vector (all gens)
        bus_on      % vector of indices into online buses for gens that are on
        pg_start    % initial active power (p.u.) for gens that are on
        qg_start    % initial reactive power (p.u.) for gens that are on
        vm_setpoint % generator voltage setpoint for gens that are on
        pg_lb       % active power lower bound (p.u.) for gens that are on
        pg_ub       % active power upper bound (p.u.) for gens that are on
        qg_lb       % reactive power lower bound (p.u.) for gens that are on
        qg_ub       % reactive power upper bound (p.u.) for gens that are on
    end     %% properties

    methods
        function name = name(obj)
            %
            name = 'gen';
        end

        function label = label(obj)
            %
            label = 'Generator';
        end

        function label = labels(obj)
            %
            label = 'Generators';
        end

        function name = cxn_type(obj)
            %
            name = 'bus';
        end

        function name = cxn_idx_prop(obj)
            %
            name = 'bus';
        end

        function names = main_table_var_names(obj)
            %
            names = horzcat( main_table_var_names@mp.dm_element(obj), ...
                {'bus', 'vm_setpoint', 'pg_lb', 'pg_ub', 'qg_lb', 'qg_ub', ...
                'pg', 'qg', ...
                'startup_cost_cold', ...
                'pc1', 'pc2', 'qc1_lb', 'qc1_ub', 'qc2_lb', 'qc2_ub'} );
        end

        function vars = export_vars(obj)
            %
            vars = {'pg', 'qg'};
        end

        function s = export_vars_offline_val(obj)
            %

            s = export_vars_offline_val@mp.dm_element(obj);     %% call parent
            s.pg = 0;
            s.qg = 0;
        end

        function TorF = have_cost(obj)
            %
            TorF = 0;
        end

        function obj = initialize(obj, dm)
            %
            initialize@mp.dm_element(obj, dm);    %% call parent

            %% get bus mapping info
            b2i = dm.elements.bus.ID2i;     %% bus num to idx mapping

            %% set bus index vectors for port connectivity
            obj.bus = b2i(obj.tab.bus);
        end

        function obj = update_status(obj, dm)
            %

            %% get bus status info
            bus_dme = dm.elements.bus;
            bs = bus_dme.tab.status;    %% bus status

            %% update status of gens at isolated/offline buses
            obj.tab.status = obj.tab.status & bs(obj.bus);

            %% call parent to fill in on/off
            update_status@mp.dm_element(obj, dm);

            %% set bus_dme.vm_control for any online gens at PV buses
            obj.bus_on = bus_dme.i2on(obj.bus(obj.on));
            bt = bus_dme.type(obj.bus_on);  %% bus type for online gens
            i = find(bt == mp.NODE_TYPE.REF | bt == mp.NODE_TYPE.PV);
            bus_dme.vm_control(obj.bus_on(i)) = 1;
        end

        function obj = apply_vm_setpoint(obj, dm)
            %

            % set starting bus voltage, if bus is voltage-controlled
            bus_dme = dm.elements.bus;
            i = find(bus_dme.vm_control(obj.bus_on));
            bus_dme.vm_start(obj.bus_on(i)) = obj.vm_setpoint(i);
        end

        function obj = build_params(obj, dm)
            %
            base_mva = dm.base_mva;

            gen = obj.tab;

            %% get generator parameters
            obj.pg_start = gen.pg(obj.on) / base_mva;
            obj.pg_lb = gen.pg_lb(obj.on) / base_mva;
            obj.pg_ub = gen.pg_ub(obj.on) / base_mva;
            obj.qg_start  = gen.qg(obj.on) / base_mva;
            obj.qg_lb = gen.qg_lb(obj.on) / base_mva;
            obj.qg_ub = gen.qg_ub(obj.on) / base_mva;
            obj.vm_setpoint = gen.vm_setpoint(obj.on);

            obj.apply_vm_setpoint(dm);
        end

        function [mn, mx, both] = violated_q_lims(obj, dm, mpopt)
            %

            %% [mn, mx, both] = obj.violated_q_lims(dm, mpopt)
            %%  indices of online gens with violated Q lims

            gen = obj.tab;
            on = obj.on;

            %% find gens with violated Q constraints
            mx = find( gen.qg(on) > gen.qg_ub(on) + mpopt.opf.violation );
            mn = find( gen.qg(on) < gen.qg_lb(on) - mpopt.opf.violation );
            both = union(mx', mn')';    %% transposes handle fact that
                                        %% union of scalars is a row vector

            if ~isempty(both)   %% we have some Q limit violations
                %% first check for INFEASIBILITY
                %% find available online gens at REF and PV buses
                bus_dme = dm.elements.bus;
                %% bus types for buses with online gens
                bt = bus_dme.type(bus_dme.i2on(obj.bus(obj.on)));
                remaining = find( bt == mp.NODE_TYPE.REF | ...
                                  bt == mp.NODE_TYPE.PV );

                if length(both) == length(remaining) && ...
                        all(both == remaining) && (isempty(mx) || isempty(mn))
                    %% all remaining PV/REF gens are violating AND all are
                    %% violating same limit (all violating qg_lb or all qg_ub)
                    mn = [];
                    mx = [];
                else
                    %% one at a time?
                    if mpopt.pf.enforce_q_lims == 2
                        %% fix largest violation, ignore the rest
                        [junk, k] = max([gen.qg(mx) - gen.qg_ub(mx);
                                         gen.qg_lb(mn) - gen.qg(mn)]);
                        if k > length(mx)
                            mn = mn(k-length(mx));
                            mx = [];
                        else
                            mx = mx(k);
                            mn = [];
                        end
                    end
                end
            end
        end

        function TorF = isload(obj, idx)
            %
            if nargin > 1
                TorF = obj.tab.pg_lb(idx) < 0 & obj.tab.pg_ub(idx) == 0;
            else
                TorF = obj.tab.pg_lb < 0 & obj.tab.pg_ub == 0;
            end
        end

        function TorF = pp_have_section_sum(obj, mpopt, pp_args)
            %
            TorF = true;
        end

        function obj = pp_data_sum(obj, dm, rows, out_e, mpopt, fd, pp_args)
            %

            %% call parent
            pp_data_sum@mp.dm_element(obj, dm, rows, out_e, mpopt, fd, pp_args);

            ac = mpopt.model(1) ~= 'D';     %% AC model?

            %% print generation summary
            fprintf(fd, '  %-29s %12.1f MW', 'Total generation', ...
                                            sum(obj.tab.pg(obj.on)));
            if ac
                fprintf(fd, ' %12.1f MVAr', sum(obj.tab.qg(obj.on)));
            end
            fprintf(fd, '\n');
            fprintf(fd, '  %-29s %12.1f MW', 'Total max generation capacity', ...
                                            sum(obj.tab.pg_ub));
            if ac
                fprintf(fd, ' %12.1f MVAr', sum(obj.tab.qg_ub));
            end
            fprintf(fd, '\n');
            if obj.n ~= obj.nr
                fprintf(fd, '  %-29s %12.1f MW', '  online', ...
                                                sum(obj.tab.pg_ub(obj.on)));
                if ac
                    fprintf(fd, ' %12.1f MVAr', sum(obj.tab.qg_ub(obj.on)));
                end
                fprintf(fd, '\n');
            end
            fprintf(fd, '  %-29s %12.1f MW', 'Total min generation capacity', ...
                                                sum(obj.tab.pg_lb));
            if ac
                fprintf(fd, ' %12.1f MVAr', sum(obj.tab.qg_lb));
            end
            fprintf(fd, '\n');
            if obj.n ~= obj.nr
                fprintf(fd, '  %-29s %12.1f MW', '  online', ...
                                                sum(obj.tab.pg_lb(obj.on)));
                if ac
                    fprintf(fd, ' %12.1f MVAr', sum(obj.tab.qg_lb(obj.on)));
                end
                fprintf(fd, '\n');
            end
        end

        function TorF = pp_have_section_det(obj, mpopt, pp_args)
            %
            TorF = true;
        end

        function h = pp_get_headers_det(obj, dm, out_e, mpopt, pp_args)
            %
            h = [ pp_get_headers_det@mp.dm_element(obj, dm, out_e, mpopt, pp_args) ...
                {   '                             Power Generation', ...
                    ' Gen ID    Bus ID   Status   P (MW)   Q (MVAr)', ...
                    '--------  --------  ------  --------  --------' } ];
            %%       1234567 123456789 -----1 12345678.0 123456.89
        end

        function f = pp_get_footers_det(obj, dm, out_e, mpopt, pp_args)
            %
            f = {'                            --------  --------',
                sprintf('%18s Total:%10.1f %9.1f', ...
                    '', sum(obj.tab.pg(obj.on)), sum(obj.tab.qg(obj.on)))};
        end

        function str = pp_data_row_det(obj, dm, k, out_e, mpopt, fd, pp_args)
            %
            if obj.tab.status(k) && abs(obj.tab.pg(k)) > 1e-5
                pg = sprintf('%10.1f', obj.tab.pg(k));
            else
                pg = '       -  ';
            end
            if obj.tab.status(k) && abs(obj.tab.qg(k)) > 1e-5
                qg = sprintf('%9.1f', obj.tab.qg(k));
            else
                qg = '      -  ';
            end
            str = sprintf('%7d %9d %6d %10s %9s', ...
                obj.tab.uid(k), obj.tab.bus(k), obj.tab.status(k), pg, qg);
        end
    end     %% methods
end         %% classdef
