classdef dme_branch < mp.dm_element
% mp.dme_branch - Data model element for branch.
%
% Implements the data element model for branch elements, including
% transmission lines and transformers.
%
% Adds the following columns in the main data table, found in the
% :attr:`tab` property:
%
%   ===========  =========  ========================================
%   Name         Type       Description
%   ===========  =========  ========================================
%   ``bus_fr``   *integer*  bus ID (``uid``) of "from" bus
%   ``bus_to``   *integer*  bus ID (``uid``) of "to" bus
%   ``r``        *double*   per unit series resistance
%   ``x``        *double*   per unit series reactance
%   ``g_fr``     *double*   per unit shunt conductance at "from" end
%   ``b_fr``     *double*   per unit shunt susceptance at "from" end
%   ``g_to``     *double*   per unit shunt conductance at "to" end
%   ``b_to``     *double*   per unit shunt susceptance at "to" end
%   ``sm_ub_a``  *double*   long term apparent power rating (MVA)
%   ``sm_ub_b``  *double*   short term apparent power rating (MVA)
%   ``sm_ub_c``  *double*   emergency apparent power rating (MVA)
%   ``cm_ub_a``  *double*   long term current magnitude rating (MVA
%                           equivalent at 1 p.u. voltage)
%   ``cm_ub_b``  *double*   short term current magnitude rating (MVA
%                           equivalent at 1 p.u. voltage)
%   ``cm_ub_c``  *double*   emergency current magnitude rating (MVA
%                           equivalent at 1 p.u. voltage)
%   ``vad_lb``   *double*   voltage angle difference lower bound
%   ``vad_ub``   *double*   voltage angle difference upper bound
%   ``tm``       *double*   transformer off-nominal turns ratio
%   ``ta``       *double*   transformer phase-shift angle (degrees)
%   ``pl_fr``    *double*   active power injection at "from" end
%   ``ql_fr``    *double*   reactive power injection at "from" end
%   ``pl_to``    *double*   active power injection at "to" end
%   ``ql_to``    *double*   reactive power injection at "to" end
%   ===========  =========  ========================================

%   MATPOWER
%   Copyright (c) 2020-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        fbus    % bus index vector for "from" port (port 1) (all branches)
        tbus    % bus index vector for "to" port (port 2) (all branches)
        r       % series resistance (p.u.) for branches that are on
        x       % series reactance (p.u.) for branches that are on
        g_fr    % shunt conductance (p.u.) at "from" end for branches that are on
        g_to    % shunt conductance (p.u.) at "to" end for branches that are on
        b_fr    % shunt susceptance (p.u.) at "from" end for branches that are on
        b_to    % shunt susceptance (p.u.) at "to" end for branches that are on
        tm      % transformer off-nominal turns ratio for branches that are on
        ta      % xformer phase-shift angle (radians) for branches that are on
        rate_a  % long term flow limit (p.u.) for branches that are on
    end     %% properties

    methods
        function name = name(obj)
            %
            name = 'branch';
        end

        function label = label(obj)
            %
            label = 'Branch';
        end

        function label = labels(obj)
            %
            label = 'Branches';
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
                {'bus_fr', 'bus_to', 'r', 'x', 'g_fr', 'b_fr', ...
                'g_to', 'b_to', 'sm_ub_a', 'sm_ub_b', 'sm_ub_c', ...
                'cm_ub_a', 'cm_ub_b', 'cm_ub_c', 'vad_lb', 'vad_ub', ...
                'tm', 'ta', ... %% remove these when we separate out xformers
                'pl_fr', 'ql_fr', 'pl_to', 'ql_to'} );
        end

        function vars = export_vars(obj)
            %
            vars = {'pl_fr', 'ql_fr', 'pl_to', 'ql_to'};
        end

        function s = export_vars_offline_val(obj)
            %

            s = export_vars_offline_val@mp.dm_element(obj);     %% call parent
            s.pl_fr = 0;
            s.ql_fr = 0;
            s.pl_to = 0;
            s.ql_to = 0;
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
            bs = dm.elements.bus.tab.status;    %% bus status

            %% update status of branches connected to isolated/offline buses
            obj.tab.status = obj.tab.status & bs(obj.fbus) & ...
                                              bs(obj.tbus);

            %% call parent to fill in on/off
            update_status@mp.dm_element(obj, dm);
        end

        function obj = build_params(obj, dm)
            %
            obj.r  = obj.tab.r(obj.on);
            obj.x  = obj.tab.x(obj.on);
            obj.g_fr  = obj.tab.g_fr(obj.on);
            obj.b_fr  = obj.tab.b_fr(obj.on);
            obj.g_to  = obj.tab.g_to(obj.on);
            obj.b_to  = obj.tab.b_to(obj.on);
            obj.tm = obj.tab.tm(obj.on);
            obj.ta = obj.tab.ta(obj.on) * pi/180;
            obj.rate_a = obj.tab.sm_ub_a(obj.on) / dm.base_mva;
        end

        function obj = pp_data_cnt(obj, dm, rows, out_e, mpopt, fd, pp_args)
            %

            %% call parent
            pp_data_cnt@mp.dm_element(obj, dm, rows, out_e, mpopt, fd, pp_args);

            num_xf = length(find(obj.tab.tm));
            num_xf_on = length(find(obj.tab.tm(obj.on)));
            num_ln = obj.nr - num_xf;
            num_ln_on = obj.n - num_xf_on;
            ln = sprintf('%d', num_ln);
            if num_ln == num_ln_on
                ln_on = ln;
                ln_off = '-';
            else
                ln_on = sprintf('%d', num_ln_on);;
                ln_off = sprintf('%d', num_ln - num_ln_on);
            end
            xf = sprintf('%d', num_xf);
            if num_xf == num_xf_on
                xf_on = xf;
                xf_off = '-';
            else
                xf_on = sprintf('%d', num_xf_on);;
                xf_off = sprintf('%d', num_xf - num_xf_on);
            end

            %% print line, transformer counts
            fprintf(fd, '  %-20s%7s %7s %7s\n', '  Lines', ln_on, ln_off, ln);
            fprintf(fd, '  %-20s%7s %7s %7s\n', '  Transformers', xf_on, xf_off, xf);
        end

        function TorF = pp_have_section_sum(obj, mpopt, pp_args)
            %
            TorF = true;
        end

        function obj = pp_data_sum(obj, dm, rows, out_e, mpopt, fd, pp_args)
            %

            %% call parent
            pp_data_sum@mp.dm_element(obj, dm, rows, out_e, mpopt, fd, pp_args);

            %% print branch summary
            fprintf(fd, '  %-29s  %12.2f MW', 'Total branch losses', ...
                sum(obj.tab.pl_fr(obj.on)) + sum(obj.tab.pl_to(obj.on)) );
            if mpopt.model(1) ~= 'D'    %% AC model
                fprintf(fd, ' %12.2f MVAr', ...
                    sum(obj.tab.ql_fr(obj.on)) + sum(obj.tab.ql_to(obj.on)) );
            end
            fprintf(fd, '\n');
        end

        function h = pp_get_headers_det(obj, dm, out_e, mpopt, pp_args)
            %
            h = [ pp_get_headers_det@mp.dm_element(obj, dm, out_e, mpopt, pp_args) ...
                {   ' Branch     From       To             From Bus Injection   To Bus Injection', ...
                    '   ID      Bus ID    Bus ID   Status   P (MW)   Q (MVAr)   P (MW)   Q (MVAr)', ...
                    '--------  --------  --------  ------  --------  --------  --------  --------' } ];
            %%       1234567 123456789 123456789 -----1 1234567.90 123456.89 123456.89 123456.89
        end

        function TorF = pp_have_section_det(obj, mpopt, pp_args)
            %
            TorF = true;
        end

        function str = pp_data_row_det(obj, dm, k, out_e, mpopt, fd, pp_args)
            %
            str = sprintf('%7d %9d %9d %6d %10.2f %9.2f %9.2f %9.2f', ...
                obj.tab.uid(k), obj.tab.bus_fr(k), obj.tab.bus_to(k), ...
                obj.tab.status(k), ...
                obj.tab.pl_fr(k), obj.tab.ql_fr(k), ...
                obj.tab.pl_to(k), obj.tab.ql_to(k) );
        end
    end     %% methods
end         %% classdef
