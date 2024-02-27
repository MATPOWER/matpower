classdef dme_legacy_dcline_opf < mp.dme_legacy_dcline & mp.dme_shared_opf
% mp.dme_legacy_dcline_opf - Data model element for legacy DC line for OPF.
%
% To parent class :class:`mp.dme_legacy_dcline`, adds costs, shadow prices
% on active and reactive flow limits, and pretty-printing for **lim**
% sections.
%
% Adds the following columns in the main data table, found in the
% :attr:`tab` property:
%
%   ==============  ===============  ===================================
%   Name            Type             Description
%   ==============  ===============  ===================================
%   ``cost_pg``     *mp.cost_table*  cost of active power flow *(u/MW)* [#]_
%   ``mu_p_fr_lb``  *double*         shadow price on MW flow lower bound at
%                                    "from" end *(u/MW)* [1]_
%   ``mu_p_fr_ub``  *double*         shadow price on MW flow upper bound at
%                                    "from" end *(u/MW)* [1]_
%   ``mu_q_fr_lb``  *double*         shadow price on lower bound of MVAr
%                                    injection at "from" bus *(u/degree)* [1]_
%   ``mu_q_fr_ub``  *double*         shadow price on upper bound of MVAr
%                                    injection at "from" bus *(u/degree)* [1]_
%   ``mu_q_to_lb``  *double*         shadow price on lower bound of MVAr
%                                    injection at "to" bus *(u/degree)* [1]_
%   ``mu_q_to_ub``  *double*         shadow price on upper bound of MVAr
%                                    injection at "to" bus *(u/degree)* [1]_
%   ==============  ===============  ===================================
%
% .. [#] Here *u* denotes the units of the objective function, e.g. USD.

%   MATPOWER
%   Copyright (c) 2020-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        function names = main_table_var_names(obj)
            %
            names = horzcat( main_table_var_names@mp.dme_legacy_dcline(obj), ...
                {   'cost', ...
                    'mu_p_fr_lb', 'mu_p_fr_ub', ...
                    'mu_q_fr_lb', 'mu_q_fr_ub', ...
                    'mu_q_to_lb', 'mu_q_to_ub'  } );
        end

        function vars = export_vars(obj)
            %
            vars = horzcat( export_vars@mp.dme_legacy_dcline(obj), ...
                {   'vm_setpoint_fr', 'vm_setpoint_to', ...
                    'mu_p_fr_lb', 'mu_p_fr_ub', ...
                    'mu_q_fr_lb', 'mu_q_fr_ub', ...
                    'mu_q_to_lb', 'mu_q_to_ub'  } );
        end

        function s = export_vars_offline_val(obj)
            %

            s = export_vars_offline_val@mp.dme_legacy_dcline(obj);  %% call parent
            s.mu_p_fr_lb = 0;
            s.mu_p_fr_ub = 0;
            s.mu_q_fr_lb = 0;
            s.mu_q_fr_ub = 0;
            s.mu_q_to_lb = 0;
            s.mu_q_to_ub = 0;
        end

        function TorF = have_cost(obj)
            %
            TorF = 1;
        end

        function cost = build_cost_params(obj, dm)
            %
            if ismember('cost', obj.tab.Properties.VariableNames)
                poly = mp.cost_table_utils.poly_params(obj.tab.cost, obj.on, dm.base_mva);
                pwl = mp.cost_table_utils.pwl_params(obj.tab.cost, obj.on, dm.base_mva);
    
                cost = struct( ...
                        'poly', poly, ...
                        'pwl',  pwl ...
                    );
            else
                cost = struct([]);
            end
        end

        function obj = pretty_print(obj, dm, section, out_e, mpopt, fd, pp_args)
            %
            switch section
                case 'lim'
                    %% pass flows and limits to parent
                    p_fr = obj.tab.p_fr;
                    lb = obj.tab.p_fr_lb;
                    ub = obj.tab.p_fr_ub;
                    pp_args.legacy_dcline.flow = struct( 'p_fr', p_fr, ...
                                                         'lb', lb, ...
                                                         'ub', ub );
            end

            pretty_print@mp.dme_legacy_dcline(obj, dm, section, out_e, mpopt, fd, pp_args);
        end

        function TorF = pp_have_section_lim(obj, mpopt, pp_args)
            %
            TorF = true;
        end

        function rows = pp_binding_rows_lim(obj, dm, out_e, mpopt, pp_args)
            %
            flow = pp_args.legacy_dcline.flow;
            rows = find( obj.tab.status & ( ...
                        flow.p_fr < flow.lb + obj.ctol | ...
                        flow.p_fr > flow.ub - obj.ctol | ...
                        obj.tab.mu_p_fr_lb > obj.ptol | ...
                        obj.tab.mu_p_fr_ub > obj.ptol ));
        end

        function h = pp_get_headers_lim(obj, dm, out_e, mpopt, pp_args)
            %
            h = [ pp_get_headers_lim@mp.dme_shared_opf(obj, dm, out_e, mpopt, pp_args) ...
                {   ' DC Line    From       To                   Active Power Flow (MW)', ...
                    '   ID      Bus ID    Bus ID     mu LB       LB      p_fr      UB      mu UB', ...
                    '--------  --------  --------  ---------  -------  -------  -------  ---------' } ];
            %%       1234567 123456789 123456789 123456.890 12345.78 12345.78 12345.78 123456.890
        end

        function str = pp_data_row_lim(obj, dm, k, out_e, mpopt, fd, pp_args)
            %
            flow = pp_args.legacy_dcline.flow;

            if (flow.p_fr(k) < flow.lb(k) + obj.ctol || ...
                    obj.tab.mu_p_fr_lb(k) > obj.ptol)
                mu_lb = sprintf('%10.3f', obj.tab.mu_p_fr_lb(k));
            else
                mu_lb = '      -   ';
            end
            if (flow.p_fr(k) > flow.ub(k) - obj.ctol || ...
                    obj.tab.mu_p_fr_ub(k) > obj.ptol)
                mu_ub = sprintf('%10.3f', obj.tab.mu_p_fr_ub(k));
            else
                mu_ub = '      -   ';
            end

            str = sprintf('%7d %9d %9d %10s %8.2f %8.2f %8.2f %10s', ...
                obj.tab.uid(k), obj.tab.bus_fr(k), obj.tab.bus_to(k), ...
                mu_lb, flow.lb(k), flow.p_fr(k), flow.ub(k), mu_ub);
        end
    end     %% methods
end         %% classdef
