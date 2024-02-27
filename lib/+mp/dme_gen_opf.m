classdef dme_gen_opf < mp.dme_gen & mp.dme_shared_opf
% mp.dme_gen_opf - Data model element for generator for OPF.
%
% To parent class :class:`mp.dme_gen`, adds costs, shadow prices on active
% and reactive generation limits, and pretty-printing for **lim** sections.
%
% Adds the following columns in the main data table, found in the
% :attr:`tab` property:
%
%   ============  ===============  =====================================
%   Name          Type             Description
%   ============  ===============  =====================================
%   ``cost_pg``   *mp.cost_table*  active power cost *(u/MW)* [#]_
%   ``cost_qg``   *mp.cost_table*  reactive power cost *(u/MVAr)* [1]_
%   ``mu_pg_lb``  *double*         shadow price on active power output lower bound *(u/MW)* [1]_
%   ``mu_pg_ub``  *double*         shadow price on active power output upper bound *(u/MW)* [1]_
%   ``mu_qg_lb``  *double*         shadow price on reactive power output lower bound *(u/MVAr)* [1]_
%   ``mu_qg_ub``  *double*         shadow price on reactive power output upper bound *(u/MVAr)* [1]_
%   ============  ===============  =====================================
%
% .. [#] Here *u* denotes the units of the objective function, e.g. USD.
%
% The cost tables ``cost_pg`` and ``cost_qg`` are defined as tables with the
% following columns:
%
% See also mp.cost_table.

%   MATPOWER
%   Copyright (c) 1996-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function names = main_table_var_names(obj)
            %
            names = horzcat( main_table_var_names@mp.dme_gen(obj), ...
                {'cost_pg', 'cost_qg', 'mu_pg_lb', 'mu_pg_ub', 'mu_qg_lb', 'mu_qg_ub'});
        end

        function vars = export_vars(obj)
            %
            vars = horzcat( export_vars@mp.dme_gen(obj), ...
                {'vm_setpoint', 'mu_pg_lb', 'mu_pg_ub', 'mu_qg_lb', 'mu_qg_ub'} );
        end

        function s = export_vars_offline_val(obj)
            %

            s = export_vars_offline_val@mp.dme_gen(obj);    %% call parent
            s.mu_pg_lb = 0;
            s.mu_pg_ub = 0;
            s.mu_qg_lb = 0;
            s.mu_qg_ub = 0;
        end

        function TorF = have_cost(obj)
            %
            TorF = 1;
        end

        function cost = build_cost_params(obj, dm, dc)
            %
            base_mva = dm.base_mva;

            poly_p = poly_params(obj.tab.cost_pg, obj.on, base_mva);
            if dc || ~ismember('cost_qg', obj.tab.Properties.VariableNames)
                pwl = pwl_params(obj.tab.cost_pg, obj.on, base_mva, obj.n, dc);
                poly_q = [];
            else
                poly_q = poly_params(obj.tab.cost_qg, obj.on, base_mva);

                %% expand cost params as needed
                polyNp = max(obj.tab.cost_pg.poly_n);
                polyNq = max(obj.tab.cost_qg.poly_n);
                pwlNp  = max(obj.tab.cost_pg.pwl_n);
                pwlNq  = max(obj.tab.cost_qg.pwl_n);
                polyN = max(polyNp, polyNq);
                pwlN = max(pwlNp, pwlNq);
                if size(obj.tab.cost_pg.poly_coef, 2) < polyN
                    obj.tab.cost_pg.poly_coef(end, polyN) = 0;
                end
                if size(obj.tab.cost_qg.poly_coef, 2) < polyN
                    obj.tab.cost_qg.poly_coef(end, polyN) = 0;
                end
                if size(obj.tab.cost_pg.pwl_cost, 2) < pwlN
                    obj.tab.cost_pg.pwl_qty(end, pwlN) = 0;
                    obj.tab.cost_pg.pwl_cost(end, pwlN) = 0;
                end
                if size(obj.tab.cost_qg.pwl_cost, 2) < pwlN
                    obj.tab.cost_qg.pwl_qty(end, pwlN) = 0;
                    obj.tab.cost_qg.pwl_cost(end, pwlN) = 0;
                end

                %% stack pg & qg cost params first
                cost = [obj.tab.cost_pg; obj.tab.cost_qg];
                pwl = pwl_params(cost, [], base_mva, obj.n, dc);
            end

            cost = struct( ...
                    'poly_p',   poly_p, ...
                    'poly_q',   poly_q, ...
                    'pwl',      pwl ...
                );
        end

        function maxgc = max_pwl_gencost(obj)
            %
            maxgc = max_pwl_cost(obj.tab.cost_pg);
            if ismember('cost_qg', obj.tab.Properties.VariableNames) && ...
                    ~isempty(obj.tab.cost_qg.pwl_cost)
                maxgc = max(maxgc, max_pwl_cost(obj.tab.cost_qg));
            end
        end

        function obj = pretty_print(obj, dm, section, out_e, mpopt, fd, pp_args)
            %
            switch section
                case 'lim'
                    pp_args.gen.pq = 'P';
                    pretty_print@mp.dme_gen(obj, dm, section, out_e, mpopt, fd, pp_args);

                    pp_args.gen.pq = 'Q';
                    pretty_print@mp.dme_gen(obj, dm, section, out_e, mpopt, fd, pp_args);
                otherwise
                    pretty_print@mp.dme_gen(obj, dm, section, out_e, mpopt, fd, pp_args);
            end
        end

        function TorF = pp_have_section_lim(obj, mpopt, pp_args)
            %
            TorF = true;
        end

        function rows = pp_binding_rows_lim(obj, dm, out_e, mpopt, pp_args)
            %
            if pp_args.gen.pq == 'P'
                rows = find( obj.tab.status & ( ...
                            obj.tab.pg < obj.tab.pg_lb + obj.ctol | ...
                            obj.tab.pg > obj.tab.pg_ub - obj.ctol | ...
                            obj.tab.mu_pg_lb > obj.ptol | ...
                            obj.tab.mu_pg_ub > obj.ptol ));
            else    % pp_args.gen.pq == 'Q'
                rows = find( obj.tab.status & ( ...
                            obj.tab.qg < obj.tab.qg_lb + obj.ctol | ...
                            obj.tab.qg > obj.tab.qg_ub - obj.ctol | ...
                            obj.tab.mu_qg_lb > obj.ptol | ...
                            obj.tab.mu_qg_ub > obj.ptol ));
            end
        end

        function h = pp_get_headers_lim(obj, dm, out_e, mpopt, pp_args)
            %
            if pp_args.gen.pq == 'P'
                h0 = pp_get_headers_lim@mp.dme_shared_opf(obj, dm, out_e, mpopt, pp_args);
                label = ' Active Power Limits';
                var = 'pg';
            else    % pp_args.gen.pq == 'Q'
                tmp = pp_args; tmp.gen.pq = 'P';
                rows = obj.pp_rows_lim(dm, out_e, mpopt, tmp);
                if isempty(rows)
                    h0 = pp_get_headers_lim@mp.dme_shared_opf(obj, dm, out_e, mpopt, pp_args);
                else
                    h0 = {''};
                end
                label = 'Reactive Power Limits';
                var = 'qg';
            end
            h = [ h0, ...
                {   sprintf(...
                    '                                 %s', label), ...
                    sprintf( ...
                    ' Gen ID    Bus ID     mu LB      LB       %s       UB       mu UB', var), ...
                    '--------  --------  ---------  -------  -------  -------   --------' } ];
            %%       1234567   1234567  12345.789 12345.78 12345.78 12345.78  12345.789
        end

        function str = pp_data_row_lim(obj, dm, k, out_e, mpopt, fd, pp_args)
            %
            if pp_args.gen.pq == 'P'
                if obj.tab.status(k) && ...
                        (obj.tab.pg(k) < obj.tab.pg_lb(k) + obj.ctol || ...
                         obj.tab.mu_pg_lb(k) > obj.ptol)
                    mu_lb = sprintf('%10.3f', obj.tab.mu_pg_lb(k));
                else
                    mu_lb = '      -   ';
                end
                if obj.tab.status(k) && ...
                        (obj.tab.pg(k) > obj.tab.pg_ub(k) - obj.ctol || ...
                         obj.tab.mu_pg_ub(k) > obj.ptol)
                    mu_ub = sprintf('%10.3f', obj.tab.mu_pg_ub(k));
                else
                    mu_ub = '      -   ';
                end
                if obj.tab.status(k)
                    pg = sprintf('%8.2f', obj.tab.pg(k));
                else
                    pg = '     -  ';
                end

                str = sprintf('%7d %9d %10s %8.2f %8s %8.2f %10s', ...
                    obj.tab.uid(k), obj.tab.bus(k), mu_lb, obj.tab.pg_lb(k), ...
                    pg, obj.tab.pg_ub(k), mu_ub);
            else    % pp_args.gen.pq == 'Q'
                if obj.tab.status(k) && ...
                        (obj.tab.qg(k) < obj.tab.qg_lb(k) + obj.ctol || ...
                         obj.tab.mu_qg_lb(k) > obj.ptol)
                    mu_lb = sprintf('%10.3f', obj.tab.mu_qg_lb(k));
                else
                    mu_lb = '      -   ';
                end
                if obj.tab.status(k) && ...
                        (obj.tab.qg(k) > obj.tab.qg_ub(k) - obj.ctol || ...
                         obj.tab.mu_qg_ub(k) > obj.ptol)
                    mu_ub = sprintf('%10.3f', obj.tab.mu_qg_ub(k));
                else
                    mu_ub = '      -   ';
                end
                if obj.tab.status(k)
                    qg = sprintf('%8.2f', obj.tab.qg(k));
                else
                    qg = '     -  ';
                end

                str = sprintf('%7d %9d %10s %8.2f %8s %8.2f %10s', ...
                    obj.tab.uid(k), obj.tab.bus(k), mu_lb, obj.tab.qg_lb(k), ...
                    qg, obj.tab.qg_ub(k), mu_ub);
            end
        end
    end     %% methods
end         %% classdef
