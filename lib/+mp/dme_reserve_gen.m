classdef dme_reserve_gen < mp.dm_element & mp.dme_shared_opf
% mp.dme_reserve_gen - Data model element for reserve generator.
%
% Implements the data element model for reserve generator elements.
%
% Adds the following columns in the main data table, found in the
% :attr:`tab` property:
%
%   ==============  =========  =====================================
%   Name            Type       Description
%   ==============  =========  =====================================
%   ``gen``         *integer*  ID (``uid``) of corresponding generator
%   ``cost``        *double*   reserve cost *(u/MW)* [#]_
%   ``qty``         *double*   available reserve quantity *(MW)*
%   ``ramp10``      *double*   10-minute ramp rate *(MW)*
%   ``r``           *double*   :math:`r`, reserve allocation *(MW)*
%   ``r_lb``        *double*   lower bound on reserve allocation *(MW)*
%   ``r_ub``        *double*   upper bound on reserve allocation *(MW)*
%   ``total_cost``  *double*   total cost of allocated reserves *(u)* [1]_
%   ``prc``         *double*   reserve price *(u/MVAr)* [1]_
%   ``mu_lb``       *double*   shadow price on :math:`r` lower bound *(u/MW)* [1]_
%   ``mu_ub``       *double*   shadow price on :math:`r` upper bound *(u/MW)* [1]_
%   ``mu_pg_ub``    *double*   shadow price on capacity constraint *(u/MW)* [1]_
%   ==============  =========  =====================================
%
% .. [#] Here *u* denotes the units of the objective function, e.g. USD.

%   MATPOWER
%   Copyright (c) 2022-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        gen     % index of online gens (for online reserve gens)
        r_ub    % upper bound on reserve qty (p.u.) for units that are on
    end     %% properties

    methods
        function name = name(obj)
            %
            name = 'reserve_gen';
        end

        function label = label(obj)
            %
            label = 'Reserve Gen';
        end

        function label = labels(obj)
            %
            label = 'Reserve Gens';
        end

        function names = main_table_var_names(obj)
            %
            names = horzcat( main_table_var_names@mp.dm_element(obj), ...
                {'gen', 'cost', 'qty', 'ramp10', 'r', 'r_lb', 'r_ub', ...
                 'total_cost', 'prc', 'mu_lb', 'mu_ub', 'mu_pg_ub'});
        end

        function vars = export_vars(obj)
            %
            vars = {'r', 'r_ub', 'total_cost', 'prc', ...
                    'mu_lb', 'mu_ub', 'mu_pg_ub'};
        end

        function s = export_vars_offline_val(obj)
            %

            s = export_vars_offline_val@mp.dm_element(obj);     %% call parent
            s.r = 0;
            s.total_cost = 0;
            s.prc = 0;
            s.mu_lb = 0;
            s.mu_ub = 0;
            s.mu_pg_ub = 0;
        end

        function obj = update_status(obj, dm)
            %

            %% get gen status info
            gs = dm.elements.gen.tab.status;    %% gen status

            %% update status of reserve gens
            obj.tab.status = obj.tab.status & gs(obj.tab.gen);

            %% call parent to fill in on/off
            update_status@mp.dm_element(obj, dm);
        end

        function obj = build_params(obj, dm)
            %
            gen_dme = dm.elements.gen;
            obj.gen = gen_dme.i2on(obj.tab.gen(obj.on));

            %% p.u. upper bound for online reserve gens
            obj.r_ub = min(obj.tab.qty(obj.on), obj.tab.ramp10(obj.on)) / dm.base_mva;
        end

        function TorF = pp_have_section_sum(obj, mpopt, pp_args)
            %
            TorF = true;
        end

        function obj = pp_data_sum(obj, dm, rows, out_e, mpopt, fd, pp_args)
            %

            %% call parent
            pp_data_sum@mp.dm_element(obj, dm, rows, out_e, mpopt, fd, pp_args);

            %% print reserve summary
            fprintf(fd, '  %-29s %12.1f MW\n', 'Total reserves', ...
                                            sum(obj.tab.r(obj.on)));
            fprintf(fd, '  %-29s %12.1f MW\n', 'Total reserve capacity', ...
                                            sum(obj.tab.r_ub));
            if obj.n ~= obj.nr
                fprintf(fd, '  %-29s %12.1f MW\n', '  online', ...
                                                sum(obj.tab.r_ub(obj.on)));
            end
        end

        function TorF = pp_have_section_det(obj, mpopt, pp_args)
            %
            TorF = true;
        end

        function h = pp_get_headers_det(obj, dm, out_e, mpopt, pp_args)
            %
            h = [ pp_get_headers_det@mp.dm_element(obj, dm, out_e, mpopt, pp_args) ...
                {   '                            Reserves    Price     Cost ', ...
                    ' Gen ID    Bus ID   Status    (MW)     ($/MW)     ($)      Included in Zones ...', ...
                    '--------  --------  ------  --------  --------  --------  ------------------------' } ];
            %%       1234567 123456789 -----1 12345678.0 123456.89 123456.89    1, 2
        end

        function str = pp_data_row_det(obj, dm, k, out_e, mpopt, fd, pp_args)
            %
            gen = dm.elements.gen.tab;
            g = obj.tab.gen(k);
            zones = dm.elements.reserve_zone.tab.zones;
            iz = find(zones(:, g));
            z = sprintf('%d', iz(1));
            for j = 2:length(iz)
                z = sprintf('%s, %d', z, iz(j));
            end
            if isempty(obj.ctol)
                obj.pp_set_tols_lim(mpopt);
            end
            if obj.tab.r(k) > obj.ctol
                r = sprintf('%10.2f', obj.tab.r(k));
            else
                r = '       -  ';
            end
            if obj.tab.total_cost(k) > obj.ctol
                c = sprintf('%9.2f', obj.tab.total_cost(k));
            else
                c = '      -  ';
            end
            str = sprintf('%7d %9d %6d %10s %9.2f %9s    %s', ...
                gen.uid(g), gen.bus(g), obj.tab.status(k), ...
                r, obj.tab.prc(k), c, z);
        end

        function TorF = pp_have_section_lim(obj, mpopt, pp_args)
            %
            TorF = true;
        end

        function rows = pp_binding_rows_lim(obj, dm, out_e, mpopt, pp_args)
            %
            rows = find( obj.tab.status & ( ...
                        obj.tab.r < obj.tab.r_lb + obj.ctol | ...
                        obj.tab.r > obj.tab.r_ub - obj.ctol | ...
                        obj.tab.mu_lb > obj.ptol | ...
                        obj.tab.mu_ub > obj.ptol | ...
                        obj.tab.mu_pg_ub > obj.ptol ));
        end

        function h = pp_get_headers_lim(obj, dm, out_e, mpopt, pp_args)
            %
            h1 = pp_get_headers_lim@mp.dme_shared_opf(obj, dm, out_e, mpopt, pp_args);
            h = [ h1 ...
                {   ' Gen ID    Bus ID     mu LB      LB        r       UB       mu UB     mu pg UB', ...
                    '--------  --------  ---------  -------  -------  -------   --------   --------' } ];
            %%       1234567   1234567  12345.789 12345.78 12345.78 12345.78  12345.789  12345.789
        end

        function str = pp_data_row_lim(obj, dm, k, out_e, mpopt, fd, pp_args)
            %
            gen = dm.elements.gen.tab;
            g = obj.tab.gen(k);
            if obj.tab.r(k) > obj.ctol
                r = sprintf('%8.2f', obj.tab.r(k));
            else
                r = '     -  ';
            end
            if obj.tab.status(k) && (obj.tab.r(k) < obj.tab.r_lb(k) + obj.ctol || ...
                    obj.tab.mu_lb(k) > obj.ptol)
                mu_lb = sprintf('%10.3f', obj.tab.mu_lb(k));
            else
                mu_lb = '      -   ';
            end
            if obj.tab.status(k) && (obj.tab.r(k) > obj.tab.r_ub(k) - obj.ctol || ...
                    obj.tab.mu_ub(k) > obj.ptol)
                mu_ub = sprintf('%10.3f', obj.tab.mu_ub(k));
            else
                mu_ub = '      -   ';
            end
            if obj.tab.status(k) && obj.tab.mu_pg_ub(k) > obj.ptol
                mu_pg_ub = sprintf('%10.3f', obj.tab.mu_pg_ub(k));
            else
                mu_pg_ub = '      -   ';
            end

            str = sprintf('%7d %9d %10s %8.2f %8s %8.2f %10s %10s', ...
                gen.uid(g), gen.bus(g), mu_lb, obj.tab.r_lb(k), ...
                r, obj.tab.r_ub(k), mu_ub, mu_pg_ub);
        end

        function f = pp_get_footers_det(obj, dm, out_e, mpopt, pp_args)
            %
            f = {'                            --------            --------',
                sprintf('%18s Total:%10.2f %9s %9.2f', ...
                    '', sum(obj.tab.r(obj.on)), '', sum(obj.tab.total_cost(obj.on)))};
        end
    end     %% methods
end         %% classdef
