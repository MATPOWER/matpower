classdef dme_gen_opf < mp.dme_gen & mp.dme_shared_opf
%MP.DME_GEN_OPF  MATPOWER data model class for gen data

%   MATPOWER
%   Copyright (c) 1996-2022, Power Systems Engineering Research Center (PSERC)
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
            names = horzcat( main_table_var_names@mp.dme_gen(obj), ...
                {'cost_pg', 'cost_qg', 'mu_pg_lb', 'mu_pg_ub', 'mu_qg_lb', 'mu_qg_ub'});
        end

        function vars = export_vars(obj)
            vars = horzcat( export_vars@mp.dme_gen(obj), ...
                {'vm_setpoint', 'mu_pg_lb', 'mu_pg_ub', 'mu_qg_lb', 'mu_qg_ub'} );
        end

        function TorF = have_cost(obj)
            TorF = 1;
        end

        function cost = build_cost_params(obj, dm, dc)
            base_mva = dm.base_mva;

            poly_p = obj.gen_cost_poly_params(base_mva, obj.tab.cost_pg(obj.on, :));
            if dc || ~ismember('cost_qg', obj.tab.Properties.VariableNames)
                pwl = obj.gen_cost_pwl_params(base_mva, obj.tab.cost_pg(obj.on, :), obj.n, dc);
                poly_q = [];
            else
                poly_q = obj.gen_cost_poly_params(base_mva, obj.tab.cost_qg(obj.on, :));

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
                pwl = obj.gen_cost_pwl_params(base_mva, cost, obj.n, dc);
            end

            cost = struct( ...
                    'poly_p',   poly_p, ...
                    'poly_q',   poly_q, ...
                    'pwl',      pwl ...
                );
        end

        function p = gen_cost_pwl_params(obj, base_mva, cost, ng, dc)
            ipwl = find(cost.pwl_n);    %% piece-wise linear costs
            ny = size(ipwl, 1);     %% number of piece-wise linear cost vars
            if dc
                nq = 0;     %% number of qg variables
                q1 = [];
            else
                nq = ng;    %% number of qg variables
                q1 = 1+ng;
            end

            %% from makeAy()
            ybas = 1+ng+nq;
            if ny == 0
                Ay = sparse([], [], [], 0, ybas+ny-1, 0);
                by = [];
            else
                %% if p(i),p(i+1),c(i),c(i+1) define one of the cost segments,
                %% then the corresponding constraint on pg (or qg) and Y is
                %%                                             c(i+1) - c(i)
                %%  Y   >=   c(i) + m * (pg - p(i)),      m = ---------------
                %%                                             p(i+1) - p(i)
                %%
                %% this becomes   m * pg - Y   <=   m*p(i) - c(i)

                %% form constraint matrix
                m = sum(cost.pwl_n(ipwl));  %% total number of cost points
                Ay = sparse([], [], [], m-ny, ybas+ny-1, 2*(m-ny));
                by = [];
                k = 1;
                j = 1;
                for i=ipwl'
                    ns = cost.pwl_n(i); %% # of cost points; segments = ns-1
                    p = cost.pwl_qty(i, 1:ns) / base_mva;
                    c = cost.pwl_cost(i, 1:ns);
                    m = diff(c) ./ diff(p);         %% slopes for pg (or qg)
                    if any(diff(p) == 0)
                        fprintf('mp.dme_gen/gen_cost_pwl_params: bad qty data in row %i of cost matrix\n', i);
                    end
                    b = m .* p(1:ns-1) - c(1:ns-1); %% and rhs
                    by = [by;  b'];
                    if i > ng
                        sidx = q1 + (i-ng) - 1;     %% for qg cost
                    else
                        sidx = i;                   %% for pg cost
                    end
                    Ay(k:k+ns-2, sidx) = m';
                    Ay(k:k+ns-2, ybas+j-1) = -ones(ns-1,1);
                    k = k + ns - 1;
                    j = j + 1;
                end
            end

            p = struct('n', ny, 'i', ipwl, 'A', Ay, 'b', by);
        end

        function p = gen_cost_poly_params(obj, base_mva, cost)
            ng = size(cost, 1);
            have_quad_cost = 0;
            kg = []; cg = []; Qg = [];
            i0 = find(cost.poly_n == 1);    %% constant
            i1 = find(cost.poly_n == 2);    %% linear
            i2 = find(cost.poly_n == 3);    %% quadratic
            i3 = find(cost.poly_n > 3);     %% cubic or greater
            if ~isempty(i2) || ~isempty(i1) || ~isempty(i0)
                have_quad_cost = 1;
                kg = zeros(ng, 1);
                cg = zeros(ng, 1);
                if ~isempty(i2)
                    Qg = zeros(ng, 1);
                    Qg(i2) = 2 * cost.poly_coef(i2, 3) * base_mva^2;
                    cg(i2) = cg(i2) + cost.poly_coef(i2, 2) * base_mva;
                    kg(i2) = kg(i2) + cost.poly_coef(i2, 1);
                end
                if ~isempty(i1)
                    cg(i1) = cg(i1) + cost.poly_coef(i1, 2) * base_mva;
                    kg(i1) = kg(i1) + cost.poly_coef(i1, 1);
                end
                if ~isempty(i0)
                    kg(i0) = kg(i0) + cost.poly_coef(i0, 1);
                end
            end
            p = struct( ...
                    'have_quad_cost', have_quad_cost, ...
                    'i0', i0, ...
                    'i1', i1, ...
                    'i2', i2, ...
                    'i3', i3, ...
                    'k', kg, ...
                    'c', cg, ...
                    'Q', Qg ...
                );
        end

        function maxgc = max_pwl_gencost(obj)
            maxgc = max(max(obj.tab.cost_pg.pwl_cost));
            if ismember('cost_qg', obj.tab.Properties.VariableNames) && ...
                    ~isempty(obj.tab.cost_qg.pwl_cost)
                maxgc = max(maxgc, max(max(obj.tab.cost_qg.pwl_cost)));
            end
        end

        function obj = pretty_print(obj, dm, section, out_e, mpopt, fd, pp_args)
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
            TorF = true;
        end

        function rows = pp_binding_rows_lim(obj, dm, out_e, mpopt, pp_args)
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
            if pp_args.gen.pq == 'P'
                if obj.tab.pg(k) < obj.tab.pg_lb(k) + obj.ctol || ...
                        obj.tab.mu_pg_lb(k) > obj.ptol
                    mu_lb = sprintf('%10.3f', obj.tab.mu_pg_lb(k));
                else
                    mu_lb = '      -   ';
                end
                if obj.tab.pg(k) > obj.tab.pg_ub(k) - obj.ctol || ...
                        obj.tab.mu_pg_ub(k) > obj.ptol
                    mu_ub = sprintf('%10.3f', obj.tab.mu_pg_ub(k));
                else
                    mu_ub = '      -   ';
                end

                str = sprintf('%7d %9d %10s %8.2f %8.2f %8.2f %10s', ...
                    obj.tab.uid(k), obj.tab.bus(k), mu_lb, obj.tab.pg_lb(k), ...
                    obj.tab.pg(k), obj.tab.pg_ub(k), mu_ub);
            else    % pp_args.gen.pq == 'Q'
                if obj.tab.qg(k) < obj.tab.qg_lb(k) + obj.ctol || ...
                        obj.tab.mu_qg_lb(k) > obj.ptol
                    mu_lb = sprintf('%10.3f', obj.tab.mu_qg_lb(k));
                else
                    mu_lb = '      -   ';
                end
                if obj.tab.qg(k) > obj.tab.qg_ub(k) - obj.ctol || ...
                        obj.tab.mu_qg_ub(k) > obj.ptol
                    mu_ub = sprintf('%10.3f', obj.tab.mu_qg_ub(k));
                else
                    mu_ub = '      -   ';
                end

                str = sprintf('%7d %9d %10s %8.2f %8.2f %8.2f %10s', ...
                    obj.tab.uid(k), obj.tab.bus(k), mu_lb, obj.tab.qg_lb(k), ...
                    obj.tab.qg(k), obj.tab.qg_ub(k), mu_ub);
            end
        end
    end     %% methods
end         %% classdef
