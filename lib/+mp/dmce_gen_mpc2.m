classdef dmce_gen_mpc2 < mp.dmc_element % & mp.dmce_gen
% mp.dmce_gen_mpc2 - Data model converter element for generator for |MATPOWER| case v2.

%   MATPOWER
%   Copyright (c) 2021-2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        % indices of single-block piecewise linear costs, all gens
        % *(automatically converted to linear cost)*
        pwl1
    end     %% properties

    methods
        function name = name(obj)
            %
            name = 'gen';
        end

        function df = data_field(obj)
            %
            df = 'gen';
        end

        function vmap = table_var_map(obj, dme, mpc)
            %
            vmap = table_var_map@mp.dmc_element(obj, dme, mpc);

            %% define named indices into data matrices
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
                MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
                QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

            sci_fcn = @(ob, mpc, spec, vn)start_cost_import(ob, mpc, spec, vn);
            gcip_fcn = @(ob, mpc, spec, vn)gen_cost_import(ob, mpc, spec, vn, 'P');
            gciq_fcn = @(ob, mpc, spec, vn)gen_cost_import(ob, mpc, spec, vn, 'Q');
            sce_fcn = @(ob, dme, mpc, spec, vn, ridx)start_cost_export(ob, dme, mpc, spec, vn, ridx);
            gcep_fcn = @(ob, dme, mpc, spec, vn, ridx)gen_cost_export(ob, dme, mpc, spec, vn, 'P', ridx);
            gceq_fcn = @(ob, dme, mpc, spec, vn, ridx)gen_cost_export(ob, dme, mpc, spec, vn, 'Q', ridx);

            %% mapping for each name, default is {'col', []}
            vmap.uid                = {'IDs'};      %% consecutive IDs, starting at 1
            vmap.name               = {'cell', ''};     %% empty char
            vmap.status{2}          = GEN_STATUS;
            vmap.source_uid         = {'cell', ''};     %% empty char
            vmap.bus{2}             = GEN_BUS;
            vmap.vm_setpoint{2}     = VG;
            vmap.pg_lb{2}           = PMIN;
            vmap.pg_ub{2}           = PMAX;
            vmap.qg_lb{2}           = QMIN;
            vmap.qg_ub{2}           = QMAX;
            vmap.pc1{2}             = PC1;
            vmap.pc2{2}             = PC2;
            vmap.qc1_lb{2}          = QC1MIN;
            vmap.qc1_ub{2}          = QC1MAX;
            vmap.qc2_lb{2}          = QC2MIN;
            vmap.qc2_ub{2}          = QC2MAX;
            vmap.pg{2}              = PG;
            vmap.qg{2}              = QG;
            vmap.startup_cost_cold  = {'fcn', sci_fcn, sce_fcn};
            if isfield(vmap, 'cost_pg')
                vmap.cost_pg        = {'fcn', gcip_fcn, gcep_fcn};
                vmap.cost_qg        = {'fcn', gciq_fcn, gceq_fcn};
                vmap.mu_pg_lb{2}    = MU_PMIN;
                vmap.mu_pg_ub{2}    = MU_PMAX;
                vmap.mu_qg_lb{2}    = MU_QMIN;
                vmap.mu_qg_ub{2}    = MU_QMAX;
            end
        end

        function dt = default_export_data_table(obj, spec)
            %

            %% define named indices into data matrices
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
                MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
                QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

            nr = obj.default_export_data_nrows(spec);
            dt = zeros(nr, APF);
        end

        function vals = start_cost_import(obj, mpc, spec, vn)
            %

            %% define named indices into data matrices
            [PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
            if isfield(mpc, 'gencost') && spec.nr
                vals = mpc.gencost(1:spec.nr, STARTUP);
            else
                vals = zeros(spec.nr, 1);
            end
        end

        function mpc = start_cost_export(obj, dme, mpc, spec, vn, ridx)
            %

            %% define named indices into data matrices
            [PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

            if dme.have_cost()
                if isempty(ridx)
                    mpc.gencost(1:spec.nr, STARTUP) = dme.tab.startup_cost_cold;
                else
                    mpc.gencost(ridx, STARTUP) = dme.tab.startup_cost_cold(ridx);
                end
            end
        end

        function val = gen_cost_import(obj, mpc, spec, vn, p_or_q)
            %
            if isfield(mpc, 'gencost') && spec.nr
                %% define named indices into data matrices
                [PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

                %% select costs for P or Q as specfied
                [pcost, qcost] = pqcost(mpc.gencost, spec.nr);
                if p_or_q == 'P'
                    %% save indices of single-block piecewise-linear costs (P & Q)
                    pwl1 = find(mpc.gencost(:, MODEL) == PW_LINEAR & mpc.gencost(:, NCOST) == 2);
                    if ~isempty(pwl1)
                        obj.pwl1 = pwl1;
                    end

                    gc = pcost;
                else    %% p_or_q = 'Q'
                    gc = qcost;
                end

                if ~isempty(gc)
                    %% reduce order where highest order has zero coefficient
                    while 1
                        ii = find(  gc(:, MODEL) == POLYNOMIAL & ...
                                    gc(:, NCOST) > 2 & ...
                                    gc(:, COST) == 0    );
                        if isempty(ii), break; end  %% exit loop if there are none
                        gc(ii, NCOST) = gc(ii, NCOST) - 1;          %% reduce order
                        gc(ii, COST:end-1) = gc(ii, COST+1:end);    %% shift coeffs
                        gc(ii, end) = 0;
                    end

                    %% convert single-block piecewise-linear costs into linear polynomial cost
                    pwl1 = find(gc(:, MODEL) == PW_LINEAR & gc(:, NCOST) == 2);
                    if ~isempty(pwl1)
                        x0 = gc(pwl1, COST);
                        y0 = gc(pwl1, COST+1);
                        x1 = gc(pwl1, COST+2);
                        y1 = gc(pwl1, COST+3);
                        m = (y1 - y0) ./ (x1 - x0);
                        b = y0 - m .* x0;
                        gc(pwl1, MODEL) = POLYNOMIAL;
                        gc(pwl1, NCOST) = 2;
                        gc(pwl1, COST:COST+1) = [m b];
                    end

                    val = mp.dmce_gen_mpc2.gencost2cost_table(gc);
                else
                    val = [];
                end
            else
                val = [];
            end
        end

        function mpc = gen_cost_export(obj, dme, mpc, spec, vn, p_or_q, ridx)
            %

            if dme.have_cost()
                [pcost, qcost] = pqcost(mpc.gencost, spec.nr);
                if p_or_q == 'P'
                    pcost = mp.dmce_gen_mpc2.cost_table2gencost( ...
                                pcost, dme.tab.cost_pg, ridx);
                    mpc.gencost(1:spec.nr, 1:size(pcost, 2)) = pcost;
                elseif ismember('cost_qg', dme.tab.Properties.VariableNames)
                    qcost = mp.dmce_gen_mpc2.cost_table2gencost( ...
                                qcost, dme.tab.cost_qg, ridx);
                    mpc.gencost(spec.nr+1:2*spec.nr, 1:size(qcost, 2)) = qcost;
                end
            end
        end
    end     %% methods

    methods (Static)
        function tab = gencost2cost_table(gencost)
            %

            %% define named indices into data matrices
            [PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

            if isempty(gencost)
                npoly = []; p = []; npwl = []; qty = []; cst = [];
            else
                ipwl = find(gencost(:, MODEL) == PW_LINEAR);
                ipol = find(gencost(:, MODEL) == POLYNOMIAL);
                nr = size(gencost, 1);

                %% get dimensions
                npoly = zeros(nr, 1);
                npwl = zeros(nr, 1);
                if isempty(ipol)
                    maxNpoly = 1;
                else
                    polycost = gencost(ipol, :);
                    npoly(ipol) = polycost(:, NCOST);
                    maxNpoly = max(npoly(ipol));
                    minNpoly = min(npoly(ipol));
                end
                if isempty(ipwl)
                    maxNpwl = 0;
                else
                    npwl(ipwl) = gencost(ipwl, NCOST);
                    maxNpwl = max(npwl);
                end

                %% initialize cost parameters
                p = zeros(nr, maxNpoly);
                qty = zeros(nr, maxNpwl);
                cst = zeros(nr, maxNpwl);

                %% polynomial costs, form coefficient matrix where 1st column
                %% is constant term, 2nd linear, etc.
                if ~isempty(ipol)
                    for n = minNpoly:maxNpoly
                        k = find(npoly(ipol) == n); %% cost with n coefficients
                        p(ipol(k), 1:n) = polycost(k, (COST+n-1):-1:COST);
                    end
                end

                %% piecewise linear costs
                if ~isempty(ipwl)
                    m = COST-1 + 2*maxNpwl;
                    qty(ipwl, :) = gencost(ipwl, COST:2:(m-1));
                    cst(ipwl, :) = gencost(ipwl, COST+1:2:m);
                end
            end
            tab = mp.cost_table(npoly, p, npwl, qty, cst);
        end

        function gencost = cost_table2gencost(gencost0, cost, ridx)
            %

            %% define named indices into data matrices
            [PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

            if isempty(cost)
                gencost = [];
            else
                gc = gencost0;
                gc(:, MODEL) = 1*(cost.pwl_n ~= 0) + 2*(cost.poly_n~=0);
                gc(:, NCOST) = cost.pwl_n + cost.poly_n;

                %% expand for polynomial or piecewise cost curve data
                n = max([cost.poly_n; 2*cost.pwl_n]);
                gc(end, COST+n-1) = 0;

                ipwl = find(cost.pwl_n);
                if ~isempty(ipwl)
                    nc = size(cost.pwl_qty, 2);
                    gc(ipwl, COST:2:COST+2*nc-1) = cost.pwl_qty(ipwl, :);
                    gc(ipwl, COST+1:2:COST+2*nc) = cost.pwl_cost(ipwl, :);
                end
                ipoly = find(cost.poly_n);
                for k = 1:length(ipoly)
                    i = ipoly(k);
                    n = cost.poly_n(i);
                    gc(i, COST:COST+n-1) = cost.poly_coef(i, n:-1:1);
                end

                if isempty(ridx)
                    gencost = gc;
                else
                    gencost = gencost0;
                    gencost(ridx, :) = gc(ridx, :);
                end
            end
        end
    end     %% methods (Static)
end         %% classdef
