classdef dmce_reserve_gen_mpc2 < mp.dmc_element
% mp.dmce_reserve_gen_mpc2 - Data model converter element for reserve generator for |MATPOWER| case v2.

%   MATPOWER
%   Copyright (c) 2022-2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function name = name(obj)
            %
            name = 'reserve_gen';
        end

        function df = data_field(obj)
            %
            df = 'reserves';
        end

        function s = data_subs(obj)
            %
            s = struct('type', {'.', '.'}, 'subs', {obj.data_field(), 'cost'});
        end

        function [nr, nc, r] = get_import_size(obj, mpc)
            %
            if isfield(mpc, 'reserves')
                tab = mpc.reserves.zones;
            else
                tab = [];
            end
            if isempty(tab)
                nr = 0;
                nc = 0;
                r = [];
            else
                res = mpc.reserves;
                nrz = size(res.req, 1);     %% number of reserve zones
                if nrz > 1
                    rgens = any(res.zones); %% mask of gens available to provide reserves
                else
                    rgens = res.zones;
                end
                r = find(rgens)';   %% indices of gens available to provide reserves
                nr = size(r, 1);
                nc = 1;
            end
        end

        function [nr, nc, r] = get_export_size(obj, dme)
            %
            [nr, nc] = size(dme.tab);   %% use size of default table
            r = dme.tab.gen;            %% rows in gen matrix
        end

        function vmap = table_var_map(obj, dme, mpc)
            %
            vmap = table_var_map@mp.dmc_element(obj, dme, mpc);

            import_cost_fcn = @(a, b, c, d)import_cost(a, b, c, d);
            import_qty_fcn  = @(a, b, c, d)import_qty(a, b, c, d);
            import_ramp_fcn = @(a, b, c, d)import_ramp(a, b, c, d);

            %% mapping for each name, default is {'col', []}
            vmap.uid        = {'IDs'};      %% consecutive IDs, starting at 1
            vmap.name       = {'cell', ''}; %% empty char
            vmap.status     = {'num', 1};   %% ones
            vmap.source_uid = {'cell', ''}; %% empty char
            vmap.gen        = {'r'};        %% row index in mpc.gen
            vmap.cost       = {'fcn', import_cost_fcn};
            vmap.qty        = {'fcn', import_qty_fcn};
            vmap.ramp10     = {'fcn', import_ramp_fcn};
            vmap.r          = {'num', 0};   %% zeros
            vmap.r_lb       = {'num', 0};   %% zeros
            vmap.r_ub       = {'num', 0};   %% zeros
            vmap.total_cost = {'num', 0};   %% zeros
            vmap.prc        = {'num', 0};   %% zeros
            vmap.mu_lb      = {'num', 0};   %% zeros
            vmap.mu_ub      = {'num', 0};   %% zeros
            vmap.mu_pg_ub   = {'num', 0};   %% zeros
        end

        function val = import_cost(obj, mpc, spec, vn)
            %
            val = mpc.reserves.cost;
            if length(val) > spec.nr
                val = val(spec.r);
            end
        end

        function val = import_qty(obj, mpc, spec, vn)
            %
            if isfield(mpc.reserves, 'qty')
                val = mpc.reserves.qty;
                if length(val) > spec.nr
                    val = val(spec.r);
                end
            else
                val = Inf(spec.nr, 1);
            end
        end

        function val = import_ramp(obj, mpc, spec, vn)
            %

            %% define named indices into data matrices
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
                MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
                QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

            val = mpc.gen(spec.r, RAMP_10);
            val(val == 0) = Inf;
        end

        function dme = import(obj, dme, mpc, varargin)
            %

            %% check consistency of requirements and zones
            r = mpc.reserves;
            ng0 = size(mpc.gen, 1);     %% number of original gens (+ disp loads)
            nrz = size(r.req, 1);       %% number of reserve zones
            if nrz > 1
                rgens = any(r.zones);   %% mask of gens available to provide reserves
            else
                rgens = r.zones;
            end
            igr = find(rgens)'; %% indices of gens available to provide reserves
            ngr = length(igr);  %% number of gens available to provide reserves

            %% check data for consistent dimensions
            if size(r.zones, 1) ~= nrz
                error('mp.dmce_reserve_gen_mpc2.import: the number of rows in mpc.reserves.req (%d) and mpc.reserves.zones (%d) must match', nrz, size(r.zones, 1));
            end
            if size(r.cost, 1) ~= ng0 && size(r.cost, 1) ~= ngr
                error('mp.dmce_reserve_gen_mpc2.import: the number of rows in mpc.reserves.cost (%d) must equal the total number of generators (%d) or the number of generators able to provide reserves (%d)', size(r.cost, 1), ng0, ngr);
            end
            if isfield(r, 'qty') && size(r.qty, 1) ~= size(r.cost, 1)
                error('mp.dmce_reserve_gen_mpc2.import: mpc.reserves.cost (%d x 1) and mpc.reserves.qty (%d x 1) must be the same dimension', size(r.cost, 1), size(r.qty, 1));
            end

            %% call parent (to create dme.tab and import uid, ..., req, zone)
            dme = import@mp.dmc_element(obj, dme, mpc, varargin{:});
        end
    end     %% methods
end         %% classdef
