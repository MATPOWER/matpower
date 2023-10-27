classdef dmce_reserve_zone_mpc2 < mp.dmc_element % & mp.dmce_reserves
% mp.dmce_reserve_zone_mpc2 - Data model converter element for reserve zone for |MATPOWER| case v2.

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
            name = 'reserve_zone';
        end

        function df = data_field(obj)
            %
            df = 'reserves';
        end

        function s = data_subs(obj)
            %
            s = struct('type', {'.', '.'}, 'subs', {obj.data_field(), 'zones'});
        end

        function vmap = table_var_map(obj, dme, mpc)
            %
            vmap = table_var_map@mp.dmc_element(obj, dme, mpc);

            import_req_fcn   = @(a, b, c, d)import_req(a, b, c, d);
            import_zones_fcn = @(a, b, c, d)import_zones(a, b, c, d);

            %% mapping for each name, default is {'col', []}
            vmap.uid        = {'IDs'};      %% consecutive IDs, starting at 1
            vmap.name       = {'cell', ''}; %% empty char
            vmap.status     = {'num', 1};   %% ones
            vmap.source_uid = {'cell', ''}; %% empty char
            vmap.req        = {'fcn', import_req_fcn};
            vmap.zones      = {'fcn', import_zones_fcn};
            vmap.prc        = {'num', 0};   %% zeros
        end

        function val = import_req(obj, mpc, spec, vn)
            %
            val = mpc.reserves.req;
        end

        function val = import_zones(obj, mpc, spec, vn)
            %
            val = mpc.reserves.zones;
        end
    end     %% methods
end         %% classdef
