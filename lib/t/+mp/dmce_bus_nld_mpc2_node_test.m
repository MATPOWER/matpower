classdef dmce_bus_nld_mpc2_node_test < mp.dmce_bus_mpc2 % & mp.dmce_bus
%MP.DMCE_BUS_NLD_MPC2_NODE_TEST  Data model converter for bus elements for MATPOWER case v2.

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus
    end     %% properties%     end     %% properties

    methods
        function name = name(obj)
            name = 'bus_nld';
        end

        function [nr, nc, r] = get_import_size(obj, mpc)
            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
               VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            tab = mpc.(obj.data_field());
            r = find(~tab(:, PD) & ~tab(:, QD) & ~tab(:, GS) & ~tab(:, BS));
            obj.bus = r;
            nr = size(r, 1);
            nc = size(tab, 2);          %% use nc of default table
        end

        function [nr, nc, r] = get_export_size(obj, dme)
            [nr, nc] = size(dme.tab);   %% use size of default table
            r = dme.tab.source_uid;     %% rows in bus matrix
        end

        function vmap = table_var_map(obj, dme, mpc)
            vmap = table_var_map@mp.dmce_bus_mpc2(obj, dme, mpc);

            %% mapping for each name, default is {'col', []}
            vmap.source_uid = {'r'};        %% row index in mpc.bus
        end
    end     %% methods
end         %% classdef
