classdef dmce_bus_ld_mpc2_node_test < mp.dmce_bus_nld_mpc2_node_test % & mp.dmce_bus
%MP.DMCE_BUS_LD_MPC2_NODE_TEST  Data model converter for bus elements for MATPOWER case v2.

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function name = name(obj)
            name = 'bus_ld';
        end

        function [nr, nc, r] = get_import_size(obj, mpc)
            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
               VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            tab = mpc.(obj.data_field());
            r = find(tab(:, PD) | tab(:, QD) | tab(:, GS) | tab(:, BS));
            obj.bus = r;
            nr = size(r, 1);
            nc = size(tab, 2);          %% use nc of default table
        end

        function vmap = table_var_map(obj, dme, mpc)
            vmap = table_var_map@mp.dmce_bus_nld_mpc2_node_test(obj, dme, mpc);

            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
               VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            %% map arguments for each name
            vmap.pd{2}   = PD;
            vmap.qd{2}   = QD;
            vmap.gs{2}   = GS;
            vmap.bs{2}   = BS;
        end
    end     %% methods
end         %% classdef
