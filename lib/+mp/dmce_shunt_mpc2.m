classdef dmce_shunt_mpc2 < mp.dmc_element % & mp.dmce_shunt
% mp.dmce_shunt_mpc2 - Data model converter element for shunt for |MATPOWER| case v2.

%   MATPOWER
%   Copyright (c) 2021-2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus     % 
    end     %% properties

    methods
        function name = name(obj)
            %
            name = 'shunt';
        end

        function df = data_field(obj)
            %
            df = 'bus';
        end

        function [nr, nc, r] = get_import_size(obj, mpc)
            %

            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
               VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            if isfield(mpc, obj.data_field())
                tab = mpc.(obj.data_field());
            else
                tab = [];
            end
            if isempty(tab)
                nr = 0;
                nc = 0;
                r = [];
            else
                r = find(tab(:, GS) | tab(:, BS));
                obj.bus = r;
                nr = size(r, 1);
                nc = size(tab, 2);          %% use nc of default table
            end
        end

        function [nr, nc, r] = get_export_size(obj, dme)
            %
            [nr, nc] = size(dme.tab);   %% use size of default table
            r = dme.tab.source_uid;     %% rows in bus matrix
        end

        function vmap = table_var_map(obj, dme, mpc)
            %
            vmap = table_var_map@mp.dmc_element(obj, dme, mpc);

            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
               VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            %% mapping for each name, default is {'col', []}
            vmap.uid        = {'IDs'};      %% consecutive IDs, starting at 1
            vmap.name       = {'cell', ''}; %% empty char
            vmap.status     = {'num', 1};   %% ones
            vmap.source_uid = {'r'};        %% row index in mpc.bus
            vmap.bus{2}     = BUS_I;
            vmap.gs{2}      = GS;
            vmap.bs{2}      = BS;
            vmap.p          = {'num', 0};   %% zeros
            vmap.q          = {'num', 0};   %% zeros
        end
    end     %% methods
end         %% classdef
