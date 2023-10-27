classdef dmce_bus3p_mpc2 < mp.dmc_element % & mp.dmce_bus3p
% mp.dmce_bus3p_mpc2 - Data model converter element for 3-phase bus for |MATPOWER| case v2.

%   MATPOWER
%   Copyright (c) 2021-2023, Power Systems Engineering Research Center (PSERC)
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
            name = 'bus3p';
        end

        function df = data_field(obj)
            %
            df = 'bus3p';
        end

        function vmap = table_var_map(obj, dme, mpc)
            %
            vmap = table_var_map@mp.dmc_element(obj, dme, mpc);

            bsi_fcn = @(ob, mpc, spec, vn)bus_status_import(ob, mpc, spec, vn, 2);

            %% mapping for each name, default is {'col', []}
            vmap.uid{2}     = 1;
            vmap.name       = {'cell', ''};     %% empty char
            vmap.status     = {'fcn', bsi_fcn}; %% fcn w/logic for mpc.bus_types
            vmap.source_uid = {'cell', ''};     %% empty char
            vmap.base_kv{2} = 3;
            vmap.type{2}    = 2;
            vmap.vm1{2}     = 4;
            vmap.vm2{2}     = 5;
            vmap.vm3{2}     = 6;
            vmap.va1{2}     = 7;
            vmap.va2{2}     = 8;
            vmap.va3{2}     = 9;
        end

        function vals = bus_status_import(obj, mpc, spec, vn, c)
            %
            if spec.nr
                if isempty(spec.r)
                    vals = mpc.bus3p(:, 2) ~= 4;
                else
                    vals = mpc.bus3p(spec.r, 2) ~= 4;
                end
            else
                vals = [];
            end
        end
    end     %% methods
end         %% classdef
