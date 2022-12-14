classdef dmce_gizmo_mpc2 < mp.dmc_element % & mp.dmce_gizmo
%MP.DMCE_GIZMO_MPC2  Data model converter for gizmo elements for MATPOWER case v2.

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
            name = 'gizmo';
        end

        function df = data_field(obj)
            df = 'gizmo';
        end

        function vmap = table_var_map(obj, dme, mpc)
            vmap = table_var_map@mp.dmc_element(obj, dme, mpc);

            %% mapping for each name, default is {'col', []}
            vmap.uid        = {'IDs'};      %% consecutive IDs, starting at 1
            vmap.name       = {'cell', ''}; %% empty char
            vmap.status     = {'num', 1};   %% ones
            vmap.source_uid = {'cell', ''}; %% empty char
            vmap.bus_1{2}   = 1;
            vmap.bus_2{2}   = 2;
            vmap.bus_3{2}   = 3;
            vmap.Y1r{2}     = 4;
            vmap.Y1i{2}     = 5;
            vmap.Y2r{2}     = 6;
            vmap.Y2i{2}     = 7;
            vmap.Lr{2}      = 8;
            vmap.Li{2}      = 9;
            vmap.Ir{2}      = 10;
            vmap.Ii{2}      = 11;
            vmap.M1r{2}     = 12;
            vmap.M1i{2}     = 13;
            vmap.M2r{2}     = 14;
            vmap.M2i{2}     = 15;
            vmap.Nr{2}      = 16;
            vmap.Ni{2}      = 17;
            vmap.Sr{2}      = 18;
            vmap.Si{2}      = 19;
            vmap.Zr1{2}     = 20;
            vmap.Zi1{2}     = 21;
            vmap.Zr2{2}     = 22;
            vmap.Zi2{2}     = 23;
        end
    end     %% methods
end         %% classdef
