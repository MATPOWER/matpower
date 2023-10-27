classdef dmce_load3p_mpc2 < mp.dmc_element % & mp.dmce_load3p
% mp.dmce_load3p_mpc2 - Data model converter element for 3-phase load for |MATPOWER| case v2.

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
            name = 'load3p';
        end

        function df = data_field(obj)
            %
            df = 'load3p';
        end

        function vmap = table_var_map(obj, dme, mpc)
            %
            vmap = table_var_map@mp.dmc_element(obj, dme, mpc);

            %% mapping for each name, default is {'col', []}
            vmap.uid{2}     = 1;
            vmap.name       = {'cell', ''};     %% empty char
            vmap.status{2}  = 3;
            vmap.source_uid = {'cell', ''};     %% empty char
            vmap.bus{2}     = 2;
            vmap.pd1{2}     = 4;
            vmap.pd2{2}     = 5;
            vmap.pd3{2}     = 6;
            vmap.pf1{2}     = 7;
            vmap.pf2{2}     = 8;
            vmap.pf3{2}     = 9;
        end
    end     %% methods
end         %% classdef
