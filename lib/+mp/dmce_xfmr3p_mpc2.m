classdef dmce_xfmr3p_mpc2 < mp.dmc_element % & mp.dmce_xfmr3p
% mp.dmce_xfmr3p_mpc2 - Data model converter element for 3-phase transformer for |MATPOWER| case v2.

%   MATPOWER
%   Copyright (c) 2021-2023s, Power Systems Engineering Research Center (PSERC)
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
            name = 'xfmr3p';
        end

        function df = data_field(obj)
            %
            df = 'xfmr3p';
        end

        function vmap = table_var_map(obj, dme, mpc)
            %
            vmap = table_var_map@mp.dmc_element(obj, dme, mpc);

            %% mapping for each name, default is {'col', []}
            vmap.uid{2}     = 1;
            vmap.name       = {'cell', ''}; %% empty char
            vmap.status{2}  = 4;
            vmap.source_uid = {'cell', ''}; %% empty char
            vmap.bus_fr{2}  = 2;
            vmap.bus_to{2}  = 3;
            vmap.r{2}       = 5;
            vmap.x{2}       = 6;
            vmap.base_kva{2}= 7;
            vmap.base_kv{2} = 8;
            vmap.pl1_fr     = {'num', 0};   %% zeros
            vmap.ql1_fr     = {'num', 0};   %% zeros
            vmap.pl2_fr     = {'num', 0};   %% zeros
            vmap.ql2_fr     = {'num', 0};   %% zeros
            vmap.pl3_fr     = {'num', 0};   %% zeros
            vmap.ql3_fr     = {'num', 0};   %% zeros
            vmap.pl1_to     = {'num', 0};   %% zeros
            vmap.ql1_to     = {'num', 0};   %% zeros
            vmap.pl2_to     = {'num', 0};   %% zeros
            vmap.ql2_to     = {'num', 0};   %% zeros
            vmap.pl3_to     = {'num', 0};   %% zeros
            vmap.ql3_to     = {'num', 0};   %% zeros
        end
    end     %% methods
end         %% classdef
