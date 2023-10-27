classdef dmce_buslink_mpc2 < mp.dmc_element % & mp.dmce_buslink
% mp.dmce_buslink_mpc2 - Data model converter element for 1-to-3-phase buslink for MATPOWER case v2.

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
            name = 'buslink';
        end

        function df = data_field(obj)
            %
            df = 'buslink';
        end

        function vmap = table_var_map(obj, dme, mpc)
            %
            vmap = table_var_map@mp.dmc_element(obj, dme, mpc);

            %% mapping for each name, default is {'col', []}
            vmap.uid{2}     = 1;
            vmap.name       = {'cell', ''};     %% empty char
            vmap.status{2}  = 4;
            vmap.source_uid = {'cell', ''};     %% empty char
            vmap.bus{2}     = 2;
            vmap.bus3p{2}   = 3;
        end
    end     %% methods
end         %% classdef
