classdef dmce_shunt3p_mpc2 < mp.dmc_element
% mp.dmce_shunt3p_mpc2 - Data model converter element for 3-phase shunt for |MATPOWER| case v2.

%   MATPOWER
%   Copyright (c) 2021-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%   and Ray Zimmerman, PSERC Cornell
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
            name = 'shunt3p';
        end

        function df = data_field(obj)
            %
            df = 'shunt3p';
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
            vmap.gs1{2}     = 4;
            vmap.gs2{2}     = 5;
            vmap.gs3{2}     = 6;
            vmap.bs1{2}     = 7;
            vmap.bs2{2}     = 8;
            vmap.bs3{2}     = 9;
            vmap.p1         = {'num', 0};   %% zeros
            vmap.p2         = {'num', 0};   %% zeros
            vmap.p3         = {'num', 0};   %% zeros
            vmap.q1         = {'num', 0};   %% zeros
            vmap.q2         = {'num', 0};   %% zeros
            vmap.q3         = {'num', 0};   %% zeros
        end
    end     %% methods
end         %% classdef
