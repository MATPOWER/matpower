classdef dmce_line3p_mpc2 < mp.dmc_element % & mp.dmce_line3p
%MP.DMCE_LINE3P_MPC2  Data model converter for 3-phase line elements for MATPOWER case v2.

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
            name = 'line3p';
        end

        function df = data_field(obj)
            df = 'line3p';
        end

        function vmap = table_var_map(obj, dme, mpc)
            vmap = table_var_map@mp.dmc_element(obj, dme, mpc);

            %% mapping for each name, default is {'col', []}
            vmap.uid{2}     = 1;
            vmap.name       = {'cell', ''}; %% empty char
            vmap.status{2}  = 4;
            vmap.source_uid = {'cell', ''}; %% empty char
            vmap.bus_fr{2}  = 2;
            vmap.bus_to{2}  = 3;
            vmap.lc{2}      = 5;
            vmap.len{2}     = 6;
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

        function tab = create_line_construction_table(obj, lc);
            id = lc(:, 1);
            r = lc(:, 2:7);
            x = lc(:, 8:13);
            c = lc(:, 14:19);
            table_class = mp_table_class();
            tab = table_class(id, r, x, c);
%             r11 = lc(:, 2);
%             r21 = lc(:, 3);
%             r31 = lc(:, 4);
%             r22 = lc(:, 5);
%             r32 = lc(:, 6);
%             r33 = lc(:, 7);
%             x11 = lc(:, 8);
%             x21 = lc(:, 9);
%             x31 = lc(:, 10);
%             x22 = lc(:, 11);
%             x32 = lc(:, 12);
%             x33 = lc(:, 13);
%             c11 = lc(:, 14);
%             c21 = lc(:, 15);
%             c31 = lc(:, 16);
%             c22 = lc(:, 17);
%             c32 = lc(:, 18);
%             c33 = lc(:, 19);
%             if have_feature('table')
%                 tab = table(id, r11, r21, r31, r22, r32, r33, ...
%                                 x11, x21, x31, x22, x32, x33, ...
%                                 c11, c21, c31, c22, c32, c33);
%             else
%                 tab = mp_table( id, r11, r21, r31, r22, r32, r33, ...
%                                     x11, x21, x31, x22, x32, x33, ...
%                                     c11, c21, c31, c22, c32, c33);
%             end
        end

        function dme = import(obj, dme, mpc, varargin)
            %% call parent
            dme = import@mp.dmc_element(obj, dme, mpc, varargin{:});

            if ~isempty(dme.tab)
                %% system frequency
                dme.freq = mpc.freq;

                %% import line construction table
                dme.lc_tab = create_line_construction_table(obj, mpc.lc);
            end
        end
    end     %% methods
end         %% classdef
