classdef dmce_line3p_mpc2 < mp.dmc_element % & mp.dmce_line3p
% mp.dmce_line3p_mpc2 - Data model converter element for 3-phase line for |MATPOWER| case v2.

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
            name = 'line3p';
        end

        function df = data_field(obj)
            %
            df = 'line3p';
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

        function tab = create_line_construction_table(obj, dme, lc)
            %
            id = lc(:, 1);
            r = lc(:, 2:7);
            x = lc(:, 8:13);
            c = lc(:, 14:19);
            tab = dme.create_line_construction_table(id, r, x, c);
        end

        function dme = import(obj, dme, mpc, varargin)
            %

            %% call parent
            dme = import@mp.dmc_element(obj, dme, mpc, varargin{:});

            if ~isempty(dme.tab)
                %% system frequency
                dme.freq = mpc.freq;

                %% import line construction table
                dme.lc_tab = create_line_construction_table(obj, dme, mpc.lc);
            end
        end
    end     %% methods
end         %% classdef
