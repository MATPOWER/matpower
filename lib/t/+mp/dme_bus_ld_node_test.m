classdef dme_bus_ld_node_test < mp.dme_bus_nld_node_test
%MP.DME_BUS_LD_NODE_TEST  MATPOWER data model class for load bus data for T_NODE_TEST

%   MATPOWER
%   Copyright (c) 2020-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        pd      %% active power demand (p.u.) for loads at buses that are on
        qd      %% reactive power demand (p.u.) for loads at buses that are on
        gs      %% shunt conductance (p.u. active power demanded at
                %% V = 1.0 p.u.) for shunts at buses that are on
        bs      %% shunt susceptance (p.u. reactive power injected at
                %% V = 1.0 p.u.) for shunts at buses that are on
    end     %% properties

    methods
        function name = name(obj)
            name = 'bus_ld';
        end

        function label = label(obj)
            label = 'Bus (w/Load)';
        end

        function label = labels(obj)
            label = 'Buses (w/Loads)';
        end

        function name = cxn_type(obj)
            name = 'bus_ld';
        end

        function val = bus_eti(obj)
            val = 2;    %% element type index, 1 => bus_nld, 2 => bus_ld
        end

        function names = main_table_var_names(obj)
            names = horzcat( ...
                main_table_var_names@mp.dme_bus_nld_node_test(obj), ...
                {'pd', 'qd', 'gs', 'bs'} );
        end

        function obj = build_params(obj, dm)
            build_params@mp.dme_bus(obj, dm);   %% call parent

            obj.pd = obj.tab.pd(obj.on) / dm.base_mva;
            obj.qd = obj.tab.qd(obj.on) / dm.base_mva;
            obj.gs = obj.tab.gs(obj.on) / dm.base_mva;
            obj.bs = obj.tab.bs(obj.on) / dm.base_mva;
        end
    end     %% methods
end         %% classdef
