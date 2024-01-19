classdef dme_bus_nld_node_test < mp.dme_bus_opf
%MP.DME_BUS_NLD_NODE_TEST  MATPOWER data model class for non-load bus data for T_NODE_TEST

%   MATPOWER
%   Copyright (c) 2020-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus    %% indices into bus matrix (all rows) for this type of bus
    end     %% properties

    methods
        function name = name(obj)
            name = 'bus_nld';
        end

        function label = label(obj)
            label = 'Bus (w/o Load)';
        end

        function label = labels(obj)
            label = 'Buses (w/o Loads)';
        end

        function val = bus_eti(obj)
            val = 1;    %% element type index, 1 => bus_nld, 2 => bus_ld
        end

        function nr = count(obj, dm)
            nr = count@mp.dme_bus(obj, dm);
            if nr
                obj.bus = obj.tab.source_uid;
            end
        end
    end     %% methods
end         %% classdef
