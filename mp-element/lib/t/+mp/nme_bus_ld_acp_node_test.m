classdef nme_bus_ld_acp_node_test < mp.nme_bus_nld_acp_node_test

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        function name = name(obj)
            name = 'bus_ld';
        end

        function np = np(obj)
            np = 1;     %% this is a 1 port element
        end

        function obj = build_params(obj, nm, dm)
            build_params@mp.nme_bus_nld_acp_node_test(obj, nm, dm); %% call parent

            dme = obj.data_model_element(dm);

            %% bus shunts
            Ysh = dme.gs + 1j * dme.bs;             %% shunt admittances
            obj.Y = sparse(1:obj.nk, 1:obj.nk, Ysh, obj.nk, obj.nk);

            %% constant power loads
            obj.s = dme.pd + 1j * dme.qd;           %% complex power demand
        end
    end     %% methods
end         %% classdef
