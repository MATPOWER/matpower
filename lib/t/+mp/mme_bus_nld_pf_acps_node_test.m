classdef mme_bus_nld_pf_acps_node_test < mp.mme_bus_pf_ac

%   MATPOWER
%   Copyright (c) 2022-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        function name = name(obj)
            name = 'bus_nld';
        end

        function obj = data_model_update_on(obj, mm, nm, dm, mpopt)
            %% complex bus voltages
            nn = nm.get_idx('node');
            V = nm.soln.v(nn.i1.(obj.name):nn.iN.(obj.name));

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.va(dme.on) = angle(V) * 180/pi;
            dme.tab.vm(dme.on) = abs(V);
        end
    end     %% methods
end         %% classdef
