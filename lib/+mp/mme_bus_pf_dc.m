classdef mme_bus_pf_dc < mp.mme_bus

%   MATPOWER
%   Copyright (c) 2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        function obj = data_model_update(obj, mm, nm, dm, mpopt)
            %% bus voltage angles
            nn = nm.get_idx('node');
            va = nm.soln.v(nn.i1.bus:nn.iN.bus);

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.va(dme.on) = va * 180/pi;
            dme.tab.vm(dme.on) = 1;
        end
    end     %% methods
end         %% classdef
