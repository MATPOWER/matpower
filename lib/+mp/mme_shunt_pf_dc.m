classdef mme_shunt_pf_dc < mp.mme_shunt

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
            %% shunt complex power consumption
            pp = nm.get_idx('port');
            P = nm.soln.gp(pp.i1.shunt:pp.iN.shunt) * dm.base_mva;

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.p(dme.on) = P;
            dme.tab.q(dme.on) = 0;
        end
    end     %% methods
end         %% classdef
