classdef mme_load_pf_dc < mp.mme_load
% mp.mme_load_pf_dc - Math model element for load for DC power flow.
%
% Math model element class for load elements for DC power flow problems.
%
% Implements method for updating the output data in the corresponding data
% model element for in-service loads from the math model solution.

%   MATPOWER
%   Copyright (c) 2022-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function obj = data_model_update_on(obj, mm, nm, dm, mpopt)
            %

            %% load active power consumption
            pp = nm.get_idx('port');
            P = nm.soln.gp(pp.i1.load:pp.iN.load);

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.p(dme.on) = P * dm.base_mva;
            dme.tab.q(dme.on) = 0;
        end
    end     %% methods
end         %% classdef
