classdef mme_legacy_dcline_pf_dc < mp.mme_legacy_dcline
% mp.mme_legacy_dcline_pf_dc - Math model element for legacy DC line for DC power flow.
%
% Math model element class for legacy DC line elements for DC power flow
% problems.
%
% Implements method for updating the output data in the corresponding data
% model element for in-service DC lines from the math model solution.

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

            %% legacy DC line active power
            pp = nm.get_idx('port');
            p_fr = nm.soln.gp(pp.i1.legacy_dcline(1):pp.iN.legacy_dcline(1));
            p_to = nm.soln.gp(pp.i1.legacy_dcline(2):pp.iN.legacy_dcline(2));

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.p_fr(dme.on) = p_fr * dm.base_mva;
            dme.tab.q_fr(dme.on) = 0;
            dme.tab.p_to(dme.on) = -p_to * dm.base_mva;
            dme.tab.q_to(dme.on) = 0;
        end
    end     %% methods
end         %% classdef
