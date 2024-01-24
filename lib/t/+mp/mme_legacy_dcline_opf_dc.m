classdef mme_legacy_dcline_opf_dc < mp.mme_legacy_dcline_opf
% mp.mme_legacy_dcline_opf_dc - Math model element for legacy DC line for DC OPF.

%   MATPOWER
%   Copyright (c) 2021-2024, Power Systems Engineering Research Center (PSERC)
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

            dme = obj.data_model_element(dm);

            %% legacy DC line active power
            ss = nm.get_idx('state');
            p_fr = nm.soln.z(ss.i1.legacy_dcline(1):ss.iN.legacy_dcline(1));
            p_to = nm.soln.z(ss.i1.legacy_dcline(2):ss.iN.legacy_dcline(2));

            %% shadow prices on legacy DC line limits
            vv = mm.get_idx();
            lambda = mm.soln.lambda;
            mu_p_fr_lb = lambda.lower(vv.i1.Pdcf:vv.iN.Pdcf);
            mu_p_fr_ub = lambda.upper(vv.i1.Pdcf:vv.iN.Pdcf);

            %% update in the data model
            dme.tab.p_fr(dme.on) = p_fr * dm.base_mva;
            dme.tab.q_fr(dme.on) = 0;
            dme.tab.p_to(dme.on) = -p_to * dm.base_mva;
            dme.tab.q_to(dme.on) = 0;
            dme.tab.mu_p_fr_lb(dme.on) = mu_p_fr_lb / dm.base_mva;
            dme.tab.mu_p_fr_ub(dme.on) = mu_p_fr_ub / dm.base_mva;
        end
    end     %% methods
end         %% classdef
