classdef mme_legacy_dcline_opf_ac < mp.mme_legacy_dcline_opf
% mp.mme_legacy_dcline_opf_ac - Math model element for legacy DC line for AC OPF.
%
% Math model element class for legacy DC line elements for AC OPF problems.
%
% Implements method for updating the output data in the corresponding data
% model element for in-service DC lines from the math model solution.

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
            nme = obj.network_model_element(nm);

            %% legacy DC line active power
            pp = nm.get_idx('port');
            s_fr = nm.soln.gs_(pp.i1.legacy_dcline(1):pp.iN.legacy_dcline(1));
            s_to = nm.soln.gs_(pp.i1.legacy_dcline(2):pp.iN.legacy_dcline(2));
            vm_setpoint = abs(nme.C' * nm.soln.v);

            %% shadow prices on legacy DC line limits
            vv = mm.get_idx();
            lambda = mm.soln.lambda;
            mu_p_fr_lb = lambda.lower(vv.i1.Pdcf:vv.iN.Pdcf);
            mu_p_fr_ub = lambda.upper(vv.i1.Pdcf:vv.iN.Pdcf);
            mu_q_fr_lb = lambda.lower(vv.i1.Qdcf:vv.iN.Qdcf);
            mu_q_fr_ub = lambda.upper(vv.i1.Qdcf:vv.iN.Qdcf);
            mu_q_to_lb = lambda.lower(vv.i1.Qdct:vv.iN.Qdct);
            mu_q_to_ub = lambda.upper(vv.i1.Qdct:vv.iN.Qdct);

            %% update in the data model
            dme.tab.p_fr(dme.on) = real(s_fr) * dm.base_mva;
            dme.tab.q_fr(dme.on) = -imag(s_fr) * dm.base_mva;
            dme.tab.p_to(dme.on) = -real(s_to) * dm.base_mva;
            dme.tab.q_to(dme.on) = -imag(s_to) * dm.base_mva;
            dme.tab.vm_setpoint_fr(dme.on) = vm_setpoint(1:dme.n);
            dme.tab.vm_setpoint_to(dme.on) = vm_setpoint(dme.n+1:end);
            dme.tab.mu_p_fr_lb(dme.on) = mu_p_fr_lb / dm.base_mva;
            dme.tab.mu_p_fr_ub(dme.on) = mu_p_fr_ub / dm.base_mva;
            %% because of sign swap on reactive quantities, bounds are swapped
            dme.tab.mu_q_fr_lb(dme.on) = mu_q_fr_ub / dm.base_mva;
            dme.tab.mu_q_fr_ub(dme.on) = mu_q_fr_lb / dm.base_mva;
            dme.tab.mu_q_to_lb(dme.on) = mu_q_to_ub / dm.base_mva;
            dme.tab.mu_q_to_ub(dme.on) = mu_q_to_lb / dm.base_mva;
        end
    end     %% methods
end         %% classdef
