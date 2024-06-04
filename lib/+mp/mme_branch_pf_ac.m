classdef mme_branch_pf_ac < mp.mme_branch
% mp.mme_branch_pf_ac - Math model element for branch for AC power flow.
%
% Math model element class for branch elements, including transmission lines
% and transformers, for AC power flow problems.
%
% Implements updating the output data in the corresponding data model
% element for in-service branches from the math model solution.

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

            dme = obj.data_model_element(dm);
            nme = obj.network_model_element(nm);

            %% branch complex power flows
            pp = nm.get_idx('port');
            S_fr = nm.soln.gs_(pp.i1.branch(1):pp.iN.branch(1));
            S_to = nm.soln.gs_(pp.i1.branch(2):pp.iN.branch(2));

            %% branch shunt power losses
            v = nme.C' * nm.soln.v;
            vm2 = v .* conj(v);
            vm2_fr = vm2(1:length(v)/2);
            vm2_to = vm2(length(v)/2+1:end);
            psh_fr =  vm2_fr .* dme.tab.g_fr(dme.on);
            qsh_fr = -vm2_fr .* dme.tab.b_fr(dme.on);
            psh_to =  vm2_to .* dme.tab.g_to(dme.on);
            qsh_to = -vm2_to .* dme.tab.b_to(dme.on);

            %% update in the data model
            dme.tab.pl_fr(dme.on) = real(S_fr) * dm.base_mva;
            dme.tab.ql_fr(dme.on) = imag(S_fr) * dm.base_mva;
            dme.tab.pl_to(dme.on) = real(S_to) * dm.base_mva;
            dme.tab.ql_to(dme.on) = imag(S_to) * dm.base_mva;
            dme.tab.psh_fr(dme.on) = psh_fr * dm.base_mva;
            dme.tab.qsh_fr(dme.on) = qsh_fr * dm.base_mva;
            dme.tab.psh_to(dme.on) = psh_to * dm.base_mva;
            dme.tab.qsh_to(dme.on) = qsh_to * dm.base_mva;
        end
    end     %% methods
end         %% classdef
