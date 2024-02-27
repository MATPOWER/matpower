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

            %% branch complex power flows
            pp = nm.get_idx('port');
            S_fr = nm.soln.gs_(pp.i1.branch(1):pp.iN.branch(1));
            S_to = nm.soln.gs_(pp.i1.branch(2):pp.iN.branch(2));

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.pl_fr(dme.on) = real(S_fr) * dm.base_mva;
            dme.tab.ql_fr(dme.on) = imag(S_fr) * dm.base_mva;
            dme.tab.pl_to(dme.on) = real(S_to) * dm.base_mva;
            dme.tab.ql_to(dme.on) = imag(S_to) * dm.base_mva;
        end
    end     %% methods
end         %% classdef
