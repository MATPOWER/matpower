classdef mme_branch_pf_dc < mp.mme_branch
% mp.mme_branch_pf_dc - Math model element for branch for DC power flow.
%
% Math model element class for branch elements, including transmission lines
% and transformers, for DC power flow problems.
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

            %% branch active power flows
            pp = nm.get_idx('port');
            pl_fr = nm.soln.gp(pp.i1.branch(1):pp.iN.branch(1));
            pl_to = nm.soln.gp(pp.i1.branch(2):pp.iN.branch(2));

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.pl_fr(dme.on) = pl_fr * dm.base_mva;
            dme.tab.ql_fr(dme.on) = 0;
            dme.tab.pl_to(dme.on) = pl_to * dm.base_mva;
            dme.tab.ql_to(dme.on) = 0;
        end
    end     %% methods
end         %% classdef
