classdef mme_branch_pf_dc < mp.mme_branch

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
            %% branch active power flows
            pp = nm.get_idx('port');
            pl_fr = nm.soln.gp(pp.i1.branch(1):pp.iN.branch(1)) * dm.base_mva;
            pl_to = nm.soln.gp(pp.i1.branch(2):pp.iN.branch(2)) * dm.base_mva;

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.pl_fr(dme.on) = pl_fr;
            dme.tab.ql_fr(dme.on) = 0;
            dme.tab.pl_to(dme.on) = pl_to;
            dme.tab.ql_to(dme.on) = 0;
        end
    end     %% methods
end         %% classdef
