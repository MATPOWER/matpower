classdef mme_branch_pf_ac < mp.mme_branch

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
            %% branch complex power flows
            pp = nm.get_idx('port');
            S_fr = nm.soln.gs_(pp.i1.branch(1):pp.iN.branch(1)) * dm.base_mva;
            S_to = nm.soln.gs_(pp.i1.branch(2):pp.iN.branch(2)) * dm.base_mva;

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.pl_fr(dme.on) = real(S_fr);
            dme.tab.ql_fr(dme.on) = imag(S_fr);
            dme.tab.pl_to(dme.on) = real(S_to);
            dme.tab.ql_to(dme.on) = imag(S_to);
        end
    end     %% methods
end         %% classdef
