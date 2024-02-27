classdef mme_load_pf_ac < mp.mme_load
% mp.mme_load_pf_ac - Math model element for load for AC power flow.
%
% Math model element class for load elements for AC power flow problems.
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

            %% load complex power consumption
            pp = nm.get_idx('port');
            S = nm.soln.gs_(pp.i1.load:pp.iN.load);

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.p(dme.on) = real(S) * dm.base_mva;
            dme.tab.q(dme.on) = imag(S) * dm.base_mva;
        end
    end     %% methods
end         %% classdef
