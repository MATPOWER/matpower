classdef mme_gen_pf_ac < mp.mme_gen
% mp.mme_gen_pf_ac - Math model element for generator for AC power flow.
%
% Math model element class for generator elements for AC power flow problems.
%
% Implements method for updating the output data in the corresponding data
% model element for in-service generators from the math model solution.

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

            %% generator active power
            ss = nm.get_idx('state');
            Sg = nm.soln.z(ss.i1.gen:ss.iN.gen);

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.pg(dme.on) = real(Sg) * dm.base_mva;
            dme.tab.qg(dme.on) = imag(Sg) * dm.base_mva;
        end
    end     %% methods
end         %% classdef
