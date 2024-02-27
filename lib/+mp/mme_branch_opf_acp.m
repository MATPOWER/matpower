classdef mme_branch_opf_acp < mp.mme_branch_opf_ac
% mp.mme_branch_opf_acp - Math model element for branch for AC polar voltage OPF.
%
% Math model element class for branch elements, including transmission lines
% and transformers, for AC polar voltage OPF problems.
%
% Implements method for adding branch angle difference constraints.

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
        function obj = add_constraints(obj, mm, nm, dm, mpopt)
            %

            %% call parent
            add_constraints@mp.mme_branch_opf_ac(obj, mm, nm, dm, mpopt);

            %% branch voltage angle difference limits
            [Aang, lang, uang, iang] = obj.ang_diff_params(...
                    dm, mpopt.opf.ignore_angle_lim);
            if length(iang)
                mm.add_lin_constraint('ang', Aang, lang, uang, {'Va'});
            end
            mm.userdata.ang_diff_constrained_branch_idx = iang;
        end
    end     %% methods
end         %% classdef
