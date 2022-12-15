classdef nme_branch_acp_nln < mp.nme_branch_acp & mp.nme_wrapper_ac_nln

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        nme_class = @mp.nme_branch_acp;
    end

    methods
        function obj = nme_branch_acp_nln()
            obj@mp.nme_branch_acp();
            obj.nme_wrapper_ac_nln_init();
        end

        function build_params(obj, nm, dm)
            build_params@mp.nme_branch_acp(obj, nm, dm);
            obj.build_nln_params(nm, dm);
        end

        function nk = count(obj, dm)
            obj.nme.count(dm);
            nk = count@mp.nme_branch_acp(obj, dm);
        end
    end     %% methods
end         %% classdef
