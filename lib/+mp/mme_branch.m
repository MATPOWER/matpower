classdef (Abstract) mme_branch < mp.mm_element
% mp.mme_branch - Math model element abstract base class for branch.

%   MATPOWER
%   Copyright (c) 2021-2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function name = name(obj)
            %
            name = 'branch';
        end

        function obj = data_model_update(obj, mm, nm, dm, mpopt)
            %

            %% call parent
            data_model_update@mp.mm_element(obj, mm, nm, dm, mpopt);

            %% zero out solution values for offline elements
            dme = obj.data_model_element(dm);
            if ~isempty(dme.off)
                dme.tab.pl_fr(dme.off) = 0;
                dme.tab.ql_fr(dme.off) = 0;
                dme.tab.pl_to(dme.off) = 0;
                dme.tab.ql_to(dme.off) = 0;
            end
        end
    end     %% methods
end         %% classdef
