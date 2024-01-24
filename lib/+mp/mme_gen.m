classdef (Abstract) mme_gen < mp.mm_element
% mp.mme_gen - Math model element abstract base class for generator.

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
            name = 'gen';
        end

        function obj = data_model_update(obj, mm, nm, dm, mpopt)
            %

            %% call parent
            data_model_update@mp.mm_element(obj, mm, nm, dm, mpopt);

            %% zero out solution values for offline elements
            dme = obj.data_model_element(dm);
            if ~isempty(dme.off)
                dme.tab.pg(dme.off) = 0;
                dme.tab.qg(dme.off) = 0;
            end
        end
    end     %% methods
end         %% classdef
