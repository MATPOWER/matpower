classdef (Abstract) mme_bus < mp.mm_element
% mp.mme_bus - Math model element abstract base class for bus.
%
% Abstract math model element base class for bus elements.

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
        function name = name(obj)
            %
            name = 'bus';
        end
    end     %% methods
end         %% classdef
