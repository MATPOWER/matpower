classdef (Abstract) mme_load < mp.mm_element
% mp.mme_load - Math model element abstract base class for load.
%
% Abstract math model element base class for load elements.

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
        function name = name(obj)
            %
            name = 'load';
        end
    end     %% methods
end         %% classdef
