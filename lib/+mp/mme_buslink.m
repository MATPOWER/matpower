classdef (Abstract) mme_buslink < mp.mm_element
% mp.mme_buslink - Math model element abstract base class for 1-to-3-phase buslink.
%
% Abstract math model element base class for 1-to-3-phase buslink elements.

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
            name = 'buslink';
        end
    end     %% methods
end         %% classdef
