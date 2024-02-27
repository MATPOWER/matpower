classdef (Abstract) mme_legacy_dcline < mp.mm_element
% mp.mme_legacy_dcline - Math model element abstract base class for legacy DC line.
%
% Abstract math model element base class for legacy DC line elements.

%   MATPOWER
%   Copyright (c) 2021-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        function name = name(obj)
            %
            name = 'legacy_dcline';
        end
    end     %% methods
end         %% classdef
