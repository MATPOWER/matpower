classdef (Abstract) nme_legacy_dcline < mp.nm_element
% mp.nme_legacy_dcline - Network model element abstract base class for legacy DC line.

%   MATPOWER
%   Copyright (c) 2019-2024, Power Systems Engineering Research Center (PSERC)
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
            name = 'legacy_dcline';
        end

        function np = np(obj)
            %
            np = 2;     %% this is a 1 port element
        end

        function nz = nz(obj)
            %
            nz = 2;     %% 2 (possibly complex) non-voltage state per element
        end
    end     %% methods
end         %% classdef
