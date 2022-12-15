classdef (Abstract) nme_gen < mp.nm_element

%   MATPOWER
%   Copyright (c) 2019-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        function name = name(obj)
            name = 'gen';
        end

        function np = np(obj)
            np = 1;     %% this is a 1 port element
        end

        function nz = nz(obj)
            nz = 1;     %% 1 (possibly complex) non-voltage state per element
        end
    end     %% methods
end         %% classdef
