classdef (Abstract) nme_gizmo < mp.nm_element

%   MATPOWER
%   Copyright (c) 2019-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        function name = name(obj)
            name = 'gizmo';
        end

        function np = np(obj)
            np = 3;     %% this is a 3 port element
        end

        function nz = nz(obj)
            nz = 2;     %% 2 (possibly complex) non-voltage states per element
        end
    end     %% methods
end         %% classdef
