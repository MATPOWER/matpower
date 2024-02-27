classdef (Abstract) nme_branch < mp.nm_element
% mp.nme_branch - Network model element abstract base class for branch.
%
% Implements the network model element for branch elements, including
% transmission lines and transformers, with 2 ports per branch.

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
            name = 'branch';
        end

        function np = np(obj)
            %
            np = 2;     %% this is a 2 port element
        end
    end     %% methods
end         %% classdef
