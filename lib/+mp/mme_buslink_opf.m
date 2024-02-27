classdef (Abstract) mme_buslink_opf < mp.mme_buslink
% mp.mme_buslink_opf - Math model element abstract base class for 1-to-3-phase buslink for OPF.
%
% Abstract math model element base class for 1-to-3-phase buslink elements
% for OPF problems.
%
% Implements (currently empty) method for forming an interior initial point.

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
        function x0 = interior_x0(obj, mm, nm, dm, x0)
            %
        end
    end     %% methods
end         %% classdef
