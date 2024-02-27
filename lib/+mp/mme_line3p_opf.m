classdef mme_line3p_opf < mp.mme_line3p
% mp.mme_line3p_opf - Math model element for 3-phase line for OPF.
%
% Math model element class for 3-phase line elements for OPF problems.
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
