classdef mme_shunt3p_opf < mp.mme_shunt3p
% mp.mme_shunt3p_opf - Math model element for 3-phase shunt for OPF.
%
% Math model element class for 3-phase shunt elements for OPF problems.
%
% Implements (currently empty) method for forming an interior initial point.

%   MATPOWER
%   Copyright (c) 2022-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
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