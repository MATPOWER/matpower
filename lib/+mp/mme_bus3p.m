classdef mme_bus3p < mp.mm_element
% mp.mme_bus3p - Math model element for 3-phase bus.
%
% Math model element base class for 3-phase bus elements.
%
% Implements method for updating the output data in the corresponding data
% model element for in-service 3-phase buses from the math model solution.

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
            name = 'bus3p';
        end

        function obj = data_model_update_on(obj, mm, nm, dm, mpopt)
            %

            dme = obj.data_model_element(dm);
            nme = obj.network_model_element(nm);

            nn = nm.get_idx('node');

            for p = 1:nme.nn
                %% complex bus voltages
                v = nm.soln.v(nn.i1.bus3p(p):nn.iN.bus3p(p));

                %% update in the data model
                dme.tab.(sprintf('va%d', p))(dme.on) = angle(v) * 180/pi;
                dme.tab.(sprintf('vm%d', p))(dme.on) = abs(v);
            end
        end
    end     %% methods
end         %% classdef
