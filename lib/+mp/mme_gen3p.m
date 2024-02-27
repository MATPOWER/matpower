classdef mme_gen3p < mp.mm_element
% mp.mme_gen3p - Math model element for 3-phase generator.
%
% Math model element base class for 3-phase generator elements.
%
% Implements method for updating the output data in the corresponding data
% model element for in-service 3-phase generators from the math model solution.

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
            name = 'gen3p';
        end

        function obj = data_model_update_on(obj, mm, nm, dm, mpopt)
            %

            dme = obj.data_model_element(dm);
            nme = obj.network_model_element(nm);

            ss = nm.get_idx('state');

            for p = 1:nme.nz
                %% generator active power
                sg = nm.soln.z(ss.i1.gen3p(p):ss.iN.gen3p(p));

                %% update in the data model
                dme.tab.(sprintf('pg%d', p))(dme.on) = real(sg) * dm.base_kva;
                dme.tab.(sprintf('qg%d', p))(dme.on) = imag(sg) * dm.base_kva;
            end
        end
    end     %% methods
end         %% classdef
