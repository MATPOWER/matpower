classdef mme_shunt3p < mp.mm_element
% mp.mme_shunt3p - Math model element for 3-phase shunt.
%
% Math model element base class for 3-phase shunt elements.
%
% Implements method for updating the output data in the corresponding data
% model element for in-service 3-phase shunts from the math model solution.

%   MATPOWER
%   Copyright (c) 2021-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function name = name(obj)
            %
            name = 'shunt3p';
        end

        function obj = data_model_update_on(obj, mm, nm, dm, mpopt)
            %

            pp = nm.get_idx('port');

            dme = obj.data_model_element(dm);
            nme = obj.network_model_element(nm);

            for p = 1:nme.np
                %% shunt complex power consumption
                S = nm.soln.gs_(pp.i1.shunt3p(p):pp.iN.shunt3p(p));

                %% update in the data model
                dme.tab.(sprintf('p%d', p))(dme.on) = real(S) * dm.base_kva;
                dme.tab.(sprintf('q%d', p))(dme.on) = imag(S) * dm.base_kva;
            end
        end
    end     %% methods
end         %% classdef
