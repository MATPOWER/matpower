classdef nme_shunt_dc < mp.nme_shunt & mp.form_dc
% mp.nme_shunt_dc - Network model element for shunt for DC formulations.
%
% Builds the parameter :math:`\pv` and inherits from mp.form_dc.

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
        function obj = build_params(obj, nm, dm)
            %
            build_params@mp.nme_shunt(obj, nm, dm); %% call parent

            dme = obj.data_model_element(dm);
            obj.p = dme.gs;                         %% shunt conductances
        end
    end     %% methods
end         %% classdef
