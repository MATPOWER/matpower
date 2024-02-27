classdef (Abstract) nme_shunt_ac < mp.nme_shunt% & mp.form_ac
% mp.nme_shunt_ac - Network model element abstract base class for shunt for AC formulations.
%
% Builds the parameter :math:`\YY`.

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
            Ysh = dme.gs + 1j * dme.bs;             %% shunt admittances
            nsh = obj.nk;
            obj.Y = sparse(1:nsh, 1:nsh, Ysh, nsh, nsh);
        end
    end     %% methods
end         %% classdef
