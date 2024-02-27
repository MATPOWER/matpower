classdef nme_load3p < mp.nm_element & mp.form_acp
% mp.nme_load3p - Network model element for 3-phase load.
%
% Implements the network model element for 3-phase load elements, with
% 3 ports per 3-phase load.
%
% Builds the parameter :math:`\sv` and inherits from mp.form_acp.

%   MATPOWER
%   Copyright (c) 2021-2024, Power Systems Engineering Research Center (PSERC)
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
            name = 'load3p';
        end

        function np = np(obj)
            %
            np = 3;     %% this is a 3 port element
        end

        function obj = build_params(obj, nm, dm)
            %
            build_params@mp.nm_element(obj, nm, dm);    %% call parent
            dme = obj.data_model_element(dm);

            %% constant complex power demand
            pd = [dme.pd1 dme.pd2 dme.pd3];
            qd = pd .* tan(acos( [dme.pf1 dme.pf2 dme.pf3] ));

            obj.s = pd(:) + 1j * qd(:);
        end
    end     %% methods
end         %% classdef
