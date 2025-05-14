classdef nme_shunt3p < mp.nm_element & mp.form_acp
% mp.nme_shunt3p - Network model element for 3-phase shunt.
%
% Implements the network model element for 3-phase shunt elements, with
% 3 ports per 3-phase shunt.
%
% Builds the parameter :math:`\YY` and inherits from mp.form_acp.

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

        function np = np(obj)
            %
            np = 3;     %% this is a 3 port element
        end

        function obj = build_params(obj, nm, dm)
            %
            build_params@mp.nm_element(obj, nm, dm);    %% call parent
            dme = obj.data_model_element(dm);

            %% constant complex power demand
            gs = [dme.gs1 dme.gs2 dme.gs3];
            bs = [dme.bs1 dme.bs2 dme.bs3];
            Ysh = gs + 1j * bs;                         %% shunt admittances
            nsh = 3 * obj.nk;

            obj.Y = sparse(1:nsh, 1:nsh, Ysh(:), nsh, nsh);
        end
    end     %% methods
end         %% classdef
