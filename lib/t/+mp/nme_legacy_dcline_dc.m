classdef nme_legacy_dcline_dc < mp.nme_legacy_dcline & mp.form_dc
% mp.nme_legacy_dcline_dc - Network model element for legacy DC line for DC formulation.
%
% Adds non-voltage state variables ``Pdcf`` and ``Pdct`` to the network
% model and builds the parameter :math:`\KK`.

%   MATPOWER
%   Copyright (c) 2019-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        function obj = add_zvars(obj, nm, dm, idx)
            %
            ndc = obj.nk;
            dme = obj.data_model_element(dm);
            switch idx{1}
                case 1      % flow at "from"
                    nm.add_var('z', 'Pdcf', ndc, dme.p_fr_start, dme.p_fr_lb, dme.p_fr_ub);
                case 2      % flow at "to"
                    nm.add_var('z', 'Pdct', ndc, dme.p_to_start, -Inf, Inf);
            end
        end

        function obj = build_params(obj, nm, dm)
            %
            build_params@mp.nme_legacy_dcline(obj, nm, dm); %% call parent
            obj.K = speye(obj.nk * obj.nz);
        end
    end     %% methods
end         %% classdef
