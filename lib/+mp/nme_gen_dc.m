classdef nme_gen_dc < mp.nme_gen & mp.form_dc
% mp.nme_gen_dc - Network model element for generator for DC formulations.
%
% Adds non-voltage state variable ``Pg`` to the network model, builds the
% parameter :math:`\KK`, and inherits from mp.form_dc.

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
        function obj = add_zvars(obj, nm, dm, idx)
            %
            ng = obj.nk;
            dme = obj.data_model_element(dm);
            nm.add_var('z', 'Pg', ng, dme.pg_start, dme.pg_lb, dme.pg_ub);
        end

        function obj = build_params(obj, nm, dm)
            %
            build_params@mp.nme_gen(obj, nm, dm);   %% call parent
            obj.K = -speye(obj.nk * obj.nz);
        end
    end     %% methods
end         %% classdef
