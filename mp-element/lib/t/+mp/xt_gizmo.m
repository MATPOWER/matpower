classdef xt_gizmo < mp.extension

%   MATPOWER
%   Copyright (c) 2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function nm_class = network_model_class(obj, nm_class, task_tag, mpopt)
            nm_class = @mp.net_model_acp_test;
        end

        function dmc_elements = dmc_element_classes(obj, dmc_class, fmt, mpopt)
            dmc_elements = { @mp.dmce_gizmo_mpc2 };
        end

        function dm_elements = dm_element_classes(obj, dm_class, task_tag, mpopt)
            dm_elements = { @mp.dme_gizmo };
        end
    end     %% methods
end         %% classdef
