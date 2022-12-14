classdef xt_node_test < mp.extension

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
        function dmc_class = dm_converter_class(obj, dmc_class, fmt, mpopt)
            dmc_class = @mp.dm_converter_mpc2_node_test;
        end

        function dm_class = data_model_class(obj, dm_class, task_tag, mpopt)
            dm_class = @mp.data_model_node_test;
        end

        function nm_class = network_model_class(obj, nm_class, task_tag, mpopt)
            nm_class = @mp.net_model_acp_node_test;
        end

        function mm_class = math_model_class(obj, mm_class, task_tag, mpopt)
            switch task_tag
                case {'PF'}
                    mm_class = @mp.math_model_pf_acps_node_test;
                case {'OPF'}
                    mm_class = @mp.math_model_opf_acps_node_test;
            end
        end
    end     %% methods
end         %% classdef
