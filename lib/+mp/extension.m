classdef (Abstract) extension < handle

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
        function task_class = task_class(obj, task_class, mpopt)
        end

        function dmc_class = dm_converter_class(obj, dmc_class, fmt, mpopt)
        end

        function dm_class = data_model_class(obj, dm_class, task_tag, mpopt)
        end

        function nm_class = network_model_class(obj, nm_class, task_tag, mpopt)
        end

        function mm_class = math_model_class(obj, mm_class, task_tag, mpopt)
        end

        function dmc_elements = dmc_element_classes(obj, dmc_class, fmt, mpopt)
            dmc_elements = {};      %% no modifications
        end

        function dm_elements = dm_element_classes(obj, dm_class, task_tag, mpopt)
            dm_elements = {};       %% no modifications
        end

        function nm_elements = nm_element_classes(obj, nm_class, task_tag, mpopt)
            nm_elements = {};       %% no modifications
        end

        function mm_elements = mm_element_classes(obj, mm_class, task_tag, mpopt)
            mm_elements = {};       %% no modifications
        end
    end     %% methods
end         %% classdef
