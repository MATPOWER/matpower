classdef xt_reserves < mp.extension

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
        function dmc_elements = dmc_element_classes(obj, dmc_class, fmt, mpopt)
            switch fmt
                case 'mpc2'
                    dmc_elements = { @mp.dmce_reserve_gen_mpc2, ...
                                     @mp.dmce_reserve_zone_mpc2 };
                otherwise
                    dmc_elements = {};      %% no modifications
            end
        end

        function dm_elements = dm_element_classes(obj, dm_class, task_tag, mpopt)
            switch task_tag
                case 'OPF'
                    dm_elements = { @mp.dme_reserve_gen, ...
                                    @mp.dme_reserve_zone };
                otherwise
                    dm_elements = {};       %% no modifications
            end
        end

        function mm_elements = mm_element_classes(obj, mm_class, task_tag, mpopt)
            switch task_tag
                case {'OPF'}
                    mm_elements = { @mp.mme_reserve_gen, ...
                                    @mp.mme_reserve_zone };
                otherwise
                    mm_elements = {};       %% no modifications
            end
        end
    end     %% methods
end         %% classdef
