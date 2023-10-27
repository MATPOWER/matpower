classdef xt_reserves < mp.extension
% mp.xt_reserves - |MATPOWER| extension for OPF with fixed zonal reserves.
%
% For OPF problems, this extension adds two types of elements to the
% data and mathematical model containers, as well as the data model converter.
%
% The ``'reserve_gen'`` element handles all of the per-generator aspects,
% such as reserve cost and quantity limit parameters, reserve variables,
% and constraints on reserve capacity.
%
% The ``'reserve_zone'`` element handles the per-zone aspects, such as
% generator/zone mappings, zonal reserve requirement parameters and
% constraints, and zonal reserve prices.
%
% mp.xt_reserves Methods:
%   * dmc_element_classes - add two classes to data model converter elements
%   * dm_element_classes - add two classes to  data model elements
%   * mm_element_classes - add two classes to mathematical model elements
%
% See the :ref:`sec_customizing` and :ref:`sec_extensions` sections in the
% |MATPOWER-Dev-Manual| for more information, and specifically the
% :ref:`sec_element_classes` section and the :ref:`tab_element_class_modifiers`
% table for details on *element class modifiers*.
%
% See also mp.extension.

%   MATPOWER
%   Copyright (c) 2022-2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function dmc_elements = dmc_element_classes(obj, dmc_class, fmt, mpopt)
            % Add two classes to data model converter elements.
            %
            % For ``'mpc2`` data formats, adds the classes:
            %
            %   - mp.dmce_reserve_gen_mpc2
            %   - mp.dmce_reserve_zone_mpc2

            switch fmt
                case 'mpc2'
                    dmc_elements = { @mp.dmce_reserve_gen_mpc2, ...
                                     @mp.dmce_reserve_zone_mpc2 };
                otherwise
                    dmc_elements = {};      %% no modifications
            end
        end

        function dm_elements = dm_element_classes(obj, dm_class, task_tag, mpopt)
            % Add two classes to data model elements.
            %
            % For ``'OPF'`` tasks, adds the classes:
            %
            %   - mp.dme_reserve_gen
            %   - mp.dme_reserve_zone

            switch task_tag
                case 'OPF'
                    dm_elements = { @mp.dme_reserve_gen, ...
                                    @mp.dme_reserve_zone };
                otherwise
                    dm_elements = {};       %% no modifications
            end
        end

        function mm_elements = mm_element_classes(obj, mm_class, task_tag, mpopt)
            % Add two classes to mathematical model elements.
            %
            % For ``'OPF'`` tasks, adds the classes:
            %
            %   - mp.mme_reserve_gen
            %   - mp.mme_reserve_zone

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
