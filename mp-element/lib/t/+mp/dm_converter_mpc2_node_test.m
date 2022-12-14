classdef dm_converter_mpc2_node_test < mp.dm_converter_mpc2
%MP.DM_CONVERTER_MPC2_NODE_TEST  MATPOWER data model converter for MATPOWER case v2.

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        %% constructor
        function obj = dm_converter_mpc2_node_test()
            %% call parent constructor
            obj@mp.dm_converter_mpc2();
            obj.element_classes = ...
                { @mp.dmce_bus_nld_mpc2_node_test, ...
                  @mp.dmce_bus_ld_mpc2_node_test, ...
                  @mp.dmce_gen_mpc2, @mp.dmce_branch_mpc2 };
        end
    end     %% methods
end         %% classdef
