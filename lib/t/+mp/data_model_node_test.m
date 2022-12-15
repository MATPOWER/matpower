classdef data_model_node_test < mp.data_model
%MP.DATA_MODEL_NODE_TEST  Data model class for node test.

%   MATPOWER
%   Copyright (c) 2020-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        %% constructor
        function obj = data_model_node_test()
            %% call parent constructor
            obj@mp.data_model();
            obj.element_classes = ...
                { @mp.dme_bus_nld_node_test, @mp.dme_bus_ld_node_test, ...
                    @mp.dme_gen_node_test, @mp.dme_branch_node_test };
        end
    end     %% methods
end         %% classdef
