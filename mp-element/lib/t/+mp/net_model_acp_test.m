classdef net_model_acp_test < mp.net_model_acp

%   MATPOWER
%   Copyright (c) 2019-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        %% constructor
        function obj = net_model_acp_test()
            obj@mp.net_model_acp();
            obj.element_classes{end+1} = @mp.nme_gizmo_acp;

            %% Due to a bug related to inheritance in constructors in
            %% Octave 5.2 and earlier (https://savannah.gnu.org/bugs/?52614),
            %% INIT_SET_TYPES() cannot be called directly in the
            %% MP_IDX_MANAGER constructor, as desired.
            %%
            %% WORKAROUND:  INIT_SET_TYPES() is called explicitly as needed
            %%              (if obj.node is empty) in BUILD() and DISPLAY(),
            %%              after object construction, but before object use.
        end
    end     %% methods
end         %% classdef
