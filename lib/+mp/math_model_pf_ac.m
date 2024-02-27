classdef (Abstract) math_model_pf_ac < mp.math_model_pf
% mp.math_model_pf_ac - Power flow (PF) **math model** for AC formulations.
%
% Provides AC-specific and PF-specific subclasses for elements.

%   MATPOWER
%   Copyright (c) 2021-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        %% constructor
        function obj = math_model_pf_ac()
            %

            obj@mp.math_model_pf();
            obj.element_classes = { @mp.mme_bus_pf_ac, @mp.mme_gen_pf_ac, ...
                @mp.mme_load_pf_ac, @mp.mme_branch_pf_ac, @mp.mme_shunt_pf_ac };
        end
    end     %% methods
end         %% classdef
