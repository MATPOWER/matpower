classdef (Abstract) math_model_cpf_acc < mp.math_model_cpf
% mp.math_model_cpf_acc - Abstract base class for AC cartesian CPF **math model** objects.
%
% Provides formulation-specific and CPF-specific subclasses for elements.

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
        function obj = math_model_cpf_acc()
            % Constructor, assign default network model element classes.
            % ::
            %
            %   mm = math_model_cpf_acc()

            obj@mp.math_model_cpf();
            obj.element_classes = { @mp.mme_bus_pf_ac, @mp.mme_gen_pf_ac, ...
                @mp.mme_branch_pf_ac, @mp.mme_load_cpf, @mp.mme_shunt_cpf };
        end
    end     %% methods
end         %% classdef
