classdef data_model_cpf < mp.data_model
%MP.DATA_MODEL_CPF  MATPOWER data model for CPF tasks

%   MATPOWER
%   Copyright (c) 2020-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
    end     %% properties

    methods
        %% constructor
        function obj = data_model_cpf()
            %% call parent constructor
            obj@mp.data_model();
            obj.element_classes = ...
                { @mp.dme_bus, @mp.dme_gen, @mp.dme_load_cpf, ...
                    @mp.dme_branch, @mp.dme_shunt_cpf };
        end
    end     %% methods
end         %% classdef
