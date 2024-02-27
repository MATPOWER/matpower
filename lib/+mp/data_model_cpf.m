classdef data_model_cpf < mp.data_model
% mp.data_model_cpf - |MATPOWER| **data model** for CPF tasks.
%
% The purpose of this class is to include CPF-specific subclasses for the
% load and shunt elements, which need to be able to provide versions of
% their model parameters that are parameterized by the continuation
% parameter :math:`\lambda`.
%
% data_model_cpf Methods:
%   * data_model_cpf - constructor, assign default data model element classes
%
% See also mp.data_model.

%   MATPOWER
%   Copyright (c) 2020-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        %% constructor
        function obj = data_model_cpf()
            % Constructor, assign default data model element classes.
            %
            % Create an empty data model object and assign the default
            % data model element classes, which are the same as those
            % defined by the base class, except for loads and shunts.
            % ::
            %
            %   dm = mp.data_model_cpf()

            %% call parent constructor
            obj@mp.data_model();
            obj.element_classes = ...
                { @mp.dme_bus, @mp.dme_gen, @mp.dme_load_cpf, ...
                    @mp.dme_branch, @mp.dme_shunt_cpf };
        end
    end     %% methods
end         %% classdef
