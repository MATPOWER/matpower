classdef net_model_dc < mp.net_model & mp.form_dc
% mp.net_model_dc - Concrete class for |MATPOWER| DC **network model** objects.
%
% This network model class and all of its network model element classes are
% specific to the DC formulation and therefore inherit from mp.form_dc.
%
% mp.net_model_dc Properties:
%   * va - vector of voltage states (voltage angles :math:`\Va`)
%   * z - vector of non-voltage states :math:`\z`
%
% mp.net_model_dc Methods:
%   * net_model_dc - constructor, assign default network model element classes
%   * def_set_types - add voltage and non-voltage variable set types for mp_idx_manager
%   * build_params - build incidence matrices, parameters, add ports for each element
%   * port_inj_soln - compute the network port injections at the solution
%
% See also mp.net_model, mp.form_dc, mp.form, mp.nm_element.

%   MATPOWER
%   Copyright (c) 2019-2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        va = [];    % *(double)* vector of voltage states (voltage angles :math:`\Va`)
        z  = [];    % *(double)* vector of non-voltage states :math:`\z`
    end

    methods
        function obj = net_model_dc()
            % Constructor, assign default network model element classes.
            % ::
            %
            %   nm = net_model_dc()
            %
            % This network model class and all of its network model element
            % classes are specific to the DC formulation and therefore
            % inherit from mp.form_dc.

            obj@mp.net_model();
            obj.element_classes = { ...
                @mp.nme_bus_dc, @mp.nme_gen_dc, @mp.nme_load_dc, ...
                    @mp.nme_branch_dc, @mp.nme_shunt_dc };

            %% Due to a bug related to inheritance in constructors in
            %% Octave 5.2 and earlier (https://savannah.gnu.org/bugs/?52614),
            %% INIT_SET_TYPES() cannot be called directly in the
            %% MP_IDX_MANAGER constructor, as desired.
            %%
            %% WORKAROUND:  INIT_SET_TYPES() is called explicitly as needed
            %%              (if obj.node is empty) in BUILD() and DISPLAY(),
            %%              after object construction, but before object use.
        end

        function obj = def_set_types(obj)
            % Add voltage and non-voltage variable set types for mp_idx_manager.
            % ::
            %
            %   nm.def_set_types()
            %
            % Add the following set types:
            %
            %   - ``'va'`` - VOLTAGE VARS (va)
            %   - ``'z'`` - NON-VOLTAGE VARS (z)
            %
            % See also mp.net_model.def_set_types, mp_idx_manager.

            def_set_types@mp.net_model(obj);        %% call parent first
            obj.set_types.va = 'VOLTAGE VARS (va)';
            obj.set_types.z  = 'NON-VOLTAGE VARS (z)';
        end

        function obj = build_params(obj, nm, dm)
            % Build incidence matrices and parameters, and add ports for each element.
            % ::
            %
            %   nm.build_params(nm, dm)
            %
            % Inputs:
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %
            % Call the parent method to do most of the work, then build
            % the aggregate network model parameters.

            %% call parent to build individual element parameters
            build_params@mp.net_model(obj, nm, dm);

            %% aggregate parameters from individual elements
            obj.B = obj.stack_matrix_params('B', 1);
            obj.K = obj.stack_matrix_params('K', 0);
            obj.p = obj.stack_vector_params('p');
        end

        function obj = port_inj_soln(obj)
            % Compute the network port injections at the solution.
            % ::
            %
            %   nm.port_inj_soln()
            %
            % Takes the solved network state, computes the port power
            % injections, and saves them in ``nm.soln.gp``.

            %% compute port injections
            obj.soln.gp = obj.port_inj_power(obj.soln.x);
        end
    end     %% methods
end         %% classdef
