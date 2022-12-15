classdef net_model_dc < mp.net_model & mp.form_dc

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        va = [];
        z  = [];
    end

    methods
        %% constructor
        function obj = net_model_dc()
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
            def_set_types@mp.net_model(obj);        %% call parent first
            obj.set_types.va = 'VOLTAGE VARS (va)';
            obj.set_types.z  = 'NON-VOLTAGE VARS (z)';
        end

        function obj = build_params(obj, nm, dm)
            %% call parent to build individual element parameters
            build_params@mp.net_model(obj, nm, dm);

            %% aggregate parameters from individual elements
            obj.B = obj.stack_matrix_params('B', 1);
            obj.K = obj.stack_matrix_params('K', 0);
            obj.p = obj.stack_vector_params('p');
        end

        function obj = port_inj_soln(obj)
            %% compute port injections
            obj.soln.gp = obj.port_inj_power(obj.soln.x);
        end
    end     %% methods
end         %% classdef
