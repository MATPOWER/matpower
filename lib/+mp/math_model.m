classdef (Abstract) math_model < mp.element_container & opt_model
% mp.math_model - Abstract base class for |MATPOWER| **mathematical model** objects.
%
% The mathematical model, or math model, formulates and defines the
% mathematical problem to be solved. That is, it determines the variables,
% constraints, and objective that define the problem. This takes on different
% forms depending on the task *(e.g. power flow, optimal power flow, etc.)*
% and the formulation *(e.g. DC, AC-polar-power, etc.)*.
%
% A math model object is a container for math model element (mp.mm_element)
% objects and it is also an MP-Opt-Model (opt_model) object. All math model
% classes inherit from mp.math_model and therefore also from
% mp.element_container, opt_model, and mp_idx_manager. Concrete math model
% classes are task and formulation specific. They also sometimes inherit from
% abstract mix-in classes that are shared across tasks or formulations.
%
% By convention, math model variables are named ``mm`` and math model class
% names begin with ``mp.math_model``.
%
% mp.math_model Properties:
%   * aux_data - auxiliary data relevant to the model
%
% mp.math_model Methods:
%   * task_tag - returns task tag, e.g. ``'PF'``, ``'OPF'``
%   * task_name - returns task name, e.g. ``'Power Flow'``, ``'Optimal Power Flow'``
%   * form_tag - returns network formulation tag, e.g. ``'dc'``, ``'acps'``
%   * form_name - returns network formulation name, e.g. ``'DC'``, ``'AC-polar-power'``
%   * build - create, add, and build math model element objects
%   * display - display the math model object
%   * add_aux_data - builds auxiliary data and adds it to the model
%   * build_base_aux_data - builds base auxiliary data, including node types & variable initial values
%   * add_vars - add variables to the model
%   * add_system_vars - add system variables to the model
%   * add_constraints - add constraints to the model
%   * add_system_constraints - add system constraints to the model
%   * add_node_balance_constraints - add node balance constraints to the model
%   * add_costs - add costs to the model
%   * add_system_costs - add system costs to the model
%   * solve_opts - return an options struct to pass to the solver
%   * update_nm_vars - update network model variables from math model solution
%   * data_model_update - update data model from math model solution
%   * network_model_x_soln - convert solved state from math model to network model solution
%
% See the :ref:`sec_math_model` section in the |MATPOWER-Dev-Manual| for more
% information.
%
% See also mp.task, mp.data_model, mp.net_model.

%   MATPOWER
%   Copyright (c) 2021-2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        % *(struct)* auxiliary data relevant to the model, e.g. can be
        % passed to model constraint functions
        aux_data
    end

    methods
        function tag = task_tag(obj)
            % Returns task tag, e.g. ``'PF'``, ``'OPF'``.
            % ::
            %
            %   tag = mm.task_tag()
            
            error('mp.math_model.task_tag: must be implemented in subclass');
        end

        function name = task_name(obj)
            % Returns task name, e.g. ``'Power Flow'``, ``'Optimal Power Flow'``.
            % ::
            %
            %   name = mm.task_name()

            error('mp.math_model.task_name: must be implemented in subclass');
        end

        function tag = form_tag(obj)
            % Returns network formulation tag, e.g. ``'dc'``, ``'acps'``.
            % ::
            %
            %   tag = mm.form_tag()

            error('mp.math_model.form_tag: must be implemented in subclass');
        end

        function name = form_name(obj)
            % Returns network formulation name, e.g. ``'DC'``, ``'AC-polar-power'``.
            % ::
            %
            %   name = mm.form_name()

            error('mp.math_model.form_name: must be implemented in subclass');
        end

        function obj = build(obj, nm, dm, mpopt)
            % Create, add, and build math model element objects.
            % ::
            %
            %   mm.build(nm, dm, mpopt);
            %
            % Inputs:
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Create and add network model objects, create and add auxiliary
            % data, and add variables, constraints, and costs.

            %% Due to a bug related to inheritance in constructors in
            %% Octave 5.2 and earlier (https://savannah.gnu.org/bugs/?52614),
            %% init_set_types() cannot be called directly in the
            %% MP_IDX_MANAGER constructor, as desired.
            %%
            %% WORKAROUND:  Initialize MP_IDX_MANAGER fields here, if needed,
            %%              after object construction, but before object use.
            if isempty(obj.var)         %% only if not already initialized
                obj.init_set_types();
            end

            %% create element objects for each class with data
            obj.elements = mp.mapped_array();
            for c = obj.element_classes
                mme = c{1}();       %% element constructor
                if dm.online(mme.name)      %% dm element exists
                    obj.elements.add_elements(mme, mme.name);
                end
            end

            obj.add_aux_data(nm, dm, mpopt);
            obj.add_vars(nm, dm, mpopt);
            obj.add_constraints(nm, dm, mpopt);
            obj.add_costs(nm, dm, mpopt);
        end

        function display(obj)
            % Display the math model object.
            %
            % This method is called automatically when omitting a semicolon
            % on a line that retuns an object of this class.
            %
            % Displays the details of the variables, constraints, costs, and
            % math model elements.
            %
            % See also mp_idx_manager.

            fprintf('MATH MODEL CLASS : %s\n', class(obj));
            fprintf('    TASK NAME               : %s\n', obj.task_name());
            fprintf('    TASK TAG                : %s\n', obj.task_tag());
            fprintf('    FORMULATION NAME        : %s\n', obj.form_name());
            fprintf('    FORMULATION TAG         : %s\n', obj.form_tag());
            display@opt_model(obj)

            %% elements
            fprintf('\nELEMENTS\n')
            fprintf('========\n')
            fprintf('  name              class\n');
            fprintf(' ----------------  --------------------\n');
            for k = 1:length(obj.elements)
                mme = obj.elements{k};
                fprintf('  %-13s     %s\n', mme.name, class(mme));
            end
        end

        function obj = add_aux_data(obj, nm, dm, mpopt)
            % Builds auxiliary data and adds it to the model.
            % ::
            %
            %   mm.add_aux_data(nm, dm, mpopt)
            %
            % Inputs:
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Calls the build_aux_data() method and assigns the
            % result to the :attr:`aux_data` property. The base build_aux_data()
            % method, which simply calls build_base_aux_data(), is defined
            % in mp.mm_shared_pfcpf (and in mp.math_model_opf) allowing it to
            % be shared across math models for different tasks (PF and CPF).

            obj.aux_data = obj.build_aux_data(nm, dm, mpopt);
        end

        function ad = build_base_aux_data(obj, nm, dm, mpopt)
            % Builds base auxiliary data, including node types & variable initial values.
            % ::
            %
            %   ad = mm.build_base_aux_data(nm, dm, mpopt)
            %
            % Inputs:
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Output:
            %   ad (struct) : struct of auxiliary data

            %% get model variables
            vvars = nm.model_vvars();
            zvars = nm.model_zvars();
            vars = {vvars{:} zvars{:}};
            vals = cellfun(@(x)nm.params_var(x), vars, 'UniformOutput', false);

            %% get node types
            [ref, pv, pq, by_elm] = nm.node_types(nm, dm);

            %% create aux_data struct
            ad = cell2struct(...
                {vals{:}, {}, ref, pv, pq, ...
                    length(ref), length(pv), length(pq), by_elm}, ...
                {vars{:}, 'var_map', 'ref', 'pv', 'pq', ...
                    'nref', 'npv', 'npq', 'node_type_by_elm'}, 2);
        end

        function obj = add_vars(obj, nm, dm, mpopt)
            % Add variables to the model.
            % ::
            %
            %   mm.add_vars(nm, dm, mpopt)
            %
            % Inputs:
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Adds system variables, then calls the
            % :meth:`add_vars() <mp.mm_element.add_vars>` method for each
            % math model element.

            obj.add_system_vars(nm, dm, mpopt);

            %% each element adds its OPF variables
            for k = 1:length(obj.elements)
                obj.elements{k}.add_vars(obj, nm, dm, mpopt);
            end
        end

        function obj = add_system_vars(obj, nm, dm, mpopt)
            % Add system variables to the model.
            % ::
            %
            %   mm.add_system_vars(nm, dm, mpopt)
            %
            % Inputs:
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Variables which correspond to a specific math model element
            % should be added by that element's
            % :meth:`add_vars() <mp.mm_element.add_vars>` method. Other
            % variables can be added by add_system_vars(). In this base class
            % this method does nothing.
        end

        function obj = add_constraints(obj, nm, dm, mpopt)
            % Add constraints to the model.
            % ::
            %
            %   mm.add_constraints(nm, dm, mpopt)
            %
            % Inputs:
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Adds system constraints, then calls the
            % :meth:`add_constraints() <mp.mm_element.add_constraints>` method
            % for each math model element.

            obj.add_system_constraints(nm, dm, mpopt);

            %% each element adds its OPF variables
            for k = 1:length(obj.elements)
                obj.elements{k}.add_constraints(obj, nm, dm, mpopt);
            end
        end

        function obj = add_system_constraints(obj, nm, dm, mpopt)
            % Add system constraints to the model.
            % ::
            %
            %   mm.add_system_constraints(nm, dm, mpopt)
            %
            % Inputs:
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Constraints which correspond to a specific math model element
            % should be added by that element's
            % :meth:`add_constraints() <mp.mm_element.add_constraints>` method.
            % Other constraints can be added by add_system_constraints(). In
            % this base class, it simply calls add_node_balance_constraints().

            %% node balance constraints
            obj.add_node_balance_constraints(nm, dm, mpopt);
        end

        function obj = add_node_balance_constraints(obj, nm, dm, mpopt)
            % Add node balance constraints to the model.
            % ::
            %
            %   mm.add_node_balance_constraints(nm, dm, mpopt)
            %
            % Inputs:
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % In this base class this method does nothing.
        end

        function obj = add_costs(obj, nm, dm, mpopt)
            % Add costs to the model.
            % ::
            %
            %   mm.add_costs(nm, dm, mpopt)
            %
            % Inputs:
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Adds system costs, then calls the
            % :meth:`add_costs() <mp.mm_element.add_costs>` method for each
            % math model element.

            obj.add_system_costs(nm, dm, mpopt);

            %% each element adds its OPF variables
            for k = 1:length(obj.elements)
                obj.elements{k}.add_costs(obj, nm, dm, mpopt);
            end
        end

        function obj = add_system_costs(obj, nm, dm, mpopt)
            % Add system costs to the model.
            % ::
            %
            %   mm.add_system_costs(nm, dm, mpopt)
            %
            % Inputs:
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Costs which correspond to a specific math model element
            % should be added by that element's
            % :meth:`add_costs() <mp.mm_element.add_costs>` method. Other
            % variables can be added by add_system_costs(). In this base class
            % this method does nothing.
        end

        function opt = solve_opts(obj, nm, dm, mpopt)
            % Return an options struct to pass to the solver.
            % ::
            %
            %   opt = mm.solve_opts(nm, dm, mpopt)
            %
            % Inputs:
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Output:
            %   opt (struct) : options struct for solver
            %
            % In this base class, returns an empty struct.

            opt = struct();
        end

        function nm_vars = update_nm_vars(obj, mmx, nm)
            % Update network model variables from math model solution.
            % ::
            %
            %   nm_vars = mm.update_nm_vars(mmx, nm)
            %
            % Inputs:
            %   mmx (double) : vector of math model variable ``x``
            %   nm (mp.net_model) : network model object
            %
            % Output:
            %   nm_vars (struct) : updated network model variables
            %
            % Returns a struct with the network model variables as fields.
            % The ``mm.aux_data.var_map`` cell array is used to track mappings
            % of math model variables back to network model variables.
            % Each entry is itself a 7-element cell array of the form
            %
            % ``{nm_var_type, nm_i1, nm_iN, nm_idx, mm_i1, mm_iN, mm_idx}``
            %
            % where
            %
            %   - ``nm_var_type`` - network model variable type (e.g. ``va``, ``vm``, ``zr``, ``zi``)
            %   - ``nm_i1`` - starting index for network model variable type
            %   - ``nm_iN`` - ending index for network model variable type
            %   - ``nm_idx`` - vector of indices for network model variable type
            %   - ``mm_i1`` - starting index for math model variable
            %   - ``mm_iN`` - ending index for math model variable
            %   - ``mm_idx`` - vector of indices for math model variable
            %
            % Uses either ``i1:iN`` (if ``i1`` is not empty) or ``idx`` as the
            % indices, unless both are empty, in which case it uses ``':'``.

            ad = obj.aux_data;
            vvars = nm.model_vvars();
            zvars = nm.model_zvars();
            vars = {vvars{:} zvars{:}};
            nm_vars = cell2struct( ...
                cellfun(@(x)ad.(x), vars, 'UniformOutput', false), ...
                vars, 2);

            %% copy mm state values to nm state based on ad.var_map
            for k = 1:length(ad.var_map)
                d = ad.var_map{k};

                %% right hand side from mmx
                if ~isempty(d{5})       %% use i1:iN
                    rhs = mmx(d{5}:d{6});
                elseif isempty(d{7})    %% use :
                    rhs = mmx;
                else                    %% use idx
                    rhs = mmx(d{7});
                end

                %% left hand side to nm_vars
                name = d{1};
                if ~isempty(rhs)
                    if ~isempty(d{2})       %% use i1:iN
                        nm_vars.(name)(d{2}:d{3}) = rhs;
                    elseif isempty(d{4})    %% use :
                        nm_vars.(name) = rhs;
                    else                    %% use idx
                        nm_vars.(name)(d{4}) = rhs;
                    end
                end
            end
        end

        function dm = data_model_update(obj, nm, dm, mpopt)
            % Update data model from math model solution.
            % ::
            %
            %   dm = mm.data_model_update(nm, dm, mpopt)
            %
            % Inputs:
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Output:
            %   dm (mp.data_model) : updated data model object
            %
            % Calls the
            % :meth:`data_model_update() <mp.mm_element.data_model_update>`
            % method for each math model element.

            %% each element updates its data model
            for k = 1:length(obj.elements)
                obj.elements{k}.data_model_update(obj, nm, dm, mpopt);
            end
        end

        function nm = network_model_x_soln(obj, nm)
            % Convert solved state from math model to network model solution.
            % ::
            %
            %   nm = mm.network_model_x_soln(nm)
            %
            % Input:
            %   nm (mp.net_model) : network model object
            %
            % Output:
            %   nm (mp.net_model) : updated network model object
            %
            % Calls convert_x_m2n() to which is defined in a subclass of
            % in mp.mm_shared_pfcpf (and of mp.math_model_opf) allowing it to
            % be shared across math models for different tasks (PF and CPF).
            
            [nm.soln.v, nm.soln.z, nm.soln.x] = ...
                obj.convert_x_m2n(obj.soln.x, nm);
        end
    end     %% methods
end         %% classdef
