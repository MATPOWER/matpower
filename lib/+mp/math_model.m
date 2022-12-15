classdef (Abstract) math_model < mp.element_container & opt_model
%MP.MATH_MODEL  MATPOWER mathematical model abstract base class.
%   ?
%
%   MP.MATH_MODEL provides properties and methods related to the specific
%   problem specification being solved (e.g. power flow, continuation
%   power flow, optimal power flow, etc.) ...
%
%   Properties
%       ? - ?
%
%   Methods
%       ?

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        aux_data    %% struct of auxiliary data relevant to the model,
                    %% e.g. can be passed to model constraint functions
    end

    methods
        function tag = task_tag(obj)
            error('mp.math_model/task_tag: must be implemented in subclass');
        end

        function name = task_name(obj)
            error('mp.math_model/task_name: must be implemented in subclass');
        end

        function tag = form_tag(obj)
            error('mp.math_model/form_tag: must be implemented in subclass');
        end

        function name = form_name(obj)
            error('mp.math_model/form_name: must be implemented in subclass');
        end

        function obj = build(obj, nm, dm, mpopt)
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
            %% create aux_data struct
            obj.aux_data = obj.build_aux_data(nm, dm, mpopt);
        end

        function ad = build_base_aux_data(obj, nm, dm, mpopt)
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
            obj.add_system_vars(nm, dm, mpopt);

            %% each element adds its OPF variables
            for k = 1:length(obj.elements)
                obj.elements{k}.add_vars(obj, nm, dm, mpopt);
            end
        end

        function obj = add_system_vars(obj, nm, dm, mpopt)
        end

        function obj = add_constraints(obj, nm, dm, mpopt)
            obj.add_system_constraints(nm, dm, mpopt);

            %% each element adds its OPF variables
            for k = 1:length(obj.elements)
                obj.elements{k}.add_constraints(obj, nm, dm, mpopt);
            end
        end

        function obj = add_system_constraints(obj, nm, dm, mpopt)
            %% node balance constraints
            obj.add_node_balance_constraints(nm, dm, mpopt);
        end

        function obj = add_node_balance_constraints(obj, nm, dm, mpopt)
        end

        function obj = add_costs(obj, nm, dm, mpopt)
            obj.add_system_costs(nm, dm, mpopt);

            %% each element adds its OPF variables
            for k = 1:length(obj.elements)
                obj.elements{k}.add_costs(obj, nm, dm, mpopt);
            end
        end

        function obj = add_system_costs(obj, nm, dm, mpopt)
        end

        function opt = solve_opts(obj, nm, dm, mpopt)
            opt = struct();
        end

        function nm_vars = update_nm_vars(obj, mmx, nm)
            %% Returns a struct with the network model variables as fields.
            %% The ad.var_map cell array is used to track mappings
            %% of math model variables back to network model variables.
            %% Each entry is a cell array:
            %%  {nm_var_type, nm_i1, nm_iN, nm_idx, mm_i1, mm_iN, mm_idx}
            %% where
            %%  nm_var_type - network model variable type (e.g. va, vm, zr, zi)
            %%  nm_i1 - starting index for network model variable type
            %%  nm_iN - ending index for network model variable type
            %%  nm_idx - vector of indices for network model variable type
            %%  mm_i1 - starting index for math model variable
            %%  mm_iN - ending index for math model variable
            %%  mm_idx - vector of indices for math model variable
            %% Uses either i1:iN (if i1 is not empty) or idx as the indices,
            %% unless both are empty, in which case it uses ':'

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
            %% each element updates its data model
            for k = 1:length(obj.elements)
                obj.elements{k}.data_model_update(obj, nm, dm, mpopt);
            end
        end

        function nm = network_model_x_soln(obj, nm)
            %% convert solved state from math model to network model soln
            [nm.soln.v, nm.soln.z, nm.soln.x] = ...
                obj.convert_x_m2n(obj.soln.x, nm);
        end
    end     %% methods
end         %% classdef
