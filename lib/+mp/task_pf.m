classdef task_pf < mp.task
% mp.task_pf - |MATPOWER| task for power flow (PF).
%
% Provides task implementation for the power flow problem.
%
% This includes the handling of iterative runs to enforce generator
% reactive power limits, if requested.
%
% mp.task_pf Properties:
%   * tag - task tag 'PF'
%   * name - task name 'Power Flow'
%   * dc - ``true`` if using DC network model
%   * iterations - total number of power flow iterations
%   * ref - current ref node indices
%   * ref0 - initial ref node indices
%   * va_ref0 - initial ref node voltage angles
%   * fixed_q_idx - indices of fixed Q gens
%   * fixed_q_qty - Q output of fixed Q gens
%
% mp.task_pf Methods:
%   * run_pre - set :attr:`dc` property
%   * next_dm - optionally iterate to enforce generator reactive limits
%   * enforce_q_lims - implementation of generator reactive limits
%   * network_model_class_default - select default network model constructor
%   * network_model_build_post - initialize properties for reactive limits
%   * network_model_x_soln - correct the voltage angles if necessary
%   * math_model_class_default - select default math model constructor
%
% See also mp.task.

%   MATPOWER
%   Copyright (c) 2020-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        tag = 'PF';             % 
        name = 'Power Flow';    % 
        dc  % ``true`` if using DC network model (from ``mpopt.model``, cached in run_pre())
        iterations              % *(integer)* total number of power flow iterations
        ref                     % *(integer)* current ref node indices
        ref0                    % *(integer)* initial ref node indices
        va_ref0                 % *(double)* initial ref node voltage angles
        fixed_q_idx             % *(integer)* indices of fixed Q gens
        fixed_q_qty             % *(double)* Q output of fixed Q gens
    end

    methods
        %%-----  task methods  -----
        function [d, mpopt] = run_pre(obj, d, mpopt)
            % Set :attr:`dc` property after calling superclass
            % :meth:`run_pre() <mp.task.run_pre>`.
            [d, mpopt] = run_pre@mp.task(obj, d, mpopt);    %% call parent

            %% cache DC model flag
            obj.dc = strcmp(upper(mpopt.model), 'DC');
        end

        function dm = next_dm(obj, mm, nm, dm, mpopt, mpx)
            % Implement optional iterations to enforce generator reactive
            % limits.
            if ~obj.dc && mpopt.pf.enforce_q_lims
                %% adjust iteration count for previous runs
                obj.iterations = obj.iterations + mm.soln.output.iterations;
                mm.soln.output.iterations = obj.iterations;

                %% enforce Q limits
                [success, dm] = obj.enforce_q_lims(nm, dm, mpopt);
                if ~success                 %% entire task fails if Q lim
                    obj.success = success;  %% enforcement indicates failure
                end
            else        %% don't enforce generator Q limits, once is enough
                dm = [];
            end
        end

        function [success, dm] = enforce_q_lims(obj, nm, dm, mpopt)
            % Used by next_dm() to implement enforcement of generator
            % reactive limits.
            gen_dme = dm.elements.gen;
            [mn, mx, both] = gen_dme.violated_q_lims(dm, mpopt);

            if ~isempty(both)   %% we have some Q limit violations
                if isempty(mn) && isempty(mx)   %% infeasible
                    if mpopt.verbose
                        fprintf('All %d remaining gens exceed their Q limits : INFEASIBLE PROBLEM\n', length(both));
                    end
                    dm = [];
                    success = 0;
                else
                    if mpopt.verbose && ~isempty(mx)
                        fprintf('Gen %d at upper Q limit, converting to PQ bus\n', gen_dme.on(mx));
                    end
                    if mpopt.verbose && ~isempty(mn)
                        fprintf('Gen %d at lower Q limit, converting to PQ bus\n', gen_dme.on(mn));
                    end

                    %% save corresponding limit values
                    obj.fixed_q_qty(mx) = gen_dme.qg_ub(mx);
                    obj.fixed_q_qty(mn) = gen_dme.qg_lb(mn);
                    mx = [mx;mn];

                    %% set qg to binding limit
                    gen_dme.tab.qg(gen_dme.on(mx)) = ...
                        obj.fixed_q_qty(mx) * dm.base_mva;

                    %% convert to PQ bus
                    bus_dme = dm.elements.bus;
                    ref0 = find(bus_dme.type == mp.NODE_TYPE.REF);
                    bidx = bus_dme.i2on(gen_dme.bus(gen_dme.on(mx)));   %% bus of mx
                    if length(ref0) > 1 && any(bus_dme.type(bidx) == mp.NODE_TYPE.REF)
                        error('mp.task_pf.enforce_q_lims: Sorry, MATPOWER cannot enforce Q limits for slack buses in systems with multiple slacks.');
                    end
                    %% set bus type to PQ
                    bus_dme.set_bus_type_pq(dm, bidx);
                    %% potentially pick new slack bus
                    ntv = nm.node_types(nm, dm);        %% node type vector
                    [i1, iN] = nm.get_node_idx('bus');  %% bus node indices
                    btv = ntv(i1:iN);                   %% bus type vector

                    %% indicate if there's been a change in slack bus
                    ref = find(btv == mp.NODE_TYPE.REF);    %% new ref bus indices
                    if mpopt.verbose && ref ~= ref0
                        fprintf('Bus %d is new slack bus\n', ...
                            bus_dme.ID(bus_dme.on(ref)));
                    end

                    %% save indices to list of Q limited gens
                    obj.fixed_q_idx = [obj.fixed_q_idx; mx];

                    %% update dm for next step
                    dm.initialize();
                    dm.update_status();
                    dm.build_params();
                    success = 1;
                end
            else                %% no more Q violations
                dm = [];
                success = 1;
            end
        end

        %%-----  data model methods  -----

        %%-----  network model methods  -----
        function nm_class = network_model_class_default(obj, dm, mpopt)
            % Implement selector for default network model constructor
            % depending on ``mpopt.model`` and ``mpopt.pf.v_cartesian``.
            switch upper(mpopt.model)
                case 'AC'
                    if mpopt.pf.v_cartesian
                        nm_class = @mp.net_model_acc;
                    else
                        nm_class = @mp.net_model_acp;
                    end
                case 'DC'
                    nm_class = @mp.net_model_dc;
            end
        end

        function nm = network_model_build_post(obj, nm, dm, mpopt)
            % Initialize mp.task_pf properties, if non-empty AC case with
            % generator reactive limits enforced.
            if ~obj.dc && mpopt.pf.enforce_q_lims ~= 0 && nm.np ~= 0
                if obj.i_nm == 1
                    [ref, ~, ~] = nm.node_types(obj, dm);
                    gen_dme =  dm.elements.gen;
                    obj.iterations = 0;
                    obj.ref0 = ref;             %% initial ref node indices
                    obj.ref = ref;              %% current ref node indices
                    obj.va_ref0 = nm.get_va(ref);%% initial ref node voltage angles
                    obj.fixed_q_idx = [];       %% indices of fixed Q gens
                    obj.fixed_q_qty = zeros(gen_dme.n, 1);  %% Q output of fixed Q gens
                else        %% update index of ref bus
                    [obj.ref, ~, ~] = nm.node_types(obj, dm);
                end
            end
        end

        function nm = network_model_x_soln(obj, mm, nm)
            % Call superclass :meth:`network_model_x_soln() <mp.task.network_model_x_soln>`
            % then correct the voltage angle if the ref node has been changed.
            nm = network_model_x_soln@mp.task(obj, mm, nm);

            if ~obj.dc && obj.i_nm > 1 && obj.ref ~= obj.ref0
                vm = abs(nm.soln.v);
                va = angle(nm.soln.v);
                va = va - va(obj.ref0) + obj.va_ref0;
                nm.soln.v = vm .* exp(1j * va);
            end
        end

        %%-----  mathematical model methods  -----
        function mm_class = math_model_class_default(obj, nm, dm, mpopt)
            % Implement selector for default mathematical model constructor
            % depending on ``mpopt.model``, ``mpopt.pf.v_cartesian``, and
            % ``mpopt.pf.current_balance``.
            switch upper(mpopt.model)
                case 'AC'
                    if mpopt.pf.v_cartesian
                        if mpopt.pf.current_balance
                            mm_class = @mp.math_model_pf_acci;
                        else
                            mm_class = @mp.math_model_pf_accs;
                        end
                    else
                        if mpopt.pf.current_balance
                            mm_class = @mp.math_model_pf_acpi;
                        else
                            mm_class = @mp.math_model_pf_acps;
                        end
                    end
                case 'DC'
                    mm_class = @mp.math_model_pf_dc;
            end
        end
    end     %% methods
end         %% classdef
