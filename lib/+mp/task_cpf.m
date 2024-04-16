classdef task_cpf < mp.task_pf
% mp.task_cpf - |MATPOWER| task for continuation power flow (CPF).
%
% Provides task implementation for the continuation power flow problem.
%
% This includes the iterative solving of the mathematical model (using
% warm restarts) after updating the problem data, e.g. when enforcing
% certain limits.
%
% mp.task_cpf Properties:
%   * warmstart - warm start data
%
% mp.task_cpf Methods:
%   * task_cpf - constructor, inherits from mp.task_pf constructor
%   * run_pre - call superclass :meth:`run_pre() <mp.task_pf.run_pre>` for base and target inputs
%   * next_mm - handle warm start of continuation iterations
%   * dm_converter_class - select data model converter class
%   * data_model_class_default - select default data model constructor
%   * data_model_build - build base and target data models
%   * network_model_build - build base and target network models
%   * network_model_x_soln - update network model solution
%   * network_model_update - evaluate port injection solution
%   * math_model_class_default - select default math model constructor
%   * math_model_opt - add warmstart parameters to math model solve options
%
% See also mp.task, mp.task_pf.

%   MATPOWER
%   Copyright (c) 2020-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        % *(struct)* warm start data, with fields:
        %
        %     - clam - corrector parameter lambda
        %     - plam - predictor parameter lambda
        %     - cV - corrector complex voltage vector
        %     - pV - predictor complex voltage vector
        warmstart
    end

    methods
        function obj = task_cpf()
            % Constructor, inherits from mp.task_pf constructor.

            %% call parent constructor
            obj@mp.task_pf();

            obj.tag = 'CPF';
            obj.name = 'Continuation Power Flow';
        end

        %%-----  task methods  -----
        function [d, mpopt] = run_pre(obj, d, mpopt)
            % Call superclass :meth:`run_pre() <mp.task_pf.run_pre>` for
            % base and target inputs.
            if ~isa(d, 'mp.data_model')
                if ~iscell(d) || length(d) < 2
                    error('mp.task_cpf.run_pre: input cases must be provided in a 2-element cell array, specifying the base and target cases, respectively')
                end
                d{1} = run_pre@mp.task_pf(obj, d{1}, mpopt);
                d{2} = run_pre@mp.task_pf(obj, d{2}, mpopt);
            end
        end

        function [mm, nm, dm] = next_mm(obj, mm, nm, dm, mpopt, mpx)
            % Handle warm start of continuation iterations, after problem
            % data update.

            % return new math model, or empty matrix if finished
            if isfield(mm.soln.output, 'warmstart')
                %% get warmstart info
                ad = mm.aux_data;
                ws = mm.soln.output.warmstart;

                %% save parameter lambda and solved voltages
                %% for current & prev step
                ws.clam = ws.x(end);
                ws.plam = ws.xp(end);
                [ws.cV, ~] = mm.convert_x_m2n_cpf(ws.x, nm);
                [ws.pV, ~] = mm.convert_x_m2n_cpf(ws.xp, nm);

                %% expand tangent z to all nodes + lambda, for cur & prev step
                [ws.z, ws.zp] = mm.expand_z_warmstart(nm, ad, ws.z, ws.zp);

                %% set warmstart for next math model
                obj.warmstart = ws;

                %% save updated target models
                nm.userdata.target = ws.nmt;
                dm.userdata.target = ws.dmt;

                %% update network model with current solution
                obj.nm = obj.network_model_update(mm, nm);

                %% update data model voltages only
                %% preserve original base/target specifications
                for k = 1:obj.nm.node.NS
                    nme = obj.nm.elements.(obj.nm.node.order(k).name);
                    mme = nme.math_model_element(mm);
                    mme.data_model_update(mm, obj.nm, obj.dm, mpopt);
                end

                %% reset var_map
                obj.nm.userdata.var_map = {};

                %% create new math model
                mm = obj.math_model_build(nm, dm, mpopt, mpx);
            else
                mm = [];
            end
        end

        %%-----  data model converter methods  -----
        function dmc_class = dm_converter_class(obj, d, mpopt, mpx)
            % Implement selector for data model converter class based on
            % superclass constructor.

            if iscell(d) && length(d) == 2
                dmc_class = dm_converter_class@mp.task_pf(obj, d{1}, mpopt, mpx);
            else
                error('mp.task_cpf.dm_converter_class: d must be 2-element cell array');
            end
        end

        %%-----  data model methods  -----
        function dm_class = data_model_class_default(obj)
            % Implement selector for default data model constructor.

            dm_class = @mp.data_model_cpf;
        end

        function dm = data_model_build(obj, d, dmc, mpopt, mpx)
            % Call superclass :meth:`data_model_build() <mp.task_pf.data_model_build>`
            % for base and target models.

            if iscell(d) && length(d) == 2
                dm  = data_model_build@mp.task_pf(obj, d{1}, dmc, mpopt, mpx);
                dmt = data_model_build@mp.task_pf(obj, d{2}, dmc, mpopt, mpx);
                dm.userdata.target = dmt;
            else
                error('mp.task_cpf.data_model_build: d must be 2-element cell array');
            end
        end

        %%-----  network model methods  -----
        function nm = network_model_build(obj, dm, mpopt, mpx)
            % Call superclass :meth:`network_model_build() <mp.task_pf.network_model_build>`
            % for base and target models.

            dmt = dm.userdata.target;
            nm  = network_model_build@mp.task_pf(obj, dm,  mpopt, mpx);
            nmt = network_model_build@mp.task_pf(obj, dmt, mpopt, mpx);
            nm.userdata.target = nmt;
        end

        function nm = network_model_x_soln(obj, mm, nm)
            % Call superclass :meth:`network_model_x_soln() <mp.task_pf.network_model_x_soln>`
            % then update solution in target network model.

            %% call parent
            nm = network_model_x_soln@mp.task_pf(obj, mm, nm);

            %% update solution in target network model
            nm.userdata.target = mm.network_model_x_soln(nm.userdata.target);
        end

        function nm = network_model_update(obj, mm, nm)
            % Call superclass :meth:`network_model_update() <mp.task_pf.network_model_update>`
            % then update port injection solution by interpolating with
            % parameter lambda.

            %% call parent
            nm = network_model_update@mp.task_pf(obj, mm, nm);

            %% update port injection solution in target network model
            nmt = nm.userdata.target;
            nmt.port_inj_soln();

            %% update port injection solution, by interpolation with lambda
            lambda = mm.soln.x(mm.get_idx('var').i1.lambda);
            nm.soln.gs_ = nm.soln.gs_ + (nmt.soln.gs_ - nm.soln.gs_) * lambda;
        end

        %%-----  mathematical model methods  -----
        function mm_class = math_model_class_default(obj, nm, dm, mpopt)
            % Implement selector for default mathematical model constructor
            % depending on ``mpopt.pf.v_cartesian`` and
            % ``mpopt.pf.current_balance``.

            switch upper(mpopt.model)
                case 'AC'
                    if mpopt.pf.v_cartesian
                        if mpopt.pf.current_balance
                            mm_class = @mp.math_model_cpf_acci;
                        else
                            mm_class = @mp.math_model_cpf_accs;
                        end
                    else
                        if mpopt.pf.current_balance
                            mm_class = @mp.math_model_cpf_acpi;
                        else
                            mm_class = @mp.math_model_cpf_acps;
                        end
                    end
                case 'DC'
                    error('mp.task_cpf.math_model_class_default: CPF not applicable for DC model');
            end
        end

        function opt = math_model_opt(obj, mm, nm, dm, mpopt)
            % Call superclass :meth:`math_model_opt() <mp.task_pf.math_model_opt>`
            % then add warmstart parameters, if available.

            opt = math_model_opt@mp.task_pf(obj, mm, nm, dm, mpopt);

            %% add the warmstart options, if available
            if ~isempty(obj.warmstart)
                opt = mm.solve_opts_warmstart(opt, obj.warmstart, nm);
                obj.warmstart = [];     %% delete warmstart data from task
            end
        end
    end     %% methods
end         %% classdef
