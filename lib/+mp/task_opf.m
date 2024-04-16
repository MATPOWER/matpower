classdef task_opf < mp.task
% mp.task_opf - |MATPOWER| task for optimal power flow (OPF).
%
% Provides task implementation for the optimal power flow problem.
%
% mp.task_opf Properties:
%   * tag - task tag 'OPF'
%   * name - task name 'Optimal Power Flow'
%   * dc - ``true`` if using DC network model
%
% mp.task_opf Methods:
%   * run_pre - set :attr:`dc` property
%   * print_soln_header - add printout of objective function value
%   * data_model_class_default - select default data model constructor
%   * data_model_build_post - adjust bus voltage limits, if requested
%   * network_model_class_default - select default network model constructor
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
        tag = 'OPF';
        name = 'Optimal Power Flow';
        dc  % ``true`` if using DC network model (from ``mpopt.model``, cached in run_pre())
    end

    methods
        %%-----  task methods  -----
        function [d, mpopt] = run_pre(obj, d, mpopt)
            % Set :attr:`dc` property after calling superclass
            % :meth:`run_pre() <mp.task.run_pre>`, then check for
            % unsupported AC OPF solver selection.
            [d, mpopt] = run_pre@mp.task(obj, d, mpopt);     %% call parent

            %% cache DC model flag
            obj.dc = strcmp(upper(mpopt.model), 'DC');

            %% check for unsupported AC OPF solver selection
            if ~obj.dc
                alg = upper(mpopt.opf.ac.solver);
                switch alg
                    case 'IPOPT'
                        if ~have_feature('ipopt')
                            error('mp.task_opf.run_pre: MPOPT.opf.ac.solver = ''%s'' requires IPOPT (see https://github.com/coin-or/Ipopt)', alg);
                        end
                    case 'FMINCON'
                        if ~have_feature('fmincon')
                            error('mp.task_opf.run_pre: MPOPT.opf.ac.solver = ''%s'' requires FMINCON (Optimization Toolbox 2.x or later)', alg);
                        end
                    case 'KNITRO'
                        if ~have_feature('knitro')
                            error('mp.task_opf.run_pre: MPOPT.opf.ac.solver = ''%s'' requires Artelys Knitro (see https://www.artelys.com/solvers/knitro/)', alg);
                        end
                    case {'MINOPF', 'PDIPM', 'TRALM', 'SDPOPF'}
                        error('mp.task_opf.run_pre: MPOPT.opf.ac.solver = ''%s'' not supported.', alg);
                end
            end
        end

        function print_soln_header(obj, mpopt, fd)
            % Call superclass :meth:`print_soln_header() <mp.task.print_soln_header>`
            % the print out the objective function value.
            if nargin < 3
                fd = 1;     %% print to stdio by default
            end

            print_soln_header@mp.task(obj, mpopt, fd);
            fprintf(fd, ...
                'Objective Function Value = %.2f $/hr\n', obj.mm.soln.f);
        end

        %%-----  data model methods  -----
        function dm_class = data_model_class_default(obj)
            % Implement selector for default data model constructor.

            dm_class = @mp.data_model_opf;
        end

        function dm = data_model_build_post(obj, dm, dmc, mpopt)
            % Call superclass :meth:`data_model_build_post() <mp.task.data_model_build_post>`
            % then adjust bus voltage magnitude limits based on generator
            % ``vm_setpoint``, if requested.

            %% call parent
            dm = data_model_build_post@mp.task(obj, dm, dmc, mpopt);

            if ~obj.dc
                %% if requested, adjust bus voltage magnitude
                %% limits based on generator vm_setpoint
                use_vg = mpopt.opf.use_vg;
                if use_vg
                    dm.set_bus_v_lims_via_vg(use_vg);
                end
            end
        end

        %%-----  network model methods  -----
        function nm_class = network_model_class_default(obj, dm, mpopt)
            % Implement selector for default network model constructor
            % depending on ``mpopt.model`` and ``mpopt.opf.v_cartesian``.
            if obj.dc
                nm_class = @mp.net_model_dc;
            else
                if mpopt.opf.v_cartesian
                    nm_class = @mp.net_model_acc;
                else
                    nm_class = @mp.net_model_acp;
                end
            end
        end

        %%-----  mathematical model methods  -----
        function mm_class = math_model_class_default(obj, nm, dm, mpopt)
            % Implement selector for default mathematical model constructor
            % depending on ``mpopt.model``, ``mpopt.opf.v_cartesian``, and
            % ``mpopt.opf.current_balance``.
            switch upper(mpopt.model)
                case 'AC'
                    if mpopt.opf.v_cartesian
                        if mpopt.opf.current_balance
                            mm_class = @mp.math_model_opf_acci;
                        else
                            mm_class = @mp.math_model_opf_accs;
                        end
                    else
                        if mpopt.opf.current_balance
                            mm_class = @mp.math_model_opf_acpi;
                        else
                            mm_class = @mp.math_model_opf_acps;
                        end
                    end
                case 'DC'
                    mm_class = @mp.math_model_opf_dc;
            end
        end
    end     %% methods
end         %% classdef
