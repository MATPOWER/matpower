classdef task_pf_legacy < mp.task_pf & mp.task_shared_legacy
% mp.task_pf_legacy - |MATPOWER| task for legacy power flow (PF).
%
% Adds functionality needed by the *legacy* |/MATPOWER/| *framework* to the
% task implementation for the power flow problem. This consists of
% pre-processing some input data and exporting and packaging result data.
%
% mp.task_pf Methods:
%   * run_pre - pre-process inputs that are for legacy framework only
%   * run_post - export results back to data model source
%   * legacy_post_run - post-process *legacy framework* outputs
%
% See also mp.task_pf, mp.task, mp.task_shared_legacy.

%   MATPOWER
%   Copyright (c) 2020-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        %%-----  task methods  -----
        function [d, mpopt] = run_pre(obj, d, mpopt)
            % Pre-process inputs that are for *legacy framework* only.
            % ::
            %
            %   [d, mpopt] = task.run_pre(d, mpopt)
            %
            % Inputs:
            %   d : data source specification, currently assumed to be a
            %       |MATPOWER| case name or case struct (``mpc``)
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Outputs:
            %   d : updated value of corresponding input
            %   mpopt (struct) : updated value of corresponding input
            %
            % Call :meth:`run_pre_legacy() <mp.task_shared_legacy.run_pre_legacy>`
            % method before calling parent.

            [d, mpopt] = obj.run_pre_legacy(d, mpopt);
            [d, mpopt] = run_pre@mp.task_pf(obj, d, mpopt);
        end

        function obj = run_post(obj, mm, nm, dm, mpopt)
            % Export results back to data model source.
            % ::
            %
            %   task.run_post(mm, nm, dm, mpopt)
            %
            % Inputs:
            %   mm (mp.math_model) : mathmatical model object
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Output:
            %   task (mp.task) : task object
            %
            % Calls mp.dm_converter.export and saves the result
            % in the data model ``source`` property.

            if obj.nm.np ~= 0
                obj.dm.source = obj.dmc.export(obj.dm, obj.dm.source);
            end
        end

        %%-----  other methods  -----
        function [results, success] = legacy_post_run(obj, mpopt)
            % Post-process *legacy framework* outputs.
            % ::
            %
            %   [results, success] = task.legacy_post_run(mpopt)
            %
            % Input:
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Outputs:
            %   results (struct) : results struct for *legacy* |/MATPOWER/|
            %       *framework*, see Table 4.1 in legacy |MUM|.
            %   success (integer) : 1 - succeeded, 0 - failed
            %
            % Extract ``results`` and ``success`` and save the
            % task object in ``results.task`` before returning.

            success = obj.success;
            results = obj.dm.source;
            results.task = obj;
        end
    end     %% methods
end         %% classdef
