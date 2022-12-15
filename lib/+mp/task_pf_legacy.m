classdef task_pf_legacy < mp.task_pf & mp.task_shared_legacy
%MP.TASK_PF_LEGACY  MATPOWER task for legacy power flow (PF).
%   MP.TASK_PF_LEGACY provides implementation for power flow problem.
%
%   Properties
%       ?
%
%   Methods
%       ?
%
%   See also MP.TASK_PF

%   MATPOWER
%   Copyright (c) 2020-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
    end

    methods
        %%-----  task methods  -----
        function [d, mpopt] = run_pre(obj, d, mpopt)
            [d, mpopt] = obj.run_pre_legacy(d, mpopt);
            [d, mpopt] = run_pre@mp.task_pf(obj, d, mpopt);
        end

        function obj = run_post(obj, mm, nm, dm, mpopt);
            if obj.nm.np ~= 0
                obj.dm.source = obj.dmc.export(obj.dm, obj.dm.source);
            end
        end

        %%-----  other methods  -----
        function [results, success] = legacy_post_run(obj, mpopt)
            success = obj.success;
            results = obj.dm.source;
            results.task = obj;
        end
    end     %% methods
end         %% classdef
