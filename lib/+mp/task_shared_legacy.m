classdef (Abstract) task_shared_legacy < handle

%   MATPOWER
%   Copyright (c) 2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        function [d, mpopt] = run_pre_legacy(obj, d, mpopt)
            if ~isa(d, 'mp.data_model')
                %% Handle experimental system-wide ZIP loads (for backward
                %% compatibility), by moving data from
                %%  mpopt.exp.sys_wide_zip_loads to
                %%  mpc.sys_wide_zip_loads
                if nargin == 3 && ~isempty(mpopt.exp) && ...
                        isfield(mpopt.exp, 'sys_wide_zip_loads') && ...
                        (~isempty(mpopt.exp.sys_wide_zip_loads.pw) || ...
                         ~isempty(mpopt.exp.sys_wide_zip_loads.qw))
                    d.sys_wide_zip_loads = mpopt.exp.sys_wide_zip_loads;
                end
            end
        end
    end     %% methods
end         %% classdef
