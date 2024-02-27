classdef (Abstract) task_shared_legacy < handle
% mp.task_shared_legacy - Shared legacy task functionality.
%
% Provides legacy task functionality shared across different tasks
% (e.g. PF, CPF, OPF), specifically, the pre-processing of input data
% for the experimental system-wide ZIP load data.
%
% mp.task_pf Methods:
%   * run_pre_legacy - handle experimental system-wide ZIP load inputs
%
% See also mp.task.

%   MATPOWER
%   Copyright (c) 2022-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        function [d, mpopt] = run_pre_legacy(obj, d, mpopt)
            % Handle experimental system-wide ZIP load inputs.
            % ::
            %
            %   [d, mpopt] = task.run_pre_legacy(d, mpopt)
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
            % Moves the legacy experimental system-wide ZIP load data from
            % ``mpopt.exp.sys_wide_zip_loads`` to ``d.sys_wide_zip_loads``
            % to make it available to the data model converter
            % (mp.dmce_load_mpc2).
            %
            % Called by :meth:`run_pre() <mp.task.run_pre>`.

            if ~isa(d, 'mp.data_model')
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
