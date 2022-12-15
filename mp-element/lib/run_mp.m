function task_rv = run_mp(task_class, d, mpopt, varargin)
%run_mp
%
%   run_mp(task_class, d, mpopt)
%   run_mp(task_class, d, mpopt, ...)
%   task = run_mp(...)
%
%   Inputs:
%       task_class - function handle to constructor of default task class for
%           type of task to be run, e.g. @mp.task_pf for power flow,
%           @mp.task_cpf for CPF, and @mp.task_opf for OPF
%       d - input data specification, e.g. MATPOWER case name, case struct, etc.
%       mpopt - MATPOWER options struct
%       additional <name>, <value> pairs, where <name> can be:
%           'print_fname' - name for file to save pretty-printed output to
%           'soln_fname' - name for file to save solved case to
%           'mpx' - MATPOWER extension or cell array of MATPOWER extensions
%               to apply
%   Output:
%       task - task object

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% assign default inputs
if nargin < 3
    mpopt = mpoption
end
print_fname = '';
soln_fname = '';
mpx = {};
if rem(length(varargin), 2)
    error('run_mp: arguments following MPOPT must appear in name/value pairs');
end

%% assign overridden inputs
for k = 1:2:length(varargin)
    val  = varargin{k+1};
    switch varargin{k}      %% arg name
        case 'print_fname'
            print_fname = val;
        case 'solved_case'
            soln_fname = val;
        case 'mpx'
            if iscell(val)
                mpx = val;
            else
                mpx = { val };
            end
    end
end

%% extract extensions from mpopt, if specified
if isfield(mpopt.exp, 'mpx') && ~isempty(mpopt.exp.mpx)
    if iscell(mpopt.exp.mpx)
        mpx = [mpx mpopt.exp.mpx];
    else
        mpx = { mpx{:}, mpopt.exp.mpx };
    end
end

%% apply extensions
for k = 1:length(mpx)
    task_class = mpx{k}.task_class(task_class, mpopt);
end

%% create task object
task = task_class();

%% run task
task.run(d, mpopt, mpx);

%% pretty-print results to console & possibly to file
task.print_soln(mpopt, print_fname);

%% save solved case
if ~isempty(soln_fname) && task.success
    task.save_soln(soln_fname);
end

if nargout
    task_rv = task;
end
