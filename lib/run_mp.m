function task_rv = run_mp(task_class, d, mpopt, varargin)
% run_mp - Run any |MATPOWER| simulation.
% ::
%
%   run_mp(task_class, d, mpopt)
%   run_mp(task_class, d, mpopt, ...)
%   task = run_mp(...)
%
% This is **the** main function in the **flexible framework** for running
% |MATPOWER|. It creates the task object, applying any specified extensions,
% runs the task, and prints or saves the solution, if desired.
%
% It is typically called from one of the wrapper functions such as run_pf,
% run_cpf, or run_opf.
%
% Inputs:
%   task_class (function handle) : handle to constructor of default task
%       class for type of task to be run, e.g. mp.task_pf for power
%       flow, mp.task_cpf for CPF, and mp.task_opf for OPF
%   d : data source specification, currently assumed to be a |MATPOWER|
%       case name or case struct (``mpc``)
%   mpopt (struct) : |MATPOWER| options struct
%
%       Additional optional inputs can be provided as *<name>, <val>* pairs,
%       with the following options:
%
%       - ``'print_fname'`` - file name for saving pretty-printed output
%       - ``'soln_fname'`` - file name for saving solved case
%       - ``'mpx'`` - |MATPOWER| extension or cell array of |MATPOWER|
%         extensions to apply
%
% Output:
%   task (mp.task) : task object containing the solved run including the
%       data, network, and mathematical model objects.
%
% Solution results are available in the data model, and its elements,
% contained in the returned task object. For example:
% ::
%
%   task = run_opf('case9');
%   lam_p = task.dm.elements.bus.tab.lam_p  % nodal price
%   pg = task.dm.elements.gen.tab.pg        % generator active dispatch
%
% See also run_pf, run_cpf, run_opf, mp.task.

%   MATPOWER
%   Copyright (c) 2021-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% assign default inputs
if nargin < 3
    mpopt = mpoption;
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
