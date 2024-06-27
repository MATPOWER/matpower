function [dm, task_rv] = load_dm(d, task_class, mpopt, varargin)
% mp.load_dm - Load a |MATPOWER| data model.
% ::
%
%   dm = mp.load_dm(d)
%   dm = mp.load_dm(d, task_class)
%   dm = mp.load_dm(d, task_class, mpopt)
%   dm = mp.load_dm(d, task_class, mpopt, ...)
%   [dm, task] = mp.load_dm(...)
%
% Uses a task object to load a |MATPOWER| data model object, optionally
% returning the task object as well.
%
% The resulting data model object can later be passed to run_pf, run_cpf,
% run_opf, or directly to the :meth:`run() <mp.task.run>` method of the
% task object.
%
% Inputs:
%   d : data source specification, currently assumed to be a |MATPOWER|
%       case name or case struct (``mpc``)
%   task_class (function handle) : *(optional)* handle to constructor of
%       task class, *(default is* mp.task_opf *)*
%   mpopt (struct) : *(optional)* |MATPOWER| options struct
%
%       Additional optional inputs can be provided as *<name>, <val>* pairs,
%       with the following options:
%
%       - ``'mpx'`` - |MATPOWER| extension or cell array of |MATPOWER|
%         extensions to apply
%
% Outputs:
%   dm (mp.data_model) : data model object
%   task (mp.task) : task object containing the solved run including the
%       data, network, and mathematical model objects.
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
    if nargin < 2
        task_class = @mp.task_opf;
    end
end
mpx = {};
if rem(length(varargin), 2)
    error('mp.load_dm: arguments following MPOPT must appear in name/value pairs');
end

%% assign overridden inputs
for k = 1:2:length(varargin)
    val  = varargin{k+1};
    switch varargin{k}      %% arg name
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

%% load data model
task.load_dm(d, mpopt, mpx);

dm = task.dm;
if nargout > 1
    task_rv = task;
end
