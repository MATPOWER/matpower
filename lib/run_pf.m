function varargout = run_pf(varargin)
% run_pf - Run a power flow.
% ::
%
%   run_pf(d, mpopt)
%   run_pf(d, mpopt, ...)
%   task = run_pf(...)
%
% This is the main function used to run power flow (PF) problems via the
% **flexible** |*MATPOWER*| **framework**.
%
% This function is a simple wrapper around run_mp, calling it
% with the first argument set to ``@mp.task_pf``.
%
% Inputs:
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
%   task (mp.task_pf) : task object containing the solved run including the
%       data, network, and mathematical model objects.
%
% Solution results are available in the data model, and its elements,
% contained in the returned task object. For example:
% ::
%
%   task = run_pf('case9');
%   va = task.dm.elements.bus.tab.va    % bus voltage angles
%   pg = task.dm.elements.gen.tab.pg    % generator active dispatch
%
% See also run_mp, mp.task_pf.

%   MATPOWER
%   Copyright (c) 2021-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

[varargout{1:nargout}] = run_mp(@mp.task_pf, varargin{:});
