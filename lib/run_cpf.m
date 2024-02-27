function varargout = run_cpf(varargin)
% run_cpf  Run a continuation power flow.
% ::
%
%   run_cpf(d, mpopt)
%   run_cpf(d, mpopt, ...)
%   task = run_cpf(...)
%
% This is the main function used to run continuation power flow (CPF)
% problems via the **flexible** |*MATPOWER*| **framework**.
%
% This function is a simple wrapper around run_mp, calling it
% with the first argument set to ``@mp.task_cpf``.
%
% Inputs:
%   d : data source specification, currently assumed to be a cell array
%       of two |MATPOWER| case names or case structs (``mpc``), the first
%       being the base case, the second the target case
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
%   task (mp.task_cpf) : task object containing the solved run including the
%       data, network, and mathematical model objects.
%
% Solution results are available in the data model, and its elements,
% contained in the returned task object. For example:
% ::
%
%   task = run_cpf({'case9', 'case9target'});
%   vm = task.dm.elements.bus.tab.vm    % bus voltage magnitudes
%   pg = task.dm.elements.gen.tab.pg    % generator active dispatch
%
% See also run_mp, mp.task_cpf.

%   MATPOWER
%   Copyright (c) 2021-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

[varargout{1:nargout}] = run_mp(@mp.task_cpf, varargin{:});
