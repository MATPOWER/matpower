function varargout = run_cpf(varargin)
%run_opf  Run a continuation power flow.
%
%   run_cpf(d, mpopt)
%   run_cpf(d, mpopt, ...)
%   task = run_cpf(...)
%
%   See also run_mp.

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

[varargout{1:nargout}] = run_mp(@mp.task_cpf, varargin{:});
