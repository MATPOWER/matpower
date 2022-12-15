function varargout = run_opf(varargin)
%run_opf  Run an optimal power flow.
%
%   run_opf(d, mpopt)
%   run_opf(d, mpopt, ...)
%   task = run_opf(...)
%
%   See also run_mp.

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

[varargout{1:nargout}] = run_mp(@mp.task_opf, varargin{:});
