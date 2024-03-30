function [varargout] = dcopf(varargin)
% dcopf - Solves a DC optimal power flow.
%
% This is a simple wrapper function around opf that sets the ``model``
% option to ``'DC'`` before calling opf.
% See opf for the details of input and output arguments.
%
% See also rundcopf.

%   MATPOWER
%   Copyright (c) 1996-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

[mpc, mpopt] = opf_args(varargin{:});
mpopt = mpoption(mpopt, 'model', 'DC');
[varargout{1:nargout}] = opf(mpc, mpopt);
