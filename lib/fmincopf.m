function [varargout] = fmincopf(varargin)
% fmincopf - Solves an AC optimal power flow using FMINCON (Opt Tbx 2.x & later).
%
% Uses algorithm 520. Please see opf for the details of input and
% output arguments.

%   MATPOWER
%   Copyright (c) 2000-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

[mpc, mpopt] = opf_args(varargin{:});
mpopt = mpoption(mpopt, 'model', 'AC', 'opf.ac.solver', 'FMINCON');
[varargout{1:nargout}] = opf(mpc, mpopt);
