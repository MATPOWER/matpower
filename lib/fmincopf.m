function [varargout] = fmincopf(varargin)
%FMINCOPF  Solves an AC optimal power flow using FMINCON (Opt Tbx 2.x & later).
%
%   Uses algorithm 520. Please see 'help opf' for the details of input and
%   output arguments.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2000-2008 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

[mpc, mpopt] = opf_args(varargin{:});
mpopt = mpoption(mpopt, 'OPF_ALG', 520);
[varargout{1:nargout}] = opf(mpc, mpopt);
