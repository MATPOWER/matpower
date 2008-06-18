function [varargout] = dcopf(varargin)
%DCOPF  Solves a DC optimal power flow.
%
%   Please see 'help opf' for the details of input and output arguments.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2008 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

[mpc, mpopt] = opf_args(varargin{:});
mpopt = mpoption(mpopt, 'PF_DC', 1);
[varargout{1:nargout}] = opf(mpc, mpopt);

return;
