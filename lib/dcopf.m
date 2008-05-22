function [varargout] = dcopf(varargin)
%DCOPF  Solves a DC optimal power flow.
%
%   Please see 'help opf' for the details of input and output arguments.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2008 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

[baseMVA, bus, gen, branch, areas, gencost, Au, lbu, ubu, mpopt, ...
    N, fparm, H, Cw, z0, zl, zu] = opf_args(varargin{:});

mpopt = mpoption(mpopt, 'PF_DC', 10);
[varargout{1:nargout}] = opf(baseMVA, bus, gen, branch, ...
                areas, gencost, Au, lbu, ubu, mpopt, ...
                N, fparm, H, Cw, z0, zl, zu);

return;
