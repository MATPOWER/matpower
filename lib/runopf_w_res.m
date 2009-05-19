function [varargout] = runopf_w_res(varargin)
%RUNOPF_W_RES  Runs an optimal power flow with fixed zonal reserves.
%
%   results = runopf_w_res(casename, mpopt, fname, solvedcase)
%   [results, success] = runopf_w_res(casename, mpopt, fname, solvedcase)
%
%   The case file or struct must define a 'reserves' field which is a struct
%   with the following fields:
%       zones   nrz x ng, zone(i, j) = 1, if gen j belongs to zone i
%                                      0, otherwise
%       req     nrz x 1, zonal reserve requirement in MW
%       cost    (ng or ngr) x 1, cost of reserves in $/MW
%       qty     (ng or ngr) x 1, max quantity of reserves in MW (optional)
%   where nrz is the number of reserve zones and ngr is the number of
%   generators belonging to at least one reserve zone and ng is the total
%   number of generators.
%
%   See 'case30_reserves.m' for an example case file with fixed reserves.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2008-2009 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

mpc = loadcase(varargin{1});
mpc = toggle_reserves(mpc, 'on');
[varargout{1:nargout}] = runopf(mpc, varargin{2:nargin});

if isstruct(varargout{1})
	varargout{1} = toggle_reserves(varargout{1}, 'off');
end

return;
