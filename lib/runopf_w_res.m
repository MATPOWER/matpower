function [varargout] = runopf_w_res(varargin)
%RUNOPF_W_RES  Runs an optimal power flow with fixed zonal reserves.
%
%   results = runopf_w_res(casedata, mpopt, fname, solvedcase)
%   [results, success] = runopf_w_res(casedata, mpopt, fname, solvedcase)
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
%	In addition to the normal OPF output, the results struct contains a
%   new reserves field with the following fields:
%       R       - ng x 1, reserves provided by each gen in MW
%       Rmin    - ng x 1, lower limit on reserves provided by each gen, (MW)
%       Rmax    - ng x 1, upper limit on reserves provided by each gen, (MW)
%       mu.l    - ng x 1, shadow price on reserve lower limit, ($/MW)
%       mu.u    - ng x 1, shadow price on reserve upper limit, ($/MW)
%       mu.Pmax - ng x 1, shadow price on Pg + R <= Pmax constraint, ($/MW)
%       prc     - ng x 1, reserve price for each gen equal to maximum of the
%                         shadow prices on the zonal requirement constraint
%                         for each zone the generator belongs to
%
%   See 't_case30_userfcns.m' for an example case file with fixed reserves,
%   and 'toggle_reserves.m' for the implementation.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2008-2010 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

mpc = loadcase(varargin{1});
mpc = toggle_reserves(mpc, 'on');
[varargout{1:nargout}] = runopf(mpc, varargin{2:nargin});

if nargout > 0 && isstruct(varargout{1})
	varargout{1} = toggle_reserves(varargout{1}, 'off');
end
