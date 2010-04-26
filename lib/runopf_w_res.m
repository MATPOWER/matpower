function [varargout] = runopf_w_res(varargin)
%RUNOPF_W_RES  Runs an optimal power flow with fixed zonal reserves.
%   RESULTS = RUNOPF_W_RES(CASEDATA, MPOPT, FNAME, SOLVEDCASE)
%   [RESULTS, SUCCESS] = RUNOPF_W_RES(CASEDATA, MPOPT, FNAME, SOLVEDCASE)
%
%   Runs an optimal power flow with the addition of reserve requirements
%   specified as a set of fixed zonal reserves. See RUNOPF for a
%   description of the input and output arguments, which are the same,
%   with the exception that the case file or struct CASEDATA must define
%   a 'reserves' field, which is a struct with the following fields:
%       zones   nrz x ng, zone(i, j) = 1, if gen j belongs to zone i
%                                      0, otherwise
%       req     nrz x 1, zonal reserve requirement in MW
%       cost    (ng or ngr) x 1, cost of reserves in $/MW
%       qty     (ng or ngr) x 1, max quantity of reserves in MW (optional)
%   where nrz is the number of reserve zones and ngr is the number of
%   generators belonging to at least one reserve zone and ng is the total
%   number of generators.
%
%   In addition to the normal OPF output, the RESULTS struct contains a
%   new 'reserves' field with the following fields, in addition to those
%   provided in the input:
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
%   See T_CASE30_USERFCNS for an example case file with fixed reserves,
%   and TOGGLE_RESERVES for the implementation.
%
%   Calling syntax options:
%       results = runopf_w_res(casedata);
%       results = runopf_w_res(casedata, mpopt);
%       results = runopf_w_res(casedata, mpopt, fname);
%       results = runopf_w_res(casedata, mpopt, fname, solvedcase);
%       [results, success] = runopf_w_res(...);
%
%   Example:
%       results = runopf_w_res('t_case30_userfcns');
%
%   See also RUNOPF, TOGGLE_RESERVES, T_CASE30_USERFCNS.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2008-2010 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

mpc = loadcase(varargin{1});
mpc = toggle_reserves(mpc, 'on');
[varargout{1:nargout}] = runopf(mpc, varargin{2:nargin});

if nargout > 0 && isstruct(varargout{1})
    varargout{1} = toggle_reserves(varargout{1}, 'off');
end
