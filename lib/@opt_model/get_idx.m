function [vv, ll, nn, cc] = get_idx(om)
%GET_IDX  Returns the idx struct for vars, lin/nln constraints, costs.
%   VV = GET_IDX(OM)
%   [VV, LL] = GET_IDX(OM)
%   [VV, LL, NN] = GET_IDX(OM)
%   [VV, LL, NN, CC] = GET_IDX(OM)
%
%   Returns a structure for each with the beginning and ending
%   index value and the number of elements for each named block.
%   The 'i1' field (that's a one) is a struct with all of the
%   starting indices, 'iN' contains all the ending indices and
%   'N' contains all the sizes. Each is a struct whose fields are
%   the named blocks.
%
%   Examples:
%       [vv, ll, nn] = get_idx(om);
%
%       For a variable block named 'z' we have ...
%           vv.i1.z - starting index for 'z' in optimization vector x
%           vv.iN.z - ending index for 'z' in optimization vector x
%           vv.N.z  - number of elements in 'z'
%
%       To extract a 'z' variable from x:
%           z = x(vv.i1.z:vv.iN.z);
%
%       To extract the multipliers on a linear constraint set
%       named 'foo', where mu_l and mu_u are the full set of
%       linear constraint multipliers:
%           mu_l_foo = mu_l(ll.i1.foo:ll.iN.foo);
%           mu_u_foo = mu_u(ll.i1.foo:ll.iN.foo);
%
%       The number of nonlinear constraints in a set named 'bar':
%           nbar = nn.N.bar;
%         (note: the following is preferable ...
%           nbar = getN(om, 'nln', 'bar');
%         ... if you haven't already called get_idx to get nn.)
%
%       If 'z', 'foo' and 'bar' are indexed sets, then you can
%       replace them with something like 'z(i,j)', 'foo(i,j,k)'
%       or 'bar(i)' in the examples above.
%
%   See also OPT_MODEL, ADD_VARS, ADD_CONSTRAINTS, ADD_COSTS.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2008-2012 by Power System Engineering Research Center (PSERC)
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

vv = om.var.idx;
if nargout > 1
    ll = om.lin.idx;
    if nargout > 2
        nn = om.nln.idx;
        if nargout > 3
            cc = om.cost.idx;
         end
    end
end
