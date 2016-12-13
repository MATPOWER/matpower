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
%   Copyright (c) 2008-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

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
