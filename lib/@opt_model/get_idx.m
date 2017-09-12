function varargout = get_idx(om, varargin)
%GET_IDX  Returns the idx struct for vars, lin/nonlin constraints, costs.
%   VV = OM.GET_IDX()
%   [VV, LL] = OM.GET_IDX()
%   [VV, LL, NNE] = OM.GET_IDX()
%   [VV, LL, NNE, NNI] = OM.GET_IDX()
%   [VV, LL, NNE, NNI, CC] = OM.GET_IDX()
%
%   Returns a structure for each with the beginning and ending
%   index value and the number of elements for each named block.
%   The 'i1' field (that's a one) is a struct with all of the
%   starting indices, 'iN' contains all the ending indices and
%   'N' contains all the sizes. Each is a struct whose fields are
%   the named blocks.
%
%   Alternatively, you can specify the type of named set(s) directly
%   as inputs ...
%
%   [IDX1, IDX2, ...] = OM.GET_IDX(SET_TYPE1, SET_TYPE2, ...);
%   [CC, VV] = OM.GET_IDX('cost', 'var');
%   [LL, NNE, NNI] = OM.GET_IDX('lin', 'nle', 'nli');
%
%   The specific type of named set being referenced is
%   given by the SET_TYPE inputs, with the following valid options:
%       SET_TYPE = 'var'   => variable set
%       SET_TYPE = 'lin'   => linear constraint set
%       SET_TYPE = 'nle'   => nonlinear equality constraint set
%       SET_TYPE = 'nli'   => nonlinear inequality constraint set
%       SET_TYPE = 'cost'  => cost set
%
%   Examples:
%       [vv, ll, nne] = om.get_idx();
%       [vv, ll, cc] = om.get_idx('var', 'lin', 'cost');
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
%       The number of nonlinear equality constraints in a set named 'bar':
%           nbar = nne.N.bar;
%         (note: the following is preferable ...
%           nbar = om.getN('nle', 'bar');
%         ... if you haven't already called get_idx to get nne.)
%
%       If 'z', 'foo' and 'bar' are indexed sets, then you can
%       replace them with something like 'z(i,j)', 'foo(i,j,k)'
%       or 'bar(i)' in the examples above.
%
%   See also OPT_MODEL, ADD_VAR, ADD_LIN_CONSTRAINT, ADD_NLN_CONSTRAINT,
%            ADD_QUAD_COST, ADD_NLN_COST and ADD_LEGACY_COST.

%   MATPOWER
%   Copyright (c) 2008-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin == 1
    varargout{1} = om.var.idx;
    if nargout > 1
        varargout{2} = om.lin.idx;
        if nargout > 2
            varargout{3} = om.nle.idx;
            if nargout > 3
                varargout{4} = om.nli.idx;
                if nargout > 4
                    varargout{5} = om.cost.idx;
                end
            end
        end
    end
else
    for k = nargout:-1:1
        varargout{k} = om.(varargin{k}).idx;
    end
end