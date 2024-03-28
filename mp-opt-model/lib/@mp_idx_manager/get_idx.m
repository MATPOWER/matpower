function varargout = get_idx(obj, varargin)
% get_idx - Returns the idx struct for the various set types.
% ::
%
%   IDX = OBJ.GET_IDX(SET_TYPE)
%   [IDX1, IDX2, ...] = OBJ.GET_IDX(SET_TYPE1, SET_TYPE2, ...)
%
%   Returns a structure for each set type with the beginning and ending
%   index value and the number of elements for each named block. The 'i1'
%   field (that's a one) is a struct with all of the starting indices, 'iN'
%   contains all the ending indices and 'N' contains all the sizes. Each is
%   a struct whose fields are the named blocks.
%
%   For example, if 'var' is the set type used for a vector x of optimization
%   variables, and 'lin' is for a set of linear constraints, then the
%   following examples illustrate how GET_IDX might be used.
%
%   Examples:
%       [vv, ll] = obj.get_idx('var', 'lin');
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
%       The number of linear constraints in a set named 'bar':
%           nbar = ll.N.bar;
%         (note: the following is preferable ...
%           nbar = obj.getN('lin', 'bar');
%         ... if you haven't already called get_idx to get ll.)
%
%       If 'z', 'foo' and 'bar' are indexed sets, then you can
%       replace them with something like 'z(i,j)', 'foo(i,j,k)'
%       or 'bar(i)' in the examples above.
%
%   The GET_IDX method can be overridden to also return the idx structs of
%   set types in a pre-specified order when called without input arguments.
%
%   E.g.
%       vv = obj.get_idx()          % for variables
%       [vv, ll] = obj.get_idx()    % for linear constraints
%
%   See OPT_MODEL, and its methods GET_IDX, ADD_VAR,
%           ADD_LIN_CONSTRAINT, ADD_NLN_CONSTRAINT, ADD_QUAD_COST and
%           ADD_NLN_COST.

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if nargin ~= 1
    for k = nargout:-1:1
        varargout{k} = obj.(varargin{k}).idx;
    end
end
