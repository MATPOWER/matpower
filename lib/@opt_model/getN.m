function N = getN(om, selector, name, idx)
%GETN  Returns the number of variables, constraints or cost rows.
%   N = OM.GETN(SELECTOR)
%   N = OM.GETN(SELECTOR, NAME)
%   N = OM.GETN(SELECTOR, NAME, IDX)
%
%   Returns either the total number of variables/constraints/cost rows
%   or the number corresponding to a specified named block.
%
%   Examples:
%       N = om.getN('var')          : total number of variables
%       N = om.getN('lin')          : total number of linear constraints
%       N = om.getN('nln')          : total number of nonlinear constraints
%       N = om.getN('cost')         : total number of cost rows (in N)
%       N = om.getN('var', name)    : number of variables in named set
%       N = om.getN('lin', name)    : number of linear constraints in named set
%       N = om.getN('nln', name)    : number of nonlinear cons. in named set
%       N = om.getN('cost', name)   : number of cost rows (in N) in named set
%       N = om.getN('var', name, idx) : number of variables in indexed named set
%
%   See also OPT_MODEL.

%   MATPOWER
%   Copyright (c) 2008-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 3
    N = om.(selector).N;
else
    if isfield(om.(selector).idx.N, name)
        if nargin < 4
            idx = {};
        end
        s1 = substruct('.', name, '()', idx);
        N = subsref(om.(selector).idx.N, s1);
    else
        N = 0;
    end
end
