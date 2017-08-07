function N = getN(om, set_type, name, idx)
%GETN  Returns the number of variables, constraints or cost rows.
%   N = OM.GETN(SET_TYPE)
%   N = OM.GETN(SET_TYPE, NAME)
%   N = OM.GETN(SET_TYPE, NAME, IDX)
%
%   Returns either the total number of variables/constraints/cost rows
%   or the number corresponding to a specified named block.
%
%   Examples:
%       N = om.getN('var')      : total number of variables
%       N = om.getN('lin')      : total number of linear constraints
%       N = om.getN('nln')      : total number of nonlinear constraints (legacy)
%       N = om.getN('nle')      : total number of nonlin equality constraints
%       N = om.getN('nli')      : total number of nonlin inequality constraints
%       N = om.getN('cost')         : total number of cost rows (in N)
%       N = om.getN('var', name)    : number of variables in named set
%       N = om.getN('lin', name)    : number of linear constraints in named set
%       N = om.getN('nln', name)    : number of nonlin con in named set (legacy)
%       N = om.getN('nle', name)    : number of nonlin eq cons. in named set
%       N = om.getN('nli', name)    : number of nonlin ineq cons. in named set
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
    N = om.(set_type).N;
else
    if isfield(om.(set_type).idx.N, name)
        if nargin < 4
            idx = {};
        end
        s1 = substruct('.', name, '()', idx);
        N = subsref(om.(set_type).idx.N, s1);
    else
        N = 0;
    end
end
