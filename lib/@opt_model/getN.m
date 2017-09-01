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
%       N = om.getN('nle')      : total number of nonlin equality constraints
%       N = om.getN('nli')      : total number of nonlin inequality constraints
%       N = om.getN('qdc')      : total number of quadratic cost rows
%       N = om.getN('nlc')      : total number of general nonlinear cost rows
%       N = om.getN('cost')     : total number of legacy cost rows (in N)
%       N = om.getN('var', name)    : # of variables in named set
%       N = om.getN('lin', name)    : # of linear constraints in named set
%       N = om.getN('nle', name)    : # of nonlin eq cons. in named set
%       N = om.getN('nli', name)    : # of nonlin ineq cons. in named set
%       N = om.getN('qdc', name)    : # of quadratic cost rows in named set
%       N = om.getN('nlc', name)    : # of gen nonlin cost rows in named set
%       N = om.getN('cost', name)   : # of legacy cost rows (in N) in named set
%       N = om.getN('var', name, idx) : # of variables in indexed named set
%
%   See also OPT_MODEL.

%   MATPOWER
%   Copyright (c) 2008-2017, Power Systems Engineering Research Center (PSERC)
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
