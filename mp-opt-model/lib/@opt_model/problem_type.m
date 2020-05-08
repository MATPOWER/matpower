function prob = problem_type(om)
%PROBLEM_TYPE  Return a string identifying the type of mathematical program
%   PROB = OM.PROBLEM_TYPE()
%
%   Outputs:
%       PROB : problem type, one of the following strings:
%           'LP'    - linear program
%           'QP'    - quadratic program
%           'NLP'   - nonlinear program
%           'MILP'  - mixed-integer linear program
%           'MIQP'  - mixed-integer quadratic program
%           'MINLP' - mixed-integer nonlinear program
%
%   See also OPT_MODEL

%   MP-Opt-Model
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if om.getN('nle') || om.getN('nli') || om.getN('nlc')
    prob = 'NLP';
else
    %% get cost parameters
    [H, C] = om.params_quad_cost();
    if isempty(H) || ~any(any(H))
        prob = 'LP';
    else
        prob = 'QP';
    end
end
if om.is_mixed_integer()
    prob = ['MI' prob];
end
