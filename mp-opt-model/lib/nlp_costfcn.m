function [f, df, d2f] = nlp_costfcn(om, x)
% nlp_costfcn - Evaluates objective function, gradient and Hessian.
% ::
%
%   [F, DF, D2F] = NLP_COSTFCN(OM, X)
%
%   Objective function evaluation routine, suitable for use with MIPS,
%   FMINCON, etc. Computes objective function value, gradient and Hessian.
%
%   Inputs:
%     OM : Opt-Model object
%     X : optimization vector
%
%   Outputs:
%     F   : value of objective function
%     DF  : (optional) gradient of objective function (column vector)
%     D2F : (optional) Hessian of objective function (sparse matrix)
%
%   Examples:
%       f = nlp_costfcn(om, x);
%       [f, df] = nlp_costfcn(om, x);
%       [f, df, d2f] = nlp_costfcn(om, x);
%
% See also nlp_consfcn, nlp_hessfcn.

%   MP-Opt-Model
%   Copyright (c) 1996-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%%----- evaluate objective function -----
%% general nonlinear costs
if nargout == 3
    [f, df, d2f]    = om.eval_nln_cost(x);
    if om.qdc.NS
        [fq, dfq, d2fq] = om.eval_quad_cost(x);
        f = f + sum(fq);
        df = df + dfq;
        d2f = d2f + d2fq;
    end
elseif nargout == 2
    [f, df]   = om.eval_nln_cost(x);
    if om.qdc.NS
        [fq, dfq] = om.eval_quad_cost(x);
        f = f + sum(fq);
        df = df + dfq;
    end
else
    f  = om.eval_nln_cost(x);
    if om.qdc.NS
        fq = om.eval_quad_cost(x);
        f = f + sum(fq);
    end
end
