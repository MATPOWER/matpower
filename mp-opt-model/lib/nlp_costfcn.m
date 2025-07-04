function [f, df, d2f] = nlp_costfcn(mm, x)
% nlp_costfcn - Evaluates objective function, gradient and Hessian for NLP solver.
% ::
%
%   [F, DF, D2F] = NLP_COSTFCN(MM, X)
%
%   Objective function evaluation routine, suitable for use with MIPS,
%   FMINCON, etc. Computes objective function value, gradient and Hessian.
%
%   Inputs:
%     MM : MP-Opt-Model object
%     X : optimization vector
%
%   Outputs:
%     F   : value of objective function
%     DF  : (optional) gradient of objective function (column vector)
%     D2F : (optional) Hessian of objective function (sparse matrix)
%
%   Examples:
%       f = nlp_costfcn(mm, x);
%       [f, df] = nlp_costfcn(mm, x);
%       [f, df, d2f] = nlp_costfcn(mm, x);
%
% See also nlp_consfcn, nlp_hessfcn.

%   MP-Opt-Model
%   Copyright (c) 1996-2025, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%%----- evaluate objective function -----
%% general nonlinear costs
if nargout == 3
    [f, df, d2f]    = mm.nlc.eval(mm.var, x);
    if mm.qdc.NS
        [fq, dfq, d2fq] = mm.qdc.eval(mm.var, x);
        f = f + sum(fq);
        df = df + dfq;
        d2f = d2f + d2fq;
    end
elseif nargout == 2
    [f, df]   = mm.nlc.eval(mm.var, x);
    if mm.qdc.NS
        [fq, dfq] = mm.qdc.eval(mm.var, x);
        f = f + sum(fq);
        df = df + dfq;
    end
else
    f  = mm.nlc.eval(mm.var, x);
    if mm.qdc.NS
        fq = mm.qdc.eval(mm.var, x);
        f = f + sum(fq);
    end
end
