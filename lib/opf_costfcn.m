function [f, df, d2f] = opf_costfcn(x, om)
%OPF_COSTFCN  Evaluates objective function, gradient and Hessian for OPF.
%   [F, DF, D2F] = OPF_COSTFCN(X, OM)
%
%   Objective function evaluation routine for AC optimal power flow,
%   suitable for use with MIPS or FMINCON. Computes objective function value,
%   gradient and Hessian.
%
%   Inputs:
%     X : optimization vector
%     OM : OPF model object
%
%   Outputs:
%     F   : value of objective function
%     DF  : (optional) gradient of objective function (column vector)
%     D2F : (optional) Hessian of objective function (sparse matrix)
%
%   Examples:
%       f = opf_costfcn(x, om);
%       [f, df] = opf_costfcn(x, om);
%       [f, df, d2f] = opf_costfcn(x, om);
%
%   See also OPF_CONSFCN, OPF_HESSFCN.

%   MATPOWER
%   Copyright (c) 1996-2018, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

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
