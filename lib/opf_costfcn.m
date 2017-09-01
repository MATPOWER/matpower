function [f, df, d2f] = opf_costfcn(x, om, varargin)
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
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialize -----
%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%% unpack data
mpc = om.get_mpc();
[baseMVA, gen, gencost] = deal(mpc.baseMVA, mpc.gen, mpc.gencost);
vv = om.get_idx();

%% problem dimensions
ng = size(gen, 1);          %% number of dispatchable injections
ny = om.getN('var', 'y');   %% number of piece-wise linear costs
nxyz = length(x);           %% total number of control vars of all types

%% grab Pg & Qg
Pg = x(vv.i1.Pg:vv.iN.Pg);  %% active generation in p.u.
Qg = x(vv.i1.Qg:vv.iN.Qg);  %% reactive generation in p.u.

%%----- evaluate objective function -----
%% polynomial cost of P and Q
% use totcost only on polynomial cost; in the minimization problem
% formulation, pwl cost is the sum of the y variables.
ipol = find(gencost(:, MODEL) == POLYNOMIAL);   %% poly MW and MVAr costs
xx = [ Pg; Qg ] * baseMVA;
if ~isempty(ipol)
  f = sum( totcost(gencost(ipol, :), xx(ipol)) );  %% cost of poly P or Q
else
  f = 0;
end

%% general nonlinear costs
[f1, df1, d2f1] = eval_nonlin_cost(om, x);
%norm(f-f1, Inf)
if norm(f-f1, Inf)
    error('Yikes! f: %g\n', norm(f-f1, Inf));
end

%% quadratic costs
if om.qdc.NS
    if nargout == 3
        [fq, dfq, d2fq] = om.eval_quad_cost(x);
    elseif nargout == 2
        [fq, dfq] = om.eval_quad_cost(x);
    else
        fq = om.eval_quad_cost(x);
    end
    f = f + fq;
end

%%----- evaluate cost gradient -----
if nargout > 1
  %% index ranges
  iPg = vv.i1.Pg:vv.iN.Pg;
  iQg = vv.i1.Qg:vv.iN.Qg;

  %% polynomial cost of P and Q
  df_dPgQg = zeros(2*ng, 1);        %% w.r.t p.u. Pg and Qg
  df_dPgQg(ipol) = baseMVA * polycost(gencost(ipol, :), xx(ipol), 1);
  df = zeros(nxyz, 1);
  df(iPg) = df_dPgQg(1:ng);
  df(iQg) = df_dPgQg((1:ng) + ng);
%norm(df-df1, Inf)
if norm(df-df1, Inf)
    error('Yikes! df: %g\n', norm(df-df1, Inf));
end

  %% quadratic costs
  if om.qdc.NS
      df = df + dfq;
  end

  %% ---- evaluate cost Hessian -----
  if nargout > 2
    pcost = gencost(1:ng, :);
    if size(gencost, 1) > ng
        qcost = gencost(ng+1:2*ng, :);
    else
        qcost = [];
    end
    
    %% polynomial generator costs
    d2f_dPg2 = sparse(ng, 1);               %% w.r.t. p.u. Pg
    d2f_dQg2 = sparse(ng, 1);               %% w.r.t. p.u. Qg
    ipolp = find(pcost(:, MODEL) == POLYNOMIAL);
    d2f_dPg2(ipolp) = baseMVA^2 * polycost(pcost(ipolp, :), Pg(ipolp)*baseMVA, 2);
    if ~isempty(qcost)          %% Qg is not free
        ipolq = find(qcost(:, MODEL) == POLYNOMIAL);
        d2f_dQg2(ipolq) = baseMVA^2 * polycost(qcost(ipolq, :), Qg(ipolq)*baseMVA, 2);
    end
    i = [iPg iQg]';
    d2f = sparse(i, i, [d2f_dPg2; d2f_dQg2], nxyz, nxyz);
%norm(d2f-d2f1, Inf)
if norm(d2f-d2f1, Inf)
    error('Yikes! d2f: %g\n', norm(d2f-d2f1, Inf));
end

    %% quadratic costs
    if om.qdc.NS
        d2f = d2f + d2fq;
    end
  end
end
