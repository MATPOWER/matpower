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
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%% unpack data
mpc = get_mpc(om);
[baseMVA, gen, gencost] = deal(mpc.baseMVA, mpc.gen, mpc.gencost);
cp = get_cost_params(om);
[N, Cw, H, dd, rh, kk, mm] = deal(cp.N, cp.Cw, cp.H, cp.dd, ...
                                    cp.rh, cp.kk, cp.mm);
vv = get_idx(om);

%% problem dimensions
ng = size(gen, 1);          %% number of dispatchable injections
ny = getN(om, 'var', 'y');  %% number of piece-wise linear costs
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

%% piecewise linear cost of P and Q
if ny > 0
  ccost = full(sparse(ones(1,ny), vv.i1.y:vv.iN.y, ones(1,ny), 1, nxyz));
  f = f + ccost * x;
else
  ccost = zeros(1, nxyz);
end

%% generalized cost term
if ~isempty(N)
    nw = size(N, 1);
    r = N * x - rh;                 %% Nx - rhat
    iLT = find(r < -kk);            %% below dead zone
    iEQ = find(r == 0 & kk == 0);   %% dead zone doesn't exist
    iGT = find(r > kk);             %% above dead zone
    iND = [iLT; iEQ; iGT];          %% rows that are Not in the Dead region
    iL = find(dd == 1);             %% rows using linear function
    iQ = find(dd == 2);             %% rows using quadratic function
    LL = sparse(iL, iL, 1, nw, nw);
    QQ = sparse(iQ, iQ, 1, nw, nw);
    kbar = sparse(iND, iND, [   ones(length(iLT), 1);
                                zeros(length(iEQ), 1);
                                -ones(length(iGT), 1)], nw, nw) * kk;
    rr = r + kbar;                  %% apply non-dead zone shift
    M = sparse(iND, iND, mm(iND), nw, nw);  %% dead zone or scale
    diagrr = sparse(1:nw, 1:nw, rr, nw, nw);
    
    %% linear rows multiplied by rr(i), quadratic rows by rr(i)^2
    w = M * (LL + QQ * diagrr) * rr;

    f = f + (w' * H * w) / 2 + Cw' * w;
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

  %% piecewise linear cost of P and Q
  df = df + ccost';  % As in MINOS, the linear cost row is additive wrt
                     % any nonlinear cost.

  %% generalized cost term
  if ~isempty(N)
    HwC = H * w + Cw;
    AA = N' * M * (LL + 2 * QQ * diagrr);
    df = df + AA * HwC;
    
    %% numerical check
    if 0    %% 1 to check, 0 to skip check
      ddff = zeros(size(df));
      step = 1e-7;
      tol  = 1e-3;
      for k = 1:length(x)
        xx = x;
        xx(k) = xx(k) + step;
        ddff(k) = (opf_costfcn(xx, om) - f) / step;
      end
      if max(abs(ddff - df)) > tol
        idx = find(abs(ddff - df) == max(abs(ddff - df)));
        fprintf('\nMismatch in gradient\n');
        fprintf('idx             df(num)         df              diff\n');
        fprintf('%4d%16g%16g%16g\n', [ 1:length(df); ddff'; df'; abs(ddff - df)' ]);
        fprintf('MAX\n');
        fprintf('%4d%16g%16g%16g\n', [ idx'; ddff(idx)'; df(idx)'; abs(ddff(idx) - df(idx))' ]);
        fprintf('\n');
      end
    end     %% numeric check
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

    %% generalized cost
    if ~isempty(N)
        d2f = d2f + AA * H * AA' + 2 * N' * M * QQ * sparse(1:nw, 1:nw, HwC, nw, nw) * N;
    end
  end
end
