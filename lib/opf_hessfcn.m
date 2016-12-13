function Lxx = opf_hessfcn(x, lambda, cost_mult, om, Ybus, Yf, Yt, mpopt, il)
%OPF_HESSFCN  Evaluates Hessian of Lagrangian for AC OPF.
%   LXX = OPF_HESSFCN(X, LAMBDA, COST_MULT, OM, YBUS, YF, YT, MPOPT, IL)
%
%   Hessian evaluation function for AC optimal power flow, suitable
%   for use with MIPS or FMINCON's interior-point algorithm.
%
%   Inputs:
%     X : optimization vector
%     LAMBDA (struct)
%       .eqnonlin : Lagrange multipliers on power balance equations
%       .ineqnonlin : Kuhn-Tucker multipliers on constrained branch flows
%     COST_MULT : (optional) Scale factor to be applied to the cost
%          (default = 1).
%     OM : OPF model object
%     YBUS : bus admittance matrix
%     YF : admittance matrix for "from" end of constrained branches
%     YT : admittance matrix for "to" end of constrained branches
%     MPOPT : MATPOWER options struct
%     IL : (optional) vector of branch indices corresponding to
%          branches with flow limits (all others are assumed to be
%          unconstrained). The default is [1:nl] (all branches).
%          YF and YT contain only the rows corresponding to IL.
%
%   Outputs:
%     LXX : Hessian of the Lagrangian.
%
%   Examples:
%       Lxx = opf_hessfcn(x, lambda, cost_mult, om, Ybus, Yf, Yt, mpopt);
%       Lxx = opf_hessfcn(x, lambda, cost_mult, om, Ybus, Yf, Yt, mpopt, il);
%
%   See also OPF_COSTFCN, OPF_CONSFCN.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
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

%% default args
if isempty(cost_mult)
    cost_mult = 1;
end

%% unpack data
mpc = get_mpc(om);
[baseMVA, bus, gen, branch, gencost] = ...
    deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.gencost);
cp = get_cost_params(om);
[N, Cw, H, dd, rh, kk, mm] = deal(cp.N, cp.Cw, cp.H, cp.dd, ...
                                    cp.rh, cp.kk, cp.mm);
vv = get_idx(om);

%% unpack needed parameters
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of branches
ng = size(gen, 1);          %% number of dispatchable injections
nxyz = length(x);           %% total number of control vars of all types

%% set default constrained lines
if nargin < 8
    il = (1:nl);            %% all lines have limits by default
end
nl2 = length(il);           %% number of constrained lines

%% grab Pg & Qg
Pg = x(vv.i1.Pg:vv.iN.Pg);  %% active generation in p.u.
Qg = x(vv.i1.Qg:vv.iN.Qg);  %% reactive generation in p.u.

%% put Pg & Qg back in gen
gen(:, PG) = Pg * baseMVA;  %% active generation in MW
gen(:, QG) = Qg * baseMVA;  %% reactive generation in MVAr
 
%% reconstruct V
Va = zeros(nb, 1);
Va = x(vv.i1.Va:vv.iN.Va);
Vm = x(vv.i1.Vm:vv.iN.Vm);
V = Vm .* exp(1j * Va);
nxtra = nxyz - 2*nb;
pcost = gencost(1:ng, :);
if size(gencost, 1) > ng
    qcost = gencost(ng+1:2*ng, :);
else
    qcost = [];
end

%% ----- evaluate d2f -----
d2f_dPg2 = sparse(ng, 1);               %% w.r.t. p.u. Pg
d2f_dQg2 = sparse(ng, 1);               %% w.r.t. p.u. Qg
ipolp = find(pcost(:, MODEL) == POLYNOMIAL);
d2f_dPg2(ipolp) = baseMVA^2 * polycost(pcost(ipolp, :), Pg(ipolp)*baseMVA, 2);
if ~isempty(qcost)          %% Qg is not free
    ipolq = find(qcost(:, MODEL) == POLYNOMIAL);
    d2f_dQg2(ipolq) = baseMVA^2 * polycost(qcost(ipolq, :), Qg(ipolq)*baseMVA, 2);
end
i = [vv.i1.Pg:vv.iN.Pg vv.i1.Qg:vv.iN.Qg]';
d2f = sparse(i, i, [d2f_dPg2; d2f_dQg2], nxyz, nxyz);

%% generalized cost
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
    HwC = H * w + Cw;
    AA = N' * M * (LL + 2 * QQ * diagrr);
    d2f = d2f + AA * H * AA' + 2 * N' * M * QQ * sparse(1:nw, 1:nw, HwC, nw, nw) * N;
end
d2f = d2f * cost_mult;

%%----- evaluate Hessian of power balance constraints -----
nlam = length(lambda.eqnonlin) / 2;
lamP = lambda.eqnonlin(1:nlam);
lamQ = lambda.eqnonlin((1:nlam)+nlam);
[Gpaa, Gpav, Gpva, Gpvv] = d2Sbus_dV2(Ybus, V, lamP);
[Gqaa, Gqav, Gqva, Gqvv] = d2Sbus_dV2(Ybus, V, lamQ);
%% constant impedance part of ZIP loads
diaglam = sparse(1:nb, 1:nb, lamP, nb, nb);
Sd = makeSdzip(baseMVA, bus, mpopt);
diagSdz = sparse(1:nb, 1:nb, Sd.z, nb, nb);
Gpvv = Gpvv + 2 * diaglam * diagSdz;
d2G = [
    real([Gpaa Gpav; Gpva Gpvv]) + imag([Gqaa Gqav; Gqva Gqvv]) sparse(2*nb, nxtra);
    sparse(nxtra, 2*nb + nxtra)
];

%%----- evaluate Hessian of flow constraints -----
nmu = length(lambda.ineqnonlin) / 2;
if nmu
    muF = lambda.ineqnonlin(1:nmu);
    muT = lambda.ineqnonlin((1:nmu)+nmu);
else    %% keep dimensions of empty matrices/vectors compatible
    muF = zeros(0,1);   %% (required to avoid problems when using Knitro
    muT = zeros(0,1);   %%  on cases with all lines unconstrained)
end
if upper(mpopt.opf.flow_lim(1)) == 'I'      %% current
    [dIf_dVa, dIf_dVm, dIt_dVa, dIt_dVm, If, It] = dIbr_dV(branch(il,:), Yf, Yt, V);
    [Hfaa, Hfav, Hfva, Hfvv] = d2AIbr_dV2(dIf_dVa, dIf_dVm, If, Yf, V, muF);
    [Htaa, Htav, Htva, Htvv] = d2AIbr_dV2(dIt_dVa, dIt_dVm, It, Yt, V, muT);
else
  f = branch(il, F_BUS);    %% list of "from" buses
  t = branch(il, T_BUS);    %% list of "to" buses
  Cf = sparse(1:nl2, f, ones(nl2, 1), nl2, nb);     %% connection matrix for line & from buses
  Ct = sparse(1:nl2, t, ones(nl2, 1), nl2, nb);     %% connection matrix for line & to buses
  [dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St] = dSbr_dV(branch(il,:), Yf, Yt, V);
  if upper(mpopt.opf.flow_lim(1)) == 'P'    %% real power
    [Hfaa, Hfav, Hfva, Hfvv] = d2ASbr_dV2(real(dSf_dVa), real(dSf_dVm), real(Sf), Cf, Yf, V, muF);
    [Htaa, Htav, Htva, Htvv] = d2ASbr_dV2(real(dSt_dVa), real(dSt_dVm), real(St), Ct, Yt, V, muT);
  else                                      %% apparent power
    [Hfaa, Hfav, Hfva, Hfvv] = d2ASbr_dV2(dSf_dVa, dSf_dVm, Sf, Cf, Yf, V, muF);
    [Htaa, Htav, Htva, Htvv] = d2ASbr_dV2(dSt_dVa, dSt_dVm, St, Ct, Yt, V, muT);
  end
end
d2H = [
    [Hfaa Hfav; Hfva Hfvv] + [Htaa Htav; Htva Htvv] sparse(2*nb, nxtra);
    sparse(nxtra, 2*nb + nxtra)
];

%%-----  do numerical check using (central) finite differences  -----
if 0
    nx = length(x);
    step = 1e-5;
    num_d2f = sparse(nx, nx);
    num_d2G = sparse(nx, nx);
    num_d2H = sparse(nx, nx);
    for i = 1:nx
        xp = x;
        xm = x;
        xp(i) = x(i) + step/2;
        xm(i) = x(i) - step/2;
        % evaluate cost & gradients
        [fp, dfp] = opf_costfcn(xp, om);
        [fm, dfm] = opf_costfcn(xm, om);
        % evaluate constraints & gradients
        [Hp, Gp, dHp, dGp] = opf_consfcn(xp, om, Ybus, Yf, Yt, mpopt, il);
        [Hm, Gm, dHm, dGm] = opf_consfcn(xm, om, Ybus, Yf, Yt, mpopt, il);
        num_d2f(:, i) = cost_mult * (dfp - dfm) / step;
        num_d2G(:, i) = (dGp - dGm) * lambda.eqnonlin   / step;
        num_d2H(:, i) = (dHp - dHm) * lambda.ineqnonlin / step;
    end
    d2f_err = full(max(max(abs(d2f - num_d2f))));
    d2G_err = full(max(max(abs(d2G - num_d2G))));
    d2H_err = full(max(max(abs(d2H - num_d2H))));
    if d2f_err > 1e-6
        fprintf('Max difference in d2f: %g\n', d2f_err);
    end
    if d2G_err > 1e-5
        fprintf('Max difference in d2G: %g\n', d2G_err);
    end
    if d2H_err > 1e-6
        fprintf('Max difference in d2H: %g\n', d2H_err);
    end
end

Lxx = d2f + d2G + d2H;
