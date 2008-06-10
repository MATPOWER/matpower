function [results, success, raw] = fmincopf_solver(om, mpopt, output)
%FMINCOPF_SOLVER  Solves an AC optimal power flow using FMINCON (Opt Tbx 2.x & later).
%
%   [results, success, raw] = fmincopf_solver(om, mpopt)
%   [results, success, raw] = fmincopf_solver(om, mpopt, output)
%
%   results
%       .bus
%       .gen
%       .branch
%       .f
%       .var
%       .mu
%           .var
%               .l
%               .u
%           .nln
%               .l
%               .u
%           .lin
%               .l
%               .u
%       .g      (optional)
%       .dg     (optional)
%       .df     (optional)
%       .d2f    (optional)
%   raw
%       .xr
%       .pimul
%       .info

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Autonoma de Manizales
%   Copyright (c) 2000-2008 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialization -----
%% optional output
if nargin < 3
    output = struct([]);
end

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% options
verbose = mpopt(31);    %% VERBOSE

%% unpack data
mpc = get_mpc(om);
[baseMVA, bus, gen, branch] = ...
    deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch);
om = build_cost_params(om);
[vv, ll, nn] = get_idx(om);

%% problem dimensions
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of branches
ny = get_var_N(om, 'y');    %% number of piece-wise linear costs

%% linear constraints
[A, l, u] = linear_constraints(om);

%% so, can we do anything good about lambda initialization?
if all(bus(:, LAM_P) == 0)
  bus(:, LAM_P) = (10)*ones(nb, 1);
end

%% split l <= A*x <= u into less than, equal to, greater than, and
%% doubly-bounded sets
ieq = find( abs(u-l) <= eps );          %% equality
igt = find( u >=  1e10 & l > -1e10 );   %% greater than, unbounded above
ilt = find( l <= -1e10 & u <  1e10 );   %% less than, unbounded below
ibx = find( (abs(u-l) > eps) & (u < 1e10) & (l > -1e10) );
Af  = [ A(ilt, :); -A(igt, :); A(ibx, :); -A(ibx, :) ];
bf  = [ u(ilt);   -l(igt);     u(ibx);    -l(ibx)];
Afeq = A(ieq, :);
bfeq = u(ieq);

%% bounds on optimization vars
[x0, LB, UB] = getv(om);

%% add constraint on ref bus angles
refs = find(bus(:, BUS_TYPE) == REF);
Varefs = bus(refs, VA) * (pi/180);
LB(vv.i1.Va-1+refs) = Varefs;
UB(vv.i1.Va-1+refs) = Varefs;

%% build admittance matrices
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

%% tolerances
if mpopt(19) == 0           %% CONSTR_MAX_IT
  mpopt(19) = 150 + 2*nb;
end

%% basic optimset options needed for fmincon
fmoptions = optimset('GradObj', 'on', 'GradConstr', 'on', ...
            'MaxIter', mpopt(19), 'TolCon', mpopt(16), ...
            'TolX', mpopt(17), 'TolFun', mpopt(18) );
fmoptions.MaxFunEvals = 4 * fmoptions.MaxIter;
if verbose == 0,
  fmoptions.Display = 'off';
elseif verbose == 1
  fmoptions.Display = 'iter';
else
  fmoptions.Display = 'testing';
end

%% select algorithm
otver = ver('optim');
if str2double(otver.Version(1)) < 4
  fmoptions = optimset(fmoptions, 'LargeScale', 'off');
  Af = full(Af);
  Afeq = full(Afeq);
else
  if mpopt(55) == 1           %% active-set
    fmoptions = optimset(fmoptions, 'Algorithm', 'active-set');
    Af = full(Af);
    Afeq = full(Afeq);
  elseif mpopt(55) == 2       %% interior-point, w/ default 'bfgs' Hessian approx
    fmoptions = optimset(fmoptions, 'Algorithm', 'interior-point');
  elseif mpopt(55) == 3       %% interior-point, w/ 'lbfgs' Hessian approx
    fmoptions = optimset(fmoptions, 'Algorithm', 'interior-point', 'Hessian','lbfgs');
  elseif mpopt(55) == 4       %% interior-point, w/ exact user-supplied Hessian
    fmc_hessian = @(x, lambda)hessfmin(x, lambda, om, Ybus, Yf, Yt, mpopt);
    fmoptions = optimset(fmoptions, 'Algorithm', 'interior-point', ...
        'Hessian', 'user-supplied', 'HessFcn', fmc_hessian);
  elseif mpopt(55) == 5       %% interior-point, w/ finite-diff Hessian
    fmoptions = optimset(fmoptions, 'Algorithm', 'interior-point', 'Hessian','fin-diff-grads', 'SubProblem', 'cg');
  else
    error('fmincopf_solver: unknown algorithm specified in FMC_ALG option');
  end
end
% fmoptions = optimset(fmoptions, 'DerivativeCheck', 'on', 'FinDiffType', 'central', 'FunValCheck', 'on');
% fmoptions = optimset(fmoptions, 'Diagnostics', 'on');

%% try to select an interior initial point for interior point solver
if str2double(otver.Version(1)) >= 4 && strcmp(optimget(fmoptions, 'Algorithm'), 'interior-point')
  x0 = zeros(get_var_N(om), 1);
  x0(vv.i1.Va:vv.iN.Va) = 0;
  x0(vv.i1.Vm:vv.iN.Vm)   = 1;
  x0(vv.i1.Pg:vv.iN.Pg) = (gen(:, PMIN) + gen(:, PMAX)) / 2 / baseMVA;
  x0(vv.i1.Qg:vv.iN.Qg) = (gen(:, QMIN) + gen(:, QMAX)) / 2 / baseMVA;
end

%%-----  run opf  -----
fmc_cost = @(x)costfmin(x, om);
fmc_cons = @(x)consfmin(x, om, Ybus, Yf, Yt, mpopt);
[x, f, info, Output, Lambda] = ...
  fmincon(fmc_cost, x0, Af, bf, Afeq, bfeq, LB, UB, fmc_cons, fmoptions);
success = (info > 0);

%% update solution data
Va = x(vv.i1.Va:vv.iN.Va);
Vm = x(vv.i1.Vm:vv.iN.Vm);
Pg = x(vv.i1.Pg:vv.iN.Pg);
Qg = x(vv.i1.Qg:vv.iN.Qg);
V = Vm .* exp(j*Va);

%%-----  calculate return values  -----
%% update voltages & generator outputs
bus(:, VA) = Va * 180/pi;
bus(:, VM) = Vm;
gen(:, PG) = Pg * baseMVA;
gen(:, QG) = Qg * baseMVA;
gen(:, VG) = Vm(gen(:, GEN_BUS));

%% compute branch flows
Sf = V(branch(:, F_BUS)) .* conj(Yf * V);  %% cplx pwr at "from" bus, p.u.
St = V(branch(:, T_BUS)) .* conj(Yt * V);  %% cplx pwr at "to" bus, p.u.
branch(:, PF) = real(Sf) * baseMVA;
branch(:, QF) = imag(Sf) * baseMVA;
branch(:, PT) = real(St) * baseMVA;
branch(:, QT) = imag(St) * baseMVA;

%% line constraint is actually on square of limit
%% so we must fix multipliers
Lambda.ineqnonlin(1:nl)      = 2 * Lambda.ineqnonlin(1:nl)      .* branch(:, RATE_A) / baseMVA;
Lambda.ineqnonlin(nl+1:2*nl) = 2 * Lambda.ineqnonlin(nl+1:2*nl) .* branch(:, RATE_A) / baseMVA;

%% update Lagrange multipliers
bus(:, MU_VMAX)  = Lambda.upper(vv.i1.Vm:vv.iN.Vm);
bus(:, MU_VMIN)  = Lambda.lower(vv.i1.Vm:vv.iN.Vm);
gen(:, MU_PMAX)  = Lambda.upper(vv.i1.Pg:vv.iN.Pg) / baseMVA;
gen(:, MU_PMIN)  = Lambda.lower(vv.i1.Pg:vv.iN.Pg) / baseMVA;
gen(:, MU_QMAX)  = Lambda.upper(vv.i1.Qg:vv.iN.Qg) / baseMVA;
gen(:, MU_QMIN)  = Lambda.lower(vv.i1.Qg:vv.iN.Qg) / baseMVA;
bus(:, LAM_P)    = Lambda.eqnonlin(nn.i1.Pmis:nn.iN.Pmis) / baseMVA;
bus(:, LAM_Q)    = Lambda.eqnonlin(nn.i1.Qmis:nn.iN.Qmis) / baseMVA;
nmis = nn.N.Pmis + nn.N.Qmis;
branch(:, MU_SF) = Lambda.ineqnonlin((nn.i1.Sf:nn.iN.Sf) - nmis) / baseMVA;
branch(:, MU_ST) = Lambda.ineqnonlin((nn.i1.St:nn.iN.St) - nmis) / baseMVA;

%% package up results
nnl = get_nln_N(om);
nlt = length(ilt);
ngt = length(igt);
nbx = length(ibx);

%% extract multipliers for non-linear constraints
kl = find(Lambda.eqnonlin < 0);
ku = find(Lambda.eqnonlin > 0);
nl_mu_l = zeros(nnl, 1);
nl_mu_u = [zeros(2*nb, 1); Lambda.ineqnonlin];
nl_mu_l(kl) = -Lambda.eqnonlin(kl);
nl_mu_u(ku) =  Lambda.eqnonlin(ku);

%% extract multipliers for linear constraints
kl = find(Lambda.eqlin < 0);
ku = find(Lambda.eqlin > 0);

mu_l = zeros(size(u));
mu_l(ieq(kl)) = -Lambda.eqlin(kl);
mu_l(igt) = Lambda.ineqlin(nlt+(1:ngt));
mu_l(ibx) = Lambda.ineqlin(nlt+ngt+nbx+(1:nbx));

mu_u = zeros(size(u));
mu_u(ieq(ku)) = Lambda.eqlin(ku);
mu_u(ilt) = Lambda.ineqlin(1:nlt);
mu_u(ibx) = Lambda.ineqlin(nlt+ngt+(1:nbx));

mu = struct( ...
  'var', struct('l', Lambda.lower, 'u', Lambda.upper), ...
  'nln', struct('l', nl_mu_l, 'u', nl_mu_u), ...
  'lin', struct('l', mu_l, 'u', mu_u) );

results = struct( ...
  'bus', bus, ...
  'gen', gen, ...
  'branch', branch, ...
  'var', x, ...
  'mu', mu, ...
  'f', f );

%% optional fields
if isfield(output, 'dg')
  [g, geq, dg, dgeq] = consfmin(x, om, Ybus, Yf, Yt, mpopt);
  results.dg = [ dgeq'; dg'];   %% true Jacobian organization
  results.g = [ geq; g];        %% include this since we computed it anyway
end
if isfield(output, 'g') && isempty(g)
  [g, geq] = consfmin(x, om, Ybus, Yf, Yt, mpopt);
  results.g = [ geq; g];
end
if isfield(output, 'df')
  results.df = [];
end
if isfield(output, 'd2f')
  results.d2f = [];
end
pimul = [ ...
  results.mu.nln.l - results.mu.nln.u;
  results.mu.lin.l - results.mu.lin.u;
  -ones(ny>0, 1);
  results.mu.var.l - results.mu.var.u;
];
raw = struct('xr', x, 'pimul', pimul, 'info', info);

return;
