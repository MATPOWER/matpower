function [results, success, raw] = dcopf_solver(om, mpopt, output)
%DCOPF_SOLVER  Solves a DC optimal power flow.
%
%   [results, success, raw] = dcopf_solver(om, mpopt)
%   [results, success, raw] = dcopf_solver(om, mpopt, output)
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
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%% options
verbose = mpopt(31);    %% VERBOSE
alg     = mpopt(26);    %% OPF_ALG_DC

if alg == 0
    if have_fcn('bpmpd')
        alg = 100;      %% BPMPD_MEX
    else
        alg = 200;      %% PDIPM (pure Matlab)
    end
end

%% unpack data
mpc = get_mpc(om);
[baseMVA, bus, gen, branch, gencost] = ...
    deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.gencost);
om = build_cost_params(om);
cp = get_cost_params(om);
[N, H, Cw] = deal(cp.N, cp.H, cp.Cw);
fparm = [cp.dd cp.rh cp.kk cp.mm];
[Bf, Pfinj] = deal(mpc.Bf, mpc.Pfinj);
[vv, ll] = get_idx(om);

%% problem dimensions
ipol = find(gencost(:, MODEL) == POLYNOMIAL); %% polynomial costs
ipwl = find(gencost(:, MODEL) == PW_LINEAR);  %% piece-wise linear costs
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of branches
nw = size(N, 1);            %% number of general cost vars, w
ny = get_var_N(om, 'y');    %% number of piece-wise linear costs
nxyz = get_var_N(om);       %% total number of control vars of all types

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
AA  = [ A(ieq, :);  A(ilt, :);  -A(igt, :);  A(ibx, :);  -A(ibx, :) ];
bb  = [ u(ieq);     u(ilt);     -l(igt);     u(ibx);     -l(ibx)    ];

il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
nl2 = length(il);           %% number of constrained lines

%% set up objective function of the form: f = 1/2 * X'*HH*X + CC'*X
%% where X = [x;y;z]. First set up as quadratic function of w,
%% f = 1/2 * w'*HHw*w + CCw'*w, where w = diag(M) * (N*X - Rhat). We
%% will be building on the (optionally present) user supplied parameters.

%% piece-wise linear costs
any_pwl = (ny > 0);
if any_pwl
    Npwl = sparse(ones(ny,1), vv.i1.y-1+ipwl, 1, 1, nxyz);       %% sum of y vars
    Hpwl = 0;
    Cpwl = 1;
    fparm_pwl = [1 0 0 1];
else
    Npwl = sparse(0, nxyz);
    Hpwl = [];
    Cpwl = [];
    fparm_pwl = [];
end

%% quadratic costs
npol = length(ipol);
if any(find(gencost(ipol, NCOST) > 3))
    error('DC opf cannot handle polynomial costs with higher than quadratic order.');
end
iqdr = find(gencost(ipol, NCOST) == 3);
ilin = find(gencost(ipol, NCOST) == 2);
polycf = zeros(npol, 3);                            %% quadratic coeffs for Pg
if ~isempty(iqdr)
  polycf(iqdr, :)   = gencost(ipol(iqdr), COST:COST+2);
end
polycf(ilin, 2:3) = gencost(ipol(ilin), COST:COST+1);
polycf = polycf * diag([ baseMVA^2 baseMVA 1]);     %% convert to p.u.
Npol = sparse(1:npol, vv.i1.Pg-1+ipol, 1, npol, nxyz);         %% Pg vars
Hpol = sparse(1:npol, 1:npol, 2*polycf(:, 1), npol, npol);
Cpol = polycf(:, 2);
fparm_pol = ones(npol,1) * [ 1 0 0 1 ];

%% combine with user costs
NN = [ Npwl; Npol; N ];
HHw = [ Hpwl, sparse(any_pwl, npol+nw);
        sparse(npol, any_pwl), Hpol, sparse(npol, nw);
        sparse(nw, any_pwl+npol), H   ];
CCw = [Cpwl; Cpol; Cw];
ffparm = [ fparm_pwl; fparm_pol; fparm ];

%% transform quadratic coefficients for w into coefficients for X
nnw = any_pwl+npol+nw;
M   = sparse(1:nnw, 1:nnw, ffparm(:, 4), nnw, nnw);
MR  = M * ffparm(:, 2);
HMR = HHw * MR;
MN  = M * NN;
HH = MN' * HHw * MN;
CC = full(MN' * (CCw - HMR));
C0 = 1/2 * MR' * HMR + sum(polycf(:, 3));   %% constant term of cost

%% run QP solver
mpopt(15) = length(ieq);            %% set number of equality constraints
if ~have_fcn('sparse_qp') || mpopt(51) == 0 %% don't use sparse matrices
    AA = full(AA);
    HH = full(HH);
end

%% bounds on optimization vars
[x0, LB, UB] = getv(om);

%%-----  run opf  -----
if any(any(HH))
  [x, lambda, how, success] = mp_qp(HH, CC, AA, bb, LB, UB, x0, mpopt(15), verbose, alg);
else
  [x, lambda, how, success] = mp_lp(CC, AA, bb, LB, UB, x0, mpopt(15), verbose, alg);
end
info = success;

%% update solution data
Va = x(vv.i1.Va:vv.iN.Va);
Pg = x(vv.i1.Pg:vv.iN.Pg);
f = 1/2 * x' * HH * x + CC' * x + C0;

%%-----  calculate return values  -----
%% update voltages & generator outputs
bus(:, VA) = Va * 180/pi;
gen(:, PG) = Pg * baseMVA;

%% compute branch flows
branch(:, [QF, QT]) = zeros(nl, 2);
branch(:, PF) = (Bf * Va + Pfinj) * baseMVA;
branch(:, PT) = -branch(:, PF);

%% package up results
nA = length(u);
neq = length(ieq);
nlt = length(ilt);
ngt = length(igt);
nbx = length(ibx);

%% extract multipliers
kl = find(lambda(1:neq) < 0);
ku = find(lambda(1:neq) > 0);

mu_l = zeros(nA, 1);
mu_l(ieq) = -lambda(1:neq);
mu_l(ieq(ku)) = 0;
mu_l(igt) = lambda(neq+nlt+(1:ngt));
mu_l(ibx) = lambda(neq+nlt+ngt+nbx+(1:nbx));

mu_u = zeros(nA, 1);
mu_u(ieq) = lambda(1:neq);
mu_u(ieq(kl)) = 0;
mu_u(ilt) = lambda(neq+(1:nlt));
mu_u(ibx) = lambda(neq+nlt+ngt+(1:nbx));

muLB = lambda(nA+(1:nxyz));
muUB = lambda(nA+nxyz+(1:nxyz));

%% update Lagrange multipliers
bus(:, [LAM_P, LAM_Q, MU_VMIN, MU_VMAX]) = zeros(nb, 4);
gen(:, [MU_PMIN, MU_PMAX, MU_QMIN, MU_QMAX]) = zeros(size(gen, 1), 4);
branch(:, [MU_SF, MU_ST]) = zeros(nl, 2);
bus(:, LAM_P)       = (mu_u(ll.i1.Pmis:ll.iN.Pmis) - mu_l(ll.i1.Pmis:ll.iN.Pmis)) / baseMVA;
branch(il, MU_SF)    = mu_u(ll.i1.Pf:ll.iN.Pf) / baseMVA;
branch(il, MU_ST)    = mu_u(ll.i1.Pt:ll.iN.Pt) / baseMVA;
gen(:, MU_PMIN)     = muLB(vv.i1.Pg:vv.iN.Pg) / baseMVA;
gen(:, MU_PMAX)     = muUB(vv.i1.Pg:vv.iN.Pg) / baseMVA;

mu = struct( ...
  'var', struct('l', muLB, 'u', muUB), ...
  'lin', struct('l', mu_l, 'u', mu_u) );

results = struct( ...
  'bus', bus, ...
  'gen', gen, ...
  'branch', branch, ...
  'x', x, ...
  'mu', mu, ...
  'f', f );

%% optional fields
%% 1st one is always computed anyway, just include it
results.dg = A;
if isfield(output, 'g')
  results.g = A * x;
end
if isfield(output, 'df')
  results.df = [];
end
if isfield(output, 'd2f')
  results.d2f = [];
end
pimul = [
  mu_l - mu_u;
 -ones(ny>0, 1);    %% dummy entry corresponding to linear cost row in A (in MINOS)
  muLB - muUB
];
raw = struct('xr', x, 'pimul', pimul, 'info', info);

return;
