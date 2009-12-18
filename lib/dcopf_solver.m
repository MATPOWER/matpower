function [results, success, raw] = dcopf_solver(om, mpopt, output)
%DCOPF_SOLVER  Solves a DC optimal power flow.
%
%   [results, success, raw] = dcopf_solver(om, mpopt)
%   [results, success, raw] = dcopf_solver(om, mpopt, output)
%
%   Inputs are an OPF model object, a MATPOWER options vector and
%   a struct containing fields (can be empty) for each of the desired
%   optional output fields.
%
%   Outputs are a results struct, success flag and raw output struct.
%
%   results is a MATPOWER case struct (mpc) with the usual baseMVA, bus
%   branch, gen, gencost fields, along with the following additional
%   fields:
%       .order      see 'help ext2int' for details of this field
%       .x          final value of optimization variables (internal order)
%       .f          final objective function value
%       .mu         shadow prices on ...
%           .var
%               .l  lower bounds on variables
%               .u  upper bounds on variables
%           .lin
%               .l  lower bounds on linear constraints
%               .u  upper bounds on linear constraints
%       .g          (optional) constraint values
%       .dg         (optional) constraint 1st derivatives
%       .df         (optional) obj fun 1st derivatives (not yet implemented)
%       .d2f        (optional) obj fun 2nd derivatives (not yet implemented)
%
%   success     1 if solver converged successfully, 0 otherwise
%
%   raw         raw output in form returned by MINOS
%       .xr     final value of optimization variables
%       .pimul  constraint multipliers
%       .info   solver specific termination code

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
    alg = 200;      %% PDIPM (pure Matlab)
end

%% unpack data
mpc = get_mpc(om);
[baseMVA, bus, gen, branch, gencost] = ...
    deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.gencost);
om = build_cost_params(om);
cp = get_cost_params(om);
[N, H, Cw] = deal(cp.N, cp.H, cp.Cw);
fparm = [cp.dd cp.rh cp.kk cp.mm];
Bf = userdata(om, 'Bf');
Pfinj = userdata(om, 'Pfinj');
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
mpopt(15) = length(ieq);    %% set number of equality constraints
if mpopt(51) == 0           %% don't use sparse matrices
    AA = full(AA);
    HH = full(HH);
end

%% bounds on optimization vars
[x0, LB, UB] = getv(om);

if alg == 200 || alg == 250
    %% try to select an interior initial point
    Varefs = bus(bus(:, BUS_TYPE) == REF, VA) * (pi/180);

    lb = LB; ub = UB;
    lb(LB == -Inf) = -1e10;   %% replace Inf with numerical proxies
    ub(UB ==  Inf) =  1e10;
    x0 = (lb + ub) / 2;
    x0(vv.i1.Va:vv.iN.Va) = Varefs(1);  %% angles set to first reference angle
    if ny > 0
        ipwl = find(gencost(:, MODEL) == PW_LINEAR);
        c = gencost(sub2ind(size(gencost), ipwl, NCOST+2*gencost(ipwl, NCOST)));    %% largest y-value in CCV data
        x0(vv.i1.y:vv.iN.y) = max(c) + 0.1 * abs(max(c));
    end

    %% set up options
    feastol = mpopt(81);    %% PDIPM_FEASTOL
    gradtol = mpopt(82);    %% PDIPM_GRADTOL
    comptol = mpopt(83);    %% PDIPM_COMPTOL
    costtol = mpopt(84);    %% PDIPM_COSTTOL
    max_it  = mpopt(85);    %% PDIPM_MAX_IT
    max_red = mpopt(86);    %% SCPDIPM_RED_IT
    if feastol == 0
        feastol = mpopt(16);    %% = OPF_VIOLATION by default
    end
    opt = struct(   'feastol', feastol, ...,
                    'gradtol', gradtol, ...,
                    'comptol', comptol, ...,
                    'costtol', costtol, ...,
                    'max_it', max_it, ...,
                    'max_red', max_red, ...,
                    'cost_mult', 1, ...,
                    'verbose', verbose  );
else
    opt = [];
end

%%-----  run opf  -----
if any(any(HH))
  [x, lambda, how, success] = mp_qp(HH, CC, AA, bb, LB, UB, x0, mpopt(15), verbose, alg, opt);
else
  [x, lambda, how, success] = mp_lp(CC, AA, bb, LB, UB, x0, mpopt(15), verbose, alg, opt);
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

results = mpc;
[results.bus, results.branch, results.gen, ...
    results.om, results.x, results.mu, results.f] = ...
        deal(bus, branch, gen, om, x, mu, f);

%% optional fields
%% 1st one is always computed anyway, just include it
if isfield(output, 'g')
  results.g = A * x;
end
results.dg = A;
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
