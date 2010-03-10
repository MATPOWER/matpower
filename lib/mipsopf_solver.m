function [results, success, raw] = mipsopf_solver(om, mpopt, out_opt)
%MIPSOPF_SOLVER  Solves AC optimal power flow using MIPS.
%
%   [RESULTS, SUCCESS, RAW] = MIPSOPF_SOLVER(OM, MPOPT)
%   [RESULTS, SUCCESS, RAW] = MIPSOPF_SOLVER(OM, MPOPT, OUT_OPT)
%
%   Inputs are an OPF model object, a MATPOWER options vector and
%   a struct containing fields (can be empty) for each of the desired
%   optional output fields.
%
%   Outputs are a RESULTS struct, SUCCESS flag and RAW output struct.
%
%   RESULTS is a MATPOWER case struct (mpc) with the usual baseMVA, bus
%   branch, gen, gencost fields, along with the following additional
%   fields:
%       .order      see 'help ext2int' for details of this field
%       .x          final value of optimization variables (internal order)
%       .f          final objective function value
%       .mu         shadow prices on ...
%           .var
%               .l  lower bounds on variables
%               .u  upper bounds on variables
%           .nln
%               .l  lower bounds on non-linear constraints
%               .u  upper bounds on non-linear constraints
%           .lin
%               .l  lower bounds on linear constraints
%               .u  upper bounds on linear constraints
%       .g          (optional) constraint values
%       .dg         (optional) constraint 1st derivatives
%       .df         (optional) obj fun 1st derivatives (not yet implemented)
%       .d2f        (optional) obj fun 2nd derivatives (not yet implemented)
%
%   SUCCESS     1 if solver converged successfully, 0 otherwise
%
%   RAW         raw output in form returned by MINOS
%       .xr     final value of optimization variables
%       .pimul  constraint multipliers
%       .info   solver specific termination code
%       .output solver specific output information
%
%   See also OPF, MIPS.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Autonoma de Manizales
%   Copyright (c) 2000-2010 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialization -----
%% optional output
if nargin < 3
    out_opt = struct([]);
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
feastol = mpopt(81);    %% PDIPM_FEASTOL
gradtol = mpopt(82);    %% PDIPM_GRADTOL
comptol = mpopt(83);    %% PDIPM_COMPTOL
costtol = mpopt(84);    %% PDIPM_COSTTOL
max_it  = mpopt(85);    %% PDIPM_MAX_IT
max_red = mpopt(86);    %% SCPDIPM_RED_IT
step_control = (mpopt(11) == 565);  %% OPF_ALG == 565, MIPS-sc
if feastol == 0
    feastol = mpopt(16);    %% = OPF_VIOLATION by default
end
opt = struct(   'feastol', feastol, ...,
                'gradtol', gradtol, ...,
                'comptol', comptol, ...,
                'costtol', costtol, ...,
                'max_it', max_it, ...,
                'max_red', max_red, ...,
                'step_control', step_control, ...,
                'cost_mult', 1e-4, ...,
                'verbose', verbose  );

%% unpack data
mpc = get_mpc(om);
[baseMVA, bus, gen, branch, gencost] = ...
    deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.gencost);
[vv, ll, nn] = get_idx(om);

%% problem dimensions
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of branches
ny = getN(om, 'var', 'y');  %% number of piece-wise linear costs

%% linear constraints
[A, l, u] = linear_constraints(om);

%% bounds on optimization vars
[x0, xmin, xmax] = getv(om);

%% build admittance matrices
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

%% try to select an interior initial point
ll = xmin; uu = xmax;
ll(xmin == -Inf) = -1e10;   %% replace Inf with numerical proxies
uu(xmax ==  Inf) =  1e10;
x0 = (ll + uu) / 2;
Varefs = bus(bus(:, BUS_TYPE) == REF, VA) * (pi/180);
x0(vv.i1.Va:vv.iN.Va) = Varefs(1);  %% angles set to first reference angle
if ny > 0
    ipwl = find(gencost(:, MODEL) == PW_LINEAR);
%     PQ = [gen(:, PMAX); gen(:, QMAX)];
%     c = totcost(gencost(ipwl, :), PQ(ipwl));
    c = gencost(sub2ind(size(gencost), ipwl, NCOST+2*gencost(ipwl, NCOST)));    %% largest y-value in CCV data
    x0(vv.i1.y:vv.iN.y) = max(c) + 0.1 * abs(max(c));
%     x0(vv.i1.y:vv.iN.y) = c + 0.1 * abs(c);
end

%% find branches with flow limits
il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
nl2 = length(il);           %% number of constrained lines

%%-----  run opf  -----
f_fcn = @(x)opf_costfcn(x, om);
gh_fcn = @(x)opf_consfcn(x, om, Ybus, Yf(il,:), Yt(il,:), mpopt, il);
hess_fcn = @(x, lambda)opf_hessfcn(x, lambda, om, Ybus, Yf(il,:), Yt(il,:), mpopt, il, opt.cost_mult);
[x, f, info, Output, Lambda] = ...
  mips(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt);
success = (info > 0);

%% update solution data
Va = x(vv.i1.Va:vv.iN.Va);
Vm = x(vv.i1.Vm:vv.iN.Vm);
Pg = x(vv.i1.Pg:vv.iN.Pg);
Qg = x(vv.i1.Qg:vv.iN.Qg);
V = Vm .* exp(1j*Va);

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
muSf = zeros(nl, 1);
muSt = zeros(nl, 1);
if ~isempty(il)
    muSf(il) = 2 * Lambda.ineqnonlin(1:nl2)       .* branch(il, RATE_A) / baseMVA;
    muSt(il) = 2 * Lambda.ineqnonlin((1:nl2)+nl2) .* branch(il, RATE_A) / baseMVA;
end

%% update Lagrange multipliers
bus(:, MU_VMAX)  = Lambda.upper(vv.i1.Vm:vv.iN.Vm);
bus(:, MU_VMIN)  = Lambda.lower(vv.i1.Vm:vv.iN.Vm);
gen(:, MU_PMAX)  = Lambda.upper(vv.i1.Pg:vv.iN.Pg) / baseMVA;
gen(:, MU_PMIN)  = Lambda.lower(vv.i1.Pg:vv.iN.Pg) / baseMVA;
gen(:, MU_QMAX)  = Lambda.upper(vv.i1.Qg:vv.iN.Qg) / baseMVA;
gen(:, MU_QMIN)  = Lambda.lower(vv.i1.Qg:vv.iN.Qg) / baseMVA;
bus(:, LAM_P)    = Lambda.eqnonlin(nn.i1.Pmis:nn.iN.Pmis) / baseMVA;
bus(:, LAM_Q)    = Lambda.eqnonlin(nn.i1.Qmis:nn.iN.Qmis) / baseMVA;
branch(:, MU_SF) = muSf / baseMVA;
branch(:, MU_ST) = muSt / baseMVA;

%% package up results
nlnN = getN(om, 'nln');

%% extract multipliers for non-linear constraints
kl = find(Lambda.eqnonlin < 0);
ku = find(Lambda.eqnonlin > 0);
nl_mu_l = zeros(nlnN, 1);
nl_mu_u = [zeros(2*nb, 1); muSf; muSt];
nl_mu_l(kl) = -Lambda.eqnonlin(kl);
nl_mu_u(ku) =  Lambda.eqnonlin(ku);

mu = struct( ...
  'var', struct('l', Lambda.lower, 'u', Lambda.upper), ...
  'nln', struct('l', nl_mu_l, 'u', nl_mu_u), ...
  'lin', struct('l', Lambda.mu_l, 'u', Lambda.mu_u) );

results = mpc;
[results.bus, results.branch, results.gen, ...
    results.om, results.x, results.mu, results.f] = ...
        deal(bus, branch, gen, om, x, mu, f);

%% optional fields
if isfield(out_opt, 'dg')
  [g, geq, dg, dgeq] = opf_consfcn(x, om, Ybus, Yf, Yt, mpopt);
  results.g = [ geq; g];        %% include this since we computed it anyway
  results.dg = [ dgeq'; dg'];   %% true Jacobian organization
end
if isfield(out_opt, 'g') && isempty(g)
  [g, geq] = opf_consfcn(x, om, Ybus, Yf, Yt, mpopt);
  results.g = [ geq; g];
end
if isfield(out_opt, 'df')
  results.df = [];
end
if isfield(out_opt, 'd2f')
  results.d2f = [];
end
pimul = [ ...
  results.mu.nln.l - results.mu.nln.u;
  results.mu.lin.l - results.mu.lin.u;
  -ones(ny>0, 1);
  results.mu.var.l - results.mu.var.u;
];
raw = struct('xr', x, 'pimul', pimul, 'info', info, 'output', Output);
