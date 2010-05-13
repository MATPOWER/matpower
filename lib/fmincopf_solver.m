function [results, success, raw] = fmincopf_solver(om, mpopt)
%FMINCOPF_SOLVER  Solves AC optimal power flow using FMINCON.
%
%   [RESULTS, SUCCESS, RAW] = FMINCOPF_SOLVER(OM, MPOPT)
%
%   Inputs are an OPF model object and a MATPOWER options vector.
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
%
%   SUCCESS     1 if solver converged successfully, 0 otherwise
%
%   RAW         raw output in form returned by MINOS
%       .xr     final value of optimization variables
%       .pimul  constraint multipliers
%       .info   solver specific termination code
%       .output solver specific output information
%
%   See also OPF, FMINCON.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Autonoma de Manizales
%   Copyright (c) 2000-2010 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

%%----- initialization -----
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
[vv, ll, nn] = get_idx(om);

%% problem dimensions
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of branches
ny = getN(om, 'var', 'y');  %% number of piece-wise linear costs

%% bounds on optimization vars
[x0, LB, UB] = getv(om);

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

%% build admittance matrices
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

%% tolerances
if mpopt(19) == 0           %% CONSTR_MAX_IT
  mpopt(19) = 150 + 2*nb;
end

%% find branches with flow limits
il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
nl2 = length(il);           %% number of constrained lines

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
    fmc_hessian = @(x, lambda)opf_hessfcn(x, lambda, 1, om, Ybus, Yf(il,:), Yt(il,:), mpopt, il);
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
  x0 = zeros(getN(om, 'var'), 1);
  x0(vv.i1.Va:vv.iN.Va) = 0;
  x0(vv.i1.Vm:vv.iN.Vm)   = 1;
  x0(vv.i1.Pg:vv.iN.Pg) = (gen(:, PMIN) + gen(:, PMAX)) / 2 / baseMVA;
  x0(vv.i1.Qg:vv.iN.Qg) = (gen(:, QMIN) + gen(:, QMAX)) / 2 / baseMVA;
end

%%-----  run opf  -----
f_fcn = @(x)opf_costfcn(x, om);
gh_fcn = @(x)opf_consfcn(x, om, Ybus, Yf(il,:), Yt(il,:), mpopt, il);
[x, f, info, Output, Lambda] = ...
  fmincon(f_fcn, x0, Af, bf, Afeq, bfeq, LB, UB, gh_fcn, fmoptions);
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
nlt = length(ilt);
ngt = length(igt);
nbx = length(ibx);

%% extract multipliers for non-linear constraints
kl = find(Lambda.eqnonlin < 0);
ku = find(Lambda.eqnonlin > 0);
nl_mu_l = zeros(nlnN, 1);
nl_mu_u = [zeros(2*nb, 1); muSf; muSt];
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

results = mpc;
[results.bus, results.branch, results.gen, ...
    results.om, results.x, results.mu, results.f] = ...
        deal(bus, branch, gen, om, x, mu, f);

pimul = [ ...
  results.mu.nln.l - results.mu.nln.u;
  results.mu.lin.l - results.mu.lin.u;
  -ones(ny>0, 1);
  results.mu.var.l - results.mu.var.u;
];
raw = struct('xr', x, 'pimul', pimul, 'info', info, 'output', Output);
