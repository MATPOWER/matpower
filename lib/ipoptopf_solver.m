function [results, success, raw] = ipoptopf_solver(om, mpopt)
%IPOPTOPF_SOLVER  Solves AC optimal power flow using MIPS.
%
%   [RESULTS, SUCCESS, RAW] = IPOPTOPF_SOLVER(OM, MPOPT)
%
%   Inputs are an OPF model object and a MATPOWER options struct.
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
%           .nln    (deprecated) 2*nb+2*nl - Pmis, Qmis, Sf, St
%               .l  lower bounds on nonlinear constraints
%               .u  upper bounds on nonlinear constraints
%           .nle    nonlinear equality constraints
%           .nli    nonlinear inequality constraints
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
%   See also OPF, IPOPT.

%   MATPOWER
%   Copyright (c) 2000-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

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
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%% unpack data
mpc = om.get_mpc();
[baseMVA, bus, gen, branch, gencost] = ...
    deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.gencost);
[vv, ll, nne, nni] = om.get_idx();

%% problem dimensions
nb = size(bus, 1);          %% number of buses
ng = size(gen, 1);          %% number of gens
nl = size(branch, 1);       %% number of branches
ny = om.getN('var', 'y');   %% number of piece-wise linear costs

%% linear constraints
[A, l, u] = om.params_lin_constraint();

%% bounds on optimization vars
[x0, xmin, xmax] = om.params_var();

%% build admittance matrices
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

%% try to select an interior initial point, unless requested not to
if mpopt.opf.start < 2
    s = 1;                      %% set init point inside bounds by s
    lb = xmin; ub = xmax;
    lb(xmin == -Inf) = -1e10;   %% replace Inf with numerical proxies
    ub(xmax ==  Inf) =  1e10;
    x0 = (lb + ub) / 2;         %% set x0 mid-way between bounds
    k = find(xmin == -Inf & xmax < Inf);    %% if only bounded above
    x0(k) = xmax(k) - s;                    %% set just below upper bound
    k = find(xmin > -Inf & xmax == Inf);    %% if only bounded below
    x0(k) = xmin(k) + s;                    %% set just above lower bound
    Varefs = bus(bus(:, BUS_TYPE) == REF, VA) * (pi/180);
    Vmax = min(bus(:, VMAX), 1.5);
    Vmin = max(bus(:, VMIN), 0.5);
    Vm = (Vmax + Vmin) / 2;
    if mpopt.opf.v_cartesian
        V = Vm * exp(1j*Varefs(1));
        x0(vv.i1.Vr:vv.iN.Vr) = real(V);
        x0(vv.i1.Vi:vv.iN.Vi) = imag(V);
    else
        x0(vv.i1.Va:vv.iN.Va) = Varefs(1);  %% angles set to first reference angle
        x0(vv.i1.Vm:vv.iN.Vm) = Vm;         %% voltage magnitudes
        if ny > 0
            ipwl = find(gencost(:, MODEL) == PW_LINEAR);
            c = gencost(sub2ind(size(gencost), ipwl, NCOST+2*gencost(ipwl, NCOST)));    %% largest y-value in CCV data
            x0(vv.i1.y:vv.iN.y) = max(c) + 0.1 * abs(max(c));
        end
    end
end

%% find branches with flow limits
il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
nl2 = length(il);           %% number of constrained lines
nx = length(x0);

%% replace equality variable bounds with an equality constraint
%% (since IPOPT does not return shadow prices on variables that it eliminates)
kk = find(xmin(nb+1:end) == xmax(nb+1:end));    %% all bounds except ref angles
nk = length(kk);
if nk
    kk = kk + nb;               %% adjust index for missing ref angles
    A = [ A; sparse((1:nk)', kk, 1, nk, nx) ];
    l = [ l; xmin(kk) ];
    u = [ u; xmax(kk) ];
    xmin(kk) = -Inf;
    xmax(kk) = Inf;
end

%%-----  run opf  -----
%% build Jacobian and Hessian structure
nA = size(A, 1);                %% number of original linear constraints
% f = branch(:, F_BUS);                           %% list of "from" buses
% t = branch(:, T_BUS);                           %% list of "to" buses
% Cf = sparse(1:nl, f, ones(nl, 1), nl, nb);      %% connection matrix for line & from buses
% Ct = sparse(1:nl, t, ones(nl, 1), nl, nb);      %% connection matrix for line & to buses
% Cl = Cf + Ct;
% Cb = Cl' * Cl + speye(nb);
% Cl2 = Cl(il, :);
% Cg = sparse(gen(:, GEN_BUS), (1:ng)', 1, nb, ng);
% nz = nx - 2*(nb+ng);
% nxtra = nx - 2*nb;
% Js = [
%     Cb      Cb      Cg              sparse(nb,ng)   sparse(nb,nz);
%     Cb      Cb      sparse(nb,ng)   Cg              sparse(nb,nz);
%     Cl2     Cl2     sparse(nl2, 2*ng)               sparse(nl2,nz);
%     Cl2     Cl2     sparse(nl2, 2*ng)               sparse(nl2,nz);
%     A;
% ];
% [f, df, d2f] = opf_costfcn(x0, om);
% Hs = tril(d2f + [
%     Cb  Cb  sparse(nb,nxtra);
%     Cb  Cb  sparse(nb,nxtra);
%     sparse(nxtra,nx);
% ]);
randx = rand(size(x0));
[h, g, dh, dg] = opf_consfcn(randx, om, Ybus, Yf(il,:), Yt(il,:), mpopt, il);
Js = [dg'; dh'; A];
lam = struct('eqnonlin', rand(size(dg,2),1), 'ineqnonlin', rand(size(dh,2),1) );
Hs = tril(opf_hessfcn(randx, lam, 1, om, Ybus, Yf(il,:), Yt(il,:), mpopt, il));
neq = length(g);
niq = length(h);

%% set options struct for IPOPT
options.ipopt = ipopt_options([], mpopt);

%% extra data to pass to functions
options.auxdata = struct( ...
    'om',       om, ...
    'Ybus',     Ybus, ...
    'Yf',       Yf(il,:), ...
    'Yt',       Yt(il,:), ...
    'mpopt',    mpopt, ...
    'il',       il, ...
    'A',        A, ...
    'nA',       nA, ...
    'neqnln',   neq, ...
    'niqnln',   niq, ...
    'Js',       Js, ...
    'Hs',       Hs    );

% %% check Jacobian and Hessian structure
% xr                  = rand(size(x0));
% lambda              = rand(neq+niq, 1);
% options.auxdata.Js  = jacobian(xr, options.auxdata);
% options.auxdata.Hs  = tril(hessian(xr, 1, lambda, options.auxdata));
% Js1 = options.auxdata.Js;
% options.auxdata.Js = Js;
% Hs1 = options.auxdata.Hs;
% [i1, j1, s] = find(Js);
% [i2, j2, s] = find(Js1);
% if length(i1) ~= length(i2) || norm(i1-i2) ~= 0 || norm(j1-j2) ~= 0
%     error('something''s wrong with the Jacobian structure');
% end
% [i1, j1, s] = find(Hs);
% [i2, j2, s] = find(Hs1);
% if length(i1) ~= length(i2) || norm(i1-i2) ~= 0 || norm(j1-j2) ~= 0
%     error('something''s wrong with the Hessian structure');
% end

%% define variable and constraint bounds
options.lb = xmin;
options.ub = xmax;
options.cl = [zeros(neq, 1);  -Inf(niq, 1); l];
options.cu = [zeros(neq, 1); zeros(niq, 1); u];

%% assign function handles
funcs.objective         = @objective;
funcs.gradient          = @gradient;
funcs.constraints       = @constraints;
funcs.jacobian          = @jacobian;
funcs.hessian           = @hessian;
funcs.jacobianstructure = @(d) Js;
funcs.hessianstructure  = @(d) Hs;
%funcs.jacobianstructure = @jacobianstructure;
%funcs.hessianstructure  = @hessianstructure;

%% run the optimization
if have_fcn('ipopt_auxdata')
    [x, info] = ipopt_auxdata(x0,funcs,options);
else
    [x, info] = ipopt(x0,funcs,options);
end

if info.status == 0 || info.status == 1
    success = 1;
else
    success = 0;
end
if isfield(info, 'iter')
    output.iterations = info.iter;
else
    output.iterations = [];
end
f = opf_costfcn(x, om);

%% update solution data
if mpopt.opf.v_cartesian
    Vi = x(vv.i1.Vi:vv.iN.Vi);
    Vr = x(vv.i1.Vr:vv.iN.Vr);
    V = Vr + 1j*Vi;
    Va = angle(V);
    Vm = abs(V);
else
    Va = x(vv.i1.Va:vv.iN.Va);
    Vm = x(vv.i1.Vm:vv.iN.Vm);
    V = Vm .* exp(1j*Va);
end
Pg = x(vv.i1.Pg:vv.iN.Pg);
Qg = x(vv.i1.Qg:vv.iN.Qg);

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

%% line constraint is typically on square of limit
%% so we must fix multipliers
muSf = zeros(nl, 1);
muSt = zeros(nl, 1);
if ~isempty(il)
    if upper(mpopt.opf.flow_lim(1)) == 'P'
        muSf(il) = info.lambda(om.nle.N+(nni.i1.Sf:nni.iN.Sf));
        muSt(il) = info.lambda(om.nle.N+(nni.i1.St:nni.iN.St));
    else
        muSf(il) = 2 * info.lambda(om.nle.N+(nni.i1.Sf:nni.iN.Sf)) .* branch(il, RATE_A) / baseMVA;
        muSt(il) = 2 * info.lambda(om.nle.N+(nni.i1.St:nni.iN.St)) .* branch(il, RATE_A) / baseMVA;
    end
end

%% extract shadow prices for equality var bounds converted to eq constraints
%% (since IPOPT does not return shadow prices on variables that it eliminates)
if nk
    offset = om.nle.N + om.nli.N + nA - nk;
    lam_tmp = info.lambda(offset+(1:nk));
    kl = find(lam_tmp < 0);             %% lower bound binding
    ku = find(lam_tmp > 0);             %% upper bound binding
    info.zl(kk(kl)) = -lam_tmp(kl);
    info.zu(kk(ku)) =  lam_tmp(ku);

    info.lambda(offset+(1:nk)) = [];    %% remove these shadow prices
    nA = nA - nk;                       %% reduce dimension accordingly
end

%% update Lagrange multipliers
if mpopt.opf.v_cartesian
    if om.userdata.veq
        lam = info.lambda(nne.i1.Veq:nne.iN.Veq);
        mu_Vmax = zeros(size(lam));
        mu_Vmin = zeros(size(lam));
        mu_Vmax(lam > 0) =  lam(lam > 0);
        mu_Vmin(lam < 0) = -lam(lam < 0);
        bus(om.userdata.veq, MU_VMAX) = mu_Vmax;
        bus(om.userdata.veq, MU_VMIN) = mu_Vmin;
    end
    bus(om.userdata.viq, MU_VMAX)  = info.lambda(om.nle.N+(nni.i1.Vmax:nni.iN.Vmax));
    bus(om.userdata.viq, MU_VMIN)  = info.lambda(om.nle.N+(nni.i1.Vmin:nni.iN.Vmin));
else
    bus(:, MU_VMAX)  = info.zu(vv.i1.Vm:vv.iN.Vm);
    bus(:, MU_VMIN)  = info.zl(vv.i1.Vm:vv.iN.Vm);
end
gen(:, MU_PMAX)  = info.zu(vv.i1.Pg:vv.iN.Pg) / baseMVA;
gen(:, MU_PMIN)  = info.zl(vv.i1.Pg:vv.iN.Pg) / baseMVA;
gen(:, MU_QMAX)  = info.zu(vv.i1.Qg:vv.iN.Qg) / baseMVA;
gen(:, MU_QMIN)  = info.zl(vv.i1.Qg:vv.iN.Qg) / baseMVA;
if mpopt.opf.current_balance
    %% convert current balance shadow prices to equivalent lamP and lamQ
    %% P + jQ = (Vr + jVi) * (M - jN)
    %% M = (Vr P + Vi Q) / (Vr^2 + Vi^2)
    %% N = (Vi P - Vr Q) / (Vr^2 + Vi^2)
    %% lamP = df/dP = df/dM * dM/dP + df/dN + dN/dP
    %% lamQ = df/dQ = df/dM * dM/dQ + df/dN + dN/dQ
    VV = V ./ (V .* conj(V));   %% V / Vm^2
    VVr = real(VV);
    VVi = imag(VV);
    lamM = info.lambda(nne.i1.rImis:nne.iN.rImis);
    lamN = info.lambda(nne.i1.iImis:nne.iN.iImis);
    bus(:, LAM_P) = (VVr.*lamM + VVi.*lamN) / baseMVA;
    bus(:, LAM_Q) = (VVi.*lamM - VVr.*lamN) / baseMVA;
else
    bus(:, LAM_P) = info.lambda(nne.i1.Pmis:nne.iN.Pmis) / baseMVA;
    bus(:, LAM_Q) = info.lambda(nne.i1.Qmis:nne.iN.Qmis) / baseMVA;
end
branch(:, MU_SF) = muSf / baseMVA;
branch(:, MU_ST) = muSt / baseMVA;

%% package up results
nlnN = 2*nb + 2*nl;     %% because muSf and muSt are nl x 1, not nl2 x 1

%% extract multipliers for nonlinear constraints
kl = find(info.lambda(1:om.nle.N) < 0);
ku = find(info.lambda(1:om.nle.N) > 0);
nl_mu_l = zeros(nlnN, 1);
nl_mu_u = [zeros(2*nb, 1); muSf; muSt];
nl_mu_l(kl) = -info.lambda(kl);
nl_mu_u(ku) =  info.lambda(ku);

%% extract multipliers for linear constraints
lam_lin = info.lambda(om.nle.N+om.nli.N+(1:nA));    %% lambda for linear constraints
kl = find(lam_lin < 0);                     %% lower bound binding
ku = find(lam_lin > 0);                     %% upper bound binding
mu_l = zeros(nA, 1);
mu_l(kl) = -lam_lin(kl);
mu_u = zeros(nA, 1);
mu_u(ku) = lam_lin(ku);

mu = struct( ...
  'var', struct('l', info.zl, 'u', info.zu), ...
  'nln', struct('l', nl_mu_l, 'u', nl_mu_u), ...
  'nle', info.lambda(1:om.nle.N), ...
  'nli', info.lambda(om.nle.N + (1:om.nli.N)), ...
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
raw = struct('xr', x, 'pimul', pimul, 'info', info.status, 'output', output);


%-----  callback functions  -----
function f = objective(x, d)
f = opf_costfcn(x, d.om);

function df = gradient(x, d)
[f, df] = opf_costfcn(x, d.om);

function c = constraints(x, d)
[hn, gn] = opf_consfcn(x, d.om, d.Ybus, d.Yf, d.Yt, d.mpopt, d.il);
if isempty(d.A)
    c = [gn; hn];
else
    c = [gn; hn; d.A*x];
end

function J = jacobian(x, d)
[hn, gn, dhn, dgn] = opf_consfcn(x, d.om, d.Ybus, d.Yf, d.Yt, d.mpopt, d.il);
J = [dgn'; dhn'; d.A];

function H = hessian(x, sigma, lambda, d)
lam.eqnonlin   = lambda(1:d.neqnln);
lam.ineqnonlin = lambda(d.neqnln+(1:d.niqnln));
H = tril(opf_hessfcn(x, lam, sigma, d.om, d.Ybus, d.Yf, d.Yt, d.mpopt, d.il));

% function Js = jacobianstructure(d)
% Js = d.Js;
% 
% function Hs = hessianstructure(d)
% Hs = d.Hs;