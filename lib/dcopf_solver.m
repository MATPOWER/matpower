function [results, success, raw] = dcopf_solver(om, mpopt)
%DCOPF_SOLVER  Solves a DC optimal power flow.
%
%   [RESULTS, SUCCESS, RAW] = DCOPF_SOLVER(OM, MPOPT)
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
%   See also OPF, QPS_MATPOWER.

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
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%% options
alg = upper(mpopt.opf.dc.solver);

%% unpack data
mpc = get_mpc(om);
[baseMVA, bus, gen, branch, gencost] = ...
    deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.gencost);
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
ny = getN(om, 'var', 'y');  %% number of piece-wise linear costs
nxyz = getN(om, 'var');     %% total number of control vars of all types

%% linear constraints & variable bounds
[A, l, u] = linear_constraints(om);
[x0, xmin, xmax] = getv(om);

%% set up objective function of the form: f = 1/2 * X'*HH*X + CC'*X
%% where X = [x;y;z]. First set up as quadratic function of w,
%% f = 1/2 * w'*HHw*w + CCw'*w, where w = diag(M) * (N*X - Rhat). We
%% will be building on the (optionally present) user supplied parameters.

%% piece-wise linear costs
any_pwl = (ny > 0);
if any_pwl
    Npwl = sparse(ones(ny,1), vv.i1.y:vv.iN.y, 1, 1, nxyz);     %% sum of y vars
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

%% default solver
if strcmp(alg, 'DEFAULT')
    if have_fcn('cplex')        %% use CPLEX by default, if available
        alg = 'CPLEX';
    elseif have_fcn('gurobi')   %% if not, then Gurobi, if available
        alg = 'GUROBI';
    elseif have_fcn('mosek')    %% if not, then MOSEK, if available
        alg = 'MOSEK';
    elseif have_fcn('bpmpd')    %% if not, then BPMPD_MEX, if available
        alg = 'BPMPD';
    elseif have_fcn('quadprog') %% if not, then Optimization Tbx, if available
        alg = 'OT';
    elseif (isempty(HH) || ~any(any(HH))) && have_fcn('glpk') %% if not, and
        alg = 'GLPK';           %% prob is LP (not QP), then GLPK, if available
    else                        %% otherwise MIPS
        alg = 'MIPS';
    end
end

%% set up input for QP solver
opt = struct('alg', alg, 'verbose', mpopt.verbose);
switch alg
    case 'MIPS'
        %% try to select an interior initial point
        Varefs = bus(bus(:, BUS_TYPE) == REF, VA) * (pi/180);

        lb = xmin; ub = xmax;
        lb(xmin == -Inf) = -1e10;   %% replace Inf with numerical proxies
        ub(xmax ==  Inf) =  1e10;
        x0 = (lb + ub) / 2;         %% set x0 mid-way between bounds
        k = find(xmin == -Inf & xmax < Inf);    %% if only bounded above
        x0(k) = xmax(k) - 1;                    %% set just below upper bound
        k = find(xmin > -Inf & xmax == Inf);    %% if only bounded below
        x0(k) = xmin(k) + 1;                    %% set just above lower bound
        x0(vv.i1.Va:vv.iN.Va) = Varefs(1);  %% angles set to first reference angle
        if ny > 0
            ipwl = find(gencost(:, MODEL) == PW_LINEAR);
            c = gencost(sub2ind(size(gencost), ipwl, NCOST+2*gencost(ipwl, NCOST)));    %% largest y-value in CCV data
            x0(vv.i1.y:vv.iN.y) = max(c) + 0.1 * abs(max(c));
        end

        %% set up options
        feastol = mpopt.mips.feastol;
        if feastol == 0
            feastol = mpopt.opf.violation;  %% = MPOPT.opf.violation by default
        end
        opt.mips_opt = struct(  'feastol', feastol, ...
                                'gradtol', mpopt.mips.gradtol, ...
                                'comptol', mpopt.mips.comptol, ...
                                'costtol', mpopt.mips.costtol, ...
                                'max_it', mpopt.mips.max_it, ...
                                'step_control', mpopt.mips.step_control, ...
                                'max_red', mpopt.mips.sc.red_it, ...
                                'cost_mult', 1  );
    case 'CPLEX'
        opt.cplex_opt = cplex_options([], mpopt);
    case 'GLPK'
        opt.glpk_opt = glpk_options([], mpopt);
    case 'GUROBI'
        opt.grb_opt = gurobi_options([], mpopt);
    case 'IPOPT'
        opt.ipopt_opt = ipopt_options([], mpopt);
    case 'MOSEK'
        opt.mosek_opt = mosek_options([], mpopt);
    case 'OT'
        if isfield(mpopt, 'linprog') && ~isempty(mpopt.linprog)
            opt.linprog_opt = mpopt.linprog;
        end
        if isfield(mpopt, 'quadprog') && ~isempty(mpopt.quadprog)
            opt.quadprog_opt = mpopt.quadprog;
        end
end

%%-----  run opf  -----
[x, f, info, output, lambda] = qps_matpower(HH, CC, A, l, u, xmin, xmax, x0, opt);
success = (info == 1);

%%-----  calculate return values  -----
if ~any(isnan(x))
    %% update solution data
    Va = x(vv.i1.Va:vv.iN.Va);
    Pg = x(vv.i1.Pg:vv.iN.Pg);
    f = f + C0;

    %% update voltages & generator outputs
    bus(:, VM) = ones(nb, 1);
    bus(:, VA) = Va * 180/pi;
    gen(:, PG) = Pg * baseMVA;

    %% compute branch flows
    branch(:, [QF, QT]) = zeros(nl, 2);
    branch(:, PF) = (Bf * Va + Pfinj) * baseMVA;
    branch(:, PT) = -branch(:, PF);
end

%% package up results
mu_l = lambda.mu_l;
mu_u = lambda.mu_u;
muLB = lambda.lower;
muUB = lambda.upper;

%% update Lagrange multipliers
il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
bus(:, [LAM_P, LAM_Q, MU_VMIN, MU_VMAX]) = zeros(nb, 4);
gen(:, [MU_PMIN, MU_PMAX, MU_QMIN, MU_QMAX]) = zeros(size(gen, 1), 4);
branch(:, [MU_SF, MU_ST]) = zeros(nl, 2);
bus(:, LAM_P)       = (mu_u(ll.i1.Pmis:ll.iN.Pmis) - mu_l(ll.i1.Pmis:ll.iN.Pmis)) / baseMVA;
branch(il, MU_SF)   = mu_u(ll.i1.Pf:ll.iN.Pf) / baseMVA;
branch(il, MU_ST)   = mu_u(ll.i1.Pt:ll.iN.Pt) / baseMVA;
gen(:, MU_PMIN)     = muLB(vv.i1.Pg:vv.iN.Pg) / baseMVA;
gen(:, MU_PMAX)     = muUB(vv.i1.Pg:vv.iN.Pg) / baseMVA;
pimul = [
  mu_l - mu_u;
 -ones(ny>0, 1);    %% dummy entry corresponding to linear cost row in A (in MINOS)
  muLB - muUB
];

mu = struct( ...
  'var', struct('l', muLB, 'u', muUB), ...
  'lin', struct('l', mu_l, 'u', mu_u) );

results = mpc;
[results.bus, results.branch, results.gen, ...
    results.om, results.x, results.mu, results.f] = ...
        deal(bus, branch, gen, om, x, mu, f);

raw = struct('xr', x, 'pimul', pimul, 'info', info, 'output', output);
