function mso = most_summary(mdo)
% most_summary - Collects and optionally prints a summary of MOST results.
% ::
%
%   ms = most_summary(mdo)
%
% .. note:: Consider this function experimental. It is included because
%   it is often better than nothing, though it is very incomplete.
%
% Input:
%   mdo (struct) : MOST data structure, output (see |MOSTman| for details)
%
% Output:
%   ms (struct) : MOST summary struct with the following fields:
%
%       ======  =======================================================
%        name   description
%       ======  =======================================================
%       f       objective function value
%       nb      :math:`n_b`, number of buses
%       ng      :math:`n_g`, number of generators *(including storage,
%               dispatchable load, etc.)*
%       nl      :math:`n_l`, number of branches
%       nt      :math:`n_t`, number of periods in planning horizon
%       nj_max  :math:`n_j^{max}`, max number of scenarios per period
%       nc_max  :math:`n_c^{max}`, max number of contingencies per scenario
%               in any period
%       Pg      :math:`n_g \times n_t \times n_j^{max} \times (n_c^{max}+1)`,
%               real power generation (MW)
%       Rup     :math:`n_g \times n_t`, upward ramping reserve quantities (MW)
%       Rdn     :math:`n_g \times n_t`, downward ramping reserve quantities
%               (MW)
%       Pf      :math:`n_l \times n_t \times n_j^{max} \times (n_c^{max}+1)`,
%               real power branch flow (MW)
%       u       :math:`n_g \times n_t \times n_j^{max} \times (n_c^{max}+1)`,
%               generator commitment status (0 or 1)
%       lamP    :math:`n_b \times n_t \times n_j^{max} \times (n_c^{max}+1)`,
%               shadow price on power balance ($/MWh)
%       muF     :math:`n_l \times n_t \times n_j^{max} \times (n_c^{max}+1)`,
%               shadow price on flow limits ($/MWh)
%       ======  =======================================================
%
% Printing to the console is currently controlled by the
% ``mdo.qp.opt.verbose`` flag.

%   MOST
%   Copyright (c) 2015-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

tol = 1e-4;
verbose = mdo.QP.opt.verbose;
% verbose = 1;
mpc = mdo.mpc;
nb = size(mpc.bus, 1);
nl = size(mpc.branch, 1);
ng = size(mpc.gen, 1);
nt = mdo.idx.nt;
nj_max = max(mdo.idx.nj);
nc_max = max(max(mdo.idx.nc));
ns = mdo.idx.ns;

%% summarize results
psi = zeros(nt, nj_max, nc_max+1);
Pg = zeros(ng, nt, nj_max, nc_max+1);
Pd = zeros(nb, nt, nj_max, nc_max+1);
if mdo.idx.ntramp
    Rup = mdo.results.Rrp;
    Rdn = mdo.results.Rrm;
else
    Rup = [];
    Rdn = [];
end
u = zeros(ng, nt);
lamP = zeros(nb, nt, nj_max, nc_max+1);
muF = zeros(nl, nt, nj_max, nc_max+1);
Pf = zeros(nl, nt, nj_max, nc_max+1);
for t = 1:nt
  for j = 1:mdo.idx.nj(t)
    for k = 1:mdo.idx.nc(t,j)+1
      rr = mdo.flow(t,j,k).mpc;
      psi(t, j, k) = mdo.CostWeightsAdj(k, j, t);
      u(:, t) = rr.gen(:, GEN_STATUS);
      Pg(:, t, j, k) = rr.gen(:, PG);
      Pd(:, t, j, k) = rr.bus(:, PD);
      lamP(:, t, j, k) = rr.bus(:, LAM_P);
      Pf(:, t, j, k) = rr.branch(:, PF);
      muF(:, t, j, k) = rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
    end
  end
end
if ns
    SoC = mdo.Storage.ExpectedStorageState;
else
    SoC = [];
end

ms = struct(...
    'f',    mdo.QP.f + mdo.QP.c1, ...
    'nb',   nb, ...
    'ng',   ng, ...
    'nl',   nl, ...
    'ns',   ns, ...
    'nt',   nt, ...
    'nj_max', nj_max, ...
    'nc_max', nc_max, ...
    'psi',  psi, ...
    'Pg',   Pg, ...
    'Pd',   Pd, ...
    'Rup',  Rup, ...
    'Rdn',  Rdn, ...
    'SoC',  SoC, ...
    'Pf',   Pf, ...
    'u',    u, ...
    'lamP', lamP, ...
    'muF',  muF ...
);

%% print results
if verbose
    fprintf('\n========== OBJECTIVE  ==========\n');
    fprintf('f = %.12g\n', ms.f);

    fprintf('\n========== GEN_STATUS ==========\n');
    fprintf(' Gen ');
    for t = 1:nt
        fprintf('   t =%2d ', t);
    end
    fprintf('\n');
    fprintf('----');
    for t = 1:nt
        fprintf('  -------');
    end
    fprintf('\n');
    for i = 1:ng
        fprintf('%4d', i);
%         fprintf('%9d', u(i, :));
        for t = 1:nt
            qty = u(i, t);
            if abs(qty) > tol
                fprintf('      1  ');
            else
                fprintf('    --0--');
            end
        end
        fprintf('\n');
    end
    fprintf('\n');

    print_most_summary_section('PG', 'Gen', nt, nj_max, nc_max, Pg);
    if mdo.idx.ntramp
        print_most_summary_section('RAMP UP', 'Gen', nt, 1, 0, Rup);
        print_most_summary_section('RAMP DOWN', 'Gen', nt, 1, 0, Rdn);
    end
    print_most_summary_section('FIXED LOAD', 'Bus', nt, nj_max, nc_max, Pd);
    if ns
        print_most_summary_section('ESS E[SoC]', 'ESS', nt, 1, 0, SoC);
    end
    if mdo.DCMODEL
        print_most_summary_section('LAM_P', 'Bus', nt, nj_max, nc_max, lamP);
        print_most_summary_section('PF',   'Brch', nt, nj_max, nc_max, Pf);
        print_most_summary_section('MU_F', 'Brch', nt, nj_max, nc_max, muF);
    end
end

if nargout
    mso = ms;
end

%%---------------------------------------------------------
function print_most_summary_section(label, section_type, nt, nj_max, nc_max, data, tol)
if nargin < 7
    tol = 1e-4;
end
n = size(data, 1);
bl = blanks(fix((12-length(label)) / 2));
fprintf('\n==========%-12s==========\n', sprintf('%s%s', bl, label));
if any(data(:))
    for j = 1:nj_max
        for k = 1:nc_max+1
            if nj_max > 1 || nc_max > 0
                fprintf('\nSCENARIO %d', j);
                if nc_max == 0
                    fprintf('\n');
                elseif k == 1
                    fprintf(', base case\n');
                else
                    fprintf(', contingency %d\n', k-1);
                end
            end
            fprintf('%4s ', section_type);
            for t = 1:nt
                fprintf('   t =%2d ', t);
            end
            fprintf('\n');
            fprintf('----');
            for t = 1:nt
                fprintf('  -------');
            end
            fprintf('\n');
            for i = 1:n
                fprintf('%4d', i);
                for t = 1:nt
                    qty = data(i, t, j, k);
                    if abs(qty) > tol
                        fprintf('%9.2f', qty);
                    else
                        fprintf('      -  ');
                    end
                end
                fprintf('\n');
            end
        end
    end
else
    fprintf('All zeros.\n');
end
fprintf('\n');
