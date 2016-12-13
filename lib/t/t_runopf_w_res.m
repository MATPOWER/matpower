function t_runopf_w_res(quiet)
%T_RUNOPF_W_RES  Tests RUNOPF_W_RES and the associated callbacks.

%   MATPOWER
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

t_begin(46, quiet);

if quiet
    verbose = 0;
else
    verbose = 0;
end

casefile = 't_case30_userfcns';
mpopt = mpoption('opf.violation', 1e-6, 'mips.gradtol', 1e-8, ...
        'mips.comptol', 1e-8, 'mips.costtol', 1e-9);
mpopt = mpoption(mpopt, 'out.all', 0, 'verbose', verbose, 'opf.ac.solver', 'MIPS');

%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

t = 'runopf_w_res(''t_case30_userfcns'') : ';
r = runopf_w_res(casefile, mpopt);
t_is(r.reserves.R, [25; 15; 0; 0; 19.3906; 0.6094], 4, [t 'R']);
t_is(r.reserves.prc, [2; 2; 2; 2; 5.5; 5.5], 6, [t 'prc']);
t_is(r.reserves.mu.l, [0; 0; 1; 2; 0; 0], 7, [t 'mu.l']);
t_is(r.reserves.mu.u, [0.1; 0; 0; 0; 0; 0], 7, [t 'mu.u']);
t_is(r.reserves.mu.Pmax, [0; 0; 0; 0; 0.5; 0], 7, [t 'mu.Pmax']);
mpc = loadcase(casefile);
t_is(r.reserves.cost, mpc.reserves.cost, 12, [t 'cost']);
t_is(r.reserves.qty, mpc.reserves.qty, 12, [t 'qty']);
t_is(r.reserves.totalcost, 177.8047, 4, [t 'totalcost']);

t = 'gen 5 no reserves : ';
mpc = loadcase(casefile);
mpc.reserves.zones(:, 5) = 0;
mpc.reserves.cost(5) = [];
mpc.reserves.qty(5) = [];
r = runopf_w_res(mpc, mpopt);
t_is(r.reserves.R, [25; 15; 0; 0; 0; 20], 4, [t 'R']);
t_is(r.reserves.prc, [2; 2; 2; 2; 0; 5.5], 6, [t 'prc']);
t_is(r.reserves.mu.l, [0; 0; 1; 2; 0; 0], 7, [t 'mu.l']);
t_is(r.reserves.mu.u, [0.1; 0; 0; 0; 0; 0], 6, [t 'mu.u']);
t_is(r.reserves.mu.Pmax, [0; 0; 0; 0; 0; 0], 7, [t 'mu.Pmax']);
t_is(r.reserves.cost, mpc.reserves.cost, 12, [t 'cost']);
t_is(r.reserves.qty, mpc.reserves.qty, 12, [t 'qty']);
t_is(r.reserves.totalcost, 187.5, 4, [t 'totalcost']);

t = 'extra offline gen : ';
mpc = loadcase(casefile);
idx = [1:3 5 4:6];
mpc.gen = mpc.gen(idx, :);
mpc.gencost = mpc.gencost(idx, :);
mpc.reserves.zones = mpc.reserves.zones(:, idx);
mpc.reserves.cost = mpc.reserves.cost(idx);
mpc.reserves.qty = mpc.reserves.qty(idx);
mpc.gen(4, GEN_STATUS) = 0;
r = runopf_w_res(mpc, mpopt);
t_is(r.reserves.R, [25; 15; 0; 0; 0; 19.3906; 0.6094], 4, [t 'R']);
t_is(r.reserves.prc, [2; 2; 2; 5.5; 2; 5.5; 5.5], 6, [t 'prc']);
t_is(r.reserves.mu.l, [0; 0; 1; 0; 2; 0; 0], 7, [t 'mu.l']);
t_is(r.reserves.mu.u, [0.1; 0; 0; 0; 0; 0; 0], 7, [t 'mu.u']);
t_is(r.reserves.mu.Pmax, [0; 0; 0; 0; 0; 0.5; 0], 7, [t 'mu.Pmax']);
t_is(r.reserves.cost, mpc.reserves.cost, 12, [t 'cost']);
t_is(r.reserves.qty, mpc.reserves.qty, 12, [t 'qty']);
t_is(r.reserves.totalcost, 177.8047, 4, [t 'totalcost']);

t = 'both extra & gen 6 no res : ';
mpc = loadcase(casefile);
idx = [1:3 5 4:6];
mpc.gen = mpc.gen(idx, :);
mpc.gencost = mpc.gencost(idx, :);
mpc.reserves.zones = mpc.reserves.zones(:, idx);
mpc.reserves.cost = mpc.reserves.cost(idx);
mpc.reserves.qty = mpc.reserves.qty(idx);
mpc.gen(4, GEN_STATUS) = 0;
mpc.reserves.zones(:, 6) = 0;
mpc.reserves.cost(6) = [];
mpc.reserves.qty(6) = [];
r = runopf_w_res(mpc, mpopt);
t_is(r.reserves.R, [25; 15; 0; 0; 0; 0; 20], 4, [t 'R']);
t_is(r.reserves.prc, [2; 2; 2; 5.5; 2; 0; 5.5], 6, [t 'prc']);
t_is(r.reserves.mu.l, [0; 0; 1; 0; 2; 0; 0], 7, [t 'mu.l']);
t_is(r.reserves.mu.u, [0.1; 0; 0; 0; 0; 0; 0], 6, [t 'mu.u']);
t_is(r.reserves.mu.Pmax, [0; 0; 0; 0; 0; 0; 0], 7, [t 'mu.Pmax']);
t_is(r.reserves.cost, mpc.reserves.cost, 12, [t 'cost']);
t_is(r.reserves.qty, mpc.reserves.qty, 12, [t 'qty']);
t_is(r.reserves.totalcost, 187.5, 4, [t 'totalcost']);

t = 'no qty (Rmax) : ';
mpc = loadcase(casefile);
mpc.reserves = rmfield(mpc.reserves, 'qty');
r = runopf_w_res(mpc, mpopt);
t_is(r.reserves.R, [39.3826; 0.6174; 0; 0; 19.3818; 0.6182], 4, [t 'R']);
t_is(r.reserves.prc, [2; 2; 2; 2; 5.5; 5.5], 5, [t 'prc']);
t_is(r.reserves.mu.l, [0; 0; 1; 2; 0; 0], 5, [t 'mu.l']);
t_is(r.reserves.mu.u, [0; 0; 0; 0; 0; 0], 7, [t 'mu.u']);
t_is(r.reserves.mu.Pmax, [0.1; 0; 0; 0; 0.5; 0], 5, [t 'mu.Pmax']);
t_is(r.reserves.cost, mpc.reserves.cost, 12, [t 'cost']);
t_is(r.reserves.totalcost, 176.3708, 4, [t 'totalcost']);

t = 'RAMP_10, no qty (Rmax) : ';
mpc = loadcase(casefile);
mpc.reserves = rmfield(mpc.reserves, 'qty');
mpc.gen(1, RAMP_10) = 25;
r = runopf_w_res(mpc, mpopt);
t_is(r.reserves.R, [25; 15; 0; 0; 19.3906; 0.6094], 4, [t 'R']);
t_is(r.reserves.prc, [2; 2; 2; 2; 5.5; 5.5], 6, [t 'prc']);
t_is(r.reserves.mu.l, [0; 0; 1; 2; 0; 0], 7, [t 'mu.l']);
t_is(r.reserves.mu.u, [0.1; 0; 0; 0; 0; 0], 7, [t 'mu.u']);
t_is(r.reserves.mu.Pmax, [0; 0; 0; 0; 0.5; 0], 7, [t 'mu.Pmax']);
t_is(r.reserves.cost, mpc.reserves.cost, 12, [t 'cost']);
t_is(r.reserves.totalcost, 177.8047, 4, [t 'totalcost']);

t_end;
