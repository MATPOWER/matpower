function t_runopf_w_res(quiet)
%T_RUNOPF_W_RES  Tests runopf_w_res and the associated callbacks.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2009 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

t_begin(34, quiet);

if quiet
    verbose = 0;
else
    verbose = 0;
end

casefile = 't_case30_userfcns';
mpopt = mpoption('OPF_VIOLATION', 1e-6, 'PDIPM_GRADTOL', 1e-8, ...
        'PDIPM_COMPTOL', 1e-8, 'PDIPM_COSTTOL', 1e-9);
mpopt = mpoption(mpopt, 'OUT_ALL', 0, 'VERBOSE', verbose, 'OPF_ALG', 560);

%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

t = 'runopf_w_res(''case30_reserves'') : ';
r = runopf_w_res(casefile, mpopt);
t_is(r.reserves.R, [25; 15; 0; 0; 13.8729; 6.1272], 4, [t 'R']);
t_is(r.reserves.prc, [2; 2; 2; 2; 4; 4], 6, [t 'prc']);
t_is(r.reserves.mu.l, [0; 0; 1; 2; 0; 0], 7, [t 'mu.l']);
t_is(r.reserves.mu.u, [1; 0; 0; 0; 0; 0], 7, [t 'mu.u']);
t_is(r.reserves.mu.Pmax, [0; 0; 0; 0; 1; 0], 7, [t 'mu.Pmax']);
mpc = loadcase(casefile);
t_is(r.reserves.cost, mpc.reserves.cost, 12, [t 'cost']);
t_is(r.reserves.qty, mpc.reserves.qty, 12, [t 'qty']);

t = 'gen 5 no reserves : ';
mpc = loadcase(casefile);
mpc.reserves.zones(:, 5) = 0;
mpc.reserves.cost(5) = [];
mpc.reserves.qty(5) = [];
r = runopf_w_res(mpc, mpopt);
t_is(r.reserves.R, [25; 15; 0; 0; 0; 20], 4, [t 'R']);
t_is(r.reserves.prc, [2; 2; 2; 2; 0; 4], 6, [t 'prc']);
t_is(r.reserves.mu.l, [0; 0; 1; 2; 0; 0], 7, [t 'mu.l']);
t_is(r.reserves.mu.u, [1; 0; 0; 0; 0; 0], 6, [t 'mu.u']);
t_is(r.reserves.mu.Pmax, [0; 0; 0; 0; 0; 0], 7, [t 'mu.Pmax']);
t_is(r.reserves.cost, mpc.reserves.cost, 12, [t 'cost']);
t_is(r.reserves.qty, mpc.reserves.qty, 12, [t 'qty']);

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
t_is(r.reserves.R, [25; 15; 0; 0; 0; 13.8729; 6.1272], 4, [t 'R']);
t_is(r.reserves.prc, [2; 2; 2; 4; 2; 4; 4], 6, [t 'prc']);
t_is(r.reserves.mu.l, [0; 0; 1; 0; 2; 0; 0], 7, [t 'mu.l']);
t_is(r.reserves.mu.u, [1; 0; 0; 0; 0; 0; 0], 7, [t 'mu.u']);
t_is(r.reserves.mu.Pmax, [0; 0; 0; 0; 0; 1; 0], 7, [t 'mu.Pmax']);
t_is(r.reserves.cost, mpc.reserves.cost, 12, [t 'cost']);
t_is(r.reserves.qty, mpc.reserves.qty, 12, [t 'qty']);

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
t_is(r.reserves.prc, [2; 2; 2; 4; 2; 0; 4], 6, [t 'prc']);
t_is(r.reserves.mu.l, [0; 0; 1; 0; 2; 0; 0], 7, [t 'mu.l']);
t_is(r.reserves.mu.u, [1; 0; 0; 0; 0; 0; 0], 6, [t 'mu.u']);
t_is(r.reserves.mu.Pmax, [0; 0; 0; 0; 0; 0; 0], 7, [t 'mu.Pmax']);
t_is(r.reserves.cost, mpc.reserves.cost, 12, [t 'cost']);
t_is(r.reserves.qty, mpc.reserves.qty, 12, [t 'qty']);

t = 'no qty (Rmax) : ';
mpc = loadcase(casefile);
mpc.reserves = rmfield(mpc.reserves, 'qty');
r = runopf_w_res(mpc, mpopt);
t_is(r.reserves.R, [38.6199; 1.3801; 0; 0; 13.8471; 6.1529], 4, [t 'R']);
t_is(r.reserves.prc, [2; 2; 2; 2; 4; 4], 5, [t 'prc']);
t_is(r.reserves.mu.l, [0; 0; 1; 2; 0; 0], 5, [t 'mu.l']);
t_is(r.reserves.mu.u, [0; 0; 0; 0; 0; 0], 7, [t 'mu.u']);
t_is(r.reserves.mu.Pmax, [1; 0; 0; 0; 1; 0], 6, [t 'mu.Pmax']);
t_is(r.reserves.cost, mpc.reserves.cost, 12, [t 'cost']);

t_end;
