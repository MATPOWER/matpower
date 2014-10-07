function t_opf_softlims(quiet)
%T_OPF_SOFTLIMS  Tests for userfcn callbacks (softlims) w/OPF.
%   Includes high-level tests of soft limits implementations.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2009-2014 by Power System Engineering Research Center (PSERC)
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

if nargin < 1
    quiet = 0;
end

casefile = 'case9';
if quiet
    verbose = 0;
else
    verbose = 0;
end

t_begin(55, quiet);

%% define constants
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%% set options
mpopt = mpoption('model', 'DC', 'opf.dc.solver', 'MIPS');
mpopt = mpoption(mpopt, 'opf.violation', 1e-6, 'mips.gradtol', 1e-8, ...
        'mips.comptol', 1e-8, 'mips.costtol', 1e-9);
if verbose <= 1
    mpopt = mpoption(mpopt, 'out.all', 0);
end
mpopt = mpoption(mpopt, 'verbose', verbose);

%% load and modify case file
mpc = loadcase(casefile);
nl = size(mpc.branch, 1);
mpc.gencost(:, NCOST) = 2;
mpc.gencost(:, COST) = [];
mpc.gencost(:, COST) = [50; 40; 25];
mpc.gencost(:, COST+1) = 0;
mpc.branch(2:6, RATE_A) = [0; 120; 0; 120; 0];
mpc.branch(3, RATE_A) = 120;
mpc.branch(5, RATE_A) = 120;
mpc.branch = [mpc.branch(1:4, :); mpc.branch(4, :);  mpc.branch(5:end, :)];
mpc.branch(5, BR_STATUS) = 0;

nl = size(mpc.branch, 1);   %% number of branches

%% create soft limit inputs
mpc.softlims.idx = (2:nl)';
mpc.softlims.cost = 100 * ones(nl-1, 1);

t = 'hard limits : ';
t_ok(~toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 0');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(~isfield(r.softlims, 'overload'), [t 'no softlims.overload']);
t_ok(~isfield(r.softlims, 'overload'), [t 'no softlims.ovl_cost']);
t_is(r.f, 9126.87, 4, [t 'f']);
t_is(r.gen(:, PG), [12.687; 62.3130; 240], 4, [t 'Pg']);
t_is(r.branch(:, PF), [12.687; -30; -120; 240; 0; 120; 20; -62.3130; 82.3130; -42.687], 4, [t 'Pf']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 35.6504; 0; 0; 7.9756; 0; 0; 0; 0], 4, [t 'mu Pf']);

t = 'soft limits (satisfied) : ';
mpc = toggle_softlims(mpc, 'on');
t_ok(toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 1');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims, 'overload'), [t 'softlims.overload exists']);
t_ok(isfield(r.softlims, 'overload'), [t 'softlims.ovl_cost exists']);
t_is(r.f, 9126.87, 4, [t 'f']);
t_is(r.gen(:, PG), [12.687; 62.3130; 240], 4, [t 'Pg']);
t_is(r.branch(:, PF), [12.687; -30; -120; 240; 0; 120; 20; -62.3130; 82.3130; -42.687], 4, [t 'Pf']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 35.6504; 0; 0; 7.9756; 0; 0; 0; 0], 4, [t 'mu Pf']);
t_is(r.softlims.overload, zeros(10, 1), 12, [t 'softlims.overload']);
t_is(r.softlims.ovl_cost, zeros(10, 1), 12, [t 'softlims.ovl_cost']);
t_is(r.order.branch.status.on(r.order.int.softlims.idx), [3; 6; 8; 9; 10], 12, [t 'mu Pf']);

t = 'savecase(fname, mpc) : ';
fn = sprintf('softlims_savecase_test_%d', fix(1e8*rand));
savecase(fn, mpc);
mpc1 = loadcase(fn);
delete([fn '.m']);
t_ok(isfield(mpc1, 'softlims'), [t 'mpc.softlims']);
t_ok(isfield(mpc1.softlims, 'idx'), [t 'mpc.softlims.idx']);
t_is(mpc1.softlims.idx, mpc.softlims.idx, 5, [t 'mpc.softlims.idx']);
t_ok(isfield(mpc1.softlims, 'cost'), [t 'mpc.softlims.cost']);
t_is(mpc1.softlims.cost, mpc.softlims.cost, 5, [t 'mpc.softlims.cost']);

t = 'savecase(fname, results) : ';
fn = sprintf('softlims_savecase_test_%d', fix(1e8*rand));
savecase(fn, r);
mpc1 = loadcase(fn);
delete([fn '.m']);
t_ok(isfield(mpc1, 'softlims'), [t 'results.softlims']);
t_ok(isfield(mpc1.softlims, 'idx'), [t 'results.softlims.idx']);
t_is(mpc1.softlims.idx, r.softlims.idx, 5, [t 'results.softlims.idx']);
t_ok(isfield(mpc1.softlims, 'cost'), [t 'results.softlims.cost']);
t_is(mpc1.softlims.cost, r.softlims.cost, 5, [t 'results.softlims.cost']);
t_ok(isfield(mpc1.softlims, 'overload'), [t 'results.softlims.overload']);
t_is(mpc1.softlims.overload, r.softlims.overload, 5, [t 'results.softlims.overload']);
t_ok(isfield(mpc1.softlims, 'ovl_cost'), [t 'results.softlims.ovl_cost']);
t_is(mpc1.softlims.ovl_cost, r.softlims.ovl_cost, 4, [t 'results.softlims.ovl_cost']);

t = 'soft limits (violated) : ';
mpc = rmfield(mpc, 'softlims');
mpc.softlims.cost = 20;
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims, 'overload'), [t 'softlims.overload exists']);
t_ok(isfield(r.softlims, 'overload'), [t 'softlims.ovl_cost exists']);
t_is(r.f, 9106.5059, 4, [t 'f']);
t_is(r.gen(:, PG), [10; 63.6988; 241.3012], 4, [t 'Pg']);
t_is(r.branch(:, PF), [10; -31.3012; -121.3012; 241.3012; 0; 120; 20; -63.6988; 83.6988; -41.3012], 4, [t 'Pf']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 20; 0; 0; 13.2992; 0; 0; 0; 0], 4, [t 'mu Pf']);
t_is(r.softlims.overload, [0; 0; 1.3011811; 0; 0; 0; 0; 0; 0; 0], 6, [t 'softlims.overload']);
t_is(r.softlims.ovl_cost, [0; 0; 26.023622; 0; 0; 0; 0; 0; 0; 0], 6, [t 'softlims.ovl_cost']);

t = 'savecase(fname, mpc) : ';
fn = sprintf('softlims_savecase_test_%d', fix(1e8*rand));
savecase(fn, mpc);
mpc1 = loadcase(fn);
delete([fn '.m']);
t_ok(isfield(mpc1, 'softlims'), [t 'mpc.softlims']);
t_ok(~isfield(mpc1.softlims, 'idx'), [t 'mpc.softlims.idx']);
t_ok(isfield(mpc1.softlims, 'cost'), [t 'mpc.softlims.cost']);
t_is(mpc1.softlims.cost, mpc.softlims.cost, 5, [t 'mpc.softlims.cost']);

t = 'savecase(fname, results) : ';
fn = sprintf('softlims_savecase_test_%d', fix(1e8*rand));
savecase(fn, r);
mpc1 = loadcase(fn);
delete([fn '.m']);
t_ok(isfield(mpc1, 'softlims'), [t 'results.softlims']);
t_ok(isfield(mpc1.softlims, 'idx'), [t 'results.softlims.idx']);
t_ok(isempty(mpc1.softlims.idx), [t 'results.softlims.idx']);
t_ok(isfield(mpc1.softlims, 'cost'), [t 'results.softlims.cost']);
t_is(mpc1.softlims.cost, r.softlims.cost, 5, [t 'results.softlims.cost']);
t_ok(isfield(mpc1.softlims, 'overload'), [t 'results.softlims.overload']);
t_is(mpc1.softlims.overload, r.softlims.overload, 5, [t 'results.softlims.overload']);
t_ok(isfield(mpc1.softlims, 'ovl_cost'), [t 'results.softlims.ovl_cost']);
t_is(mpc1.softlims.ovl_cost, r.softlims.ovl_cost, 4, [t 'results.softlims.ovl_cost']);

t_end;
