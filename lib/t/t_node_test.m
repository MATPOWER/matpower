function [pf_, opf_] = t_node_test(quiet)
% t_node_test - Tests for network model with multipe node-creating elements.

%   MATPOWER
%   Copyright (c) 2021-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 1
    quiet = 0;
end

cases = {'t_case9_gizmo', 'case4gs', 'case4_dist', 'case9', 'case14', 'case57', 'case300'};

t_begin(12*length(cases), quiet);

define_constants;
if quiet
    verbose = 0;
else
    verbose = 1;
end

mpopt = mpoption('out.all', 0, 'verbose', 0, 'pf.tol', 1e-10);
mpopt = mpoption(mpopt, 'opf.ignore_angle_lim', 1);
mpopt0 = mpopt;
mpopt0.exp.use_legacy_core = 1;
mpopt.exp.mpx = mp.xt_node_test();

for k = 1:length(cases)
    t = sprintf('PF - %s - ', cases{k});
    mpc = loadcase(cases{k});
    if strcmp(cases{k}, 't_case9_gizmo')
        mpc.bus(2, BS) = 1;
    end
    if ~isfield(mpc, 'gencost')
        ng = size(mpc.gen, 1);
        mpc.gencost = ones(ng, 1) * [2 0 0 2 1 0];
        if strcmp(cases{k}, 'case4gs')
            mpc.gen(:, PMAX) = 320;
            mpc.gen(:, QMAX) = 200;
        end
    end

    r = runpf(mpc, mpopt0);
    evm = r.bus(:, VM);
    eva = r.bus(:, VA);
    epg = r.gen(:, PG);
    eqg = r.gen(:, QG);
    t_ok(r.success, [t 'success 1']);

    pf = run_pf(mpc, mpopt);
    r2 = pf.dmc.export(pf.dm, pf.dm.source);
    va = r2.bus(:, VA);
    vm = r2.bus(:, VM);
    pg = r2.gen(:, PG);
    qg = r2.gen(:, QG);
    t_ok(pf.success, [t 'success 2']);
    t_is(va, eva, 9, [t 'va']);
    t_is(vm, evm, 9, [t 'vm']);
    t_is(pg, epg, 9, [t 'pg']);
    t_is(qg, eqg, 9, [t 'qg']);

    t = sprintf('OPF - %s - ', cases{k});
    r = runopf(mpc, mpopt0);
    t_ok(r.success, [t 'success 1']);
    evm = r.bus(:, VM);
    eva = r.bus(:, VA);
    epg = r.gen(:, PG);
    eqg = r.gen(:, QG);

    opf = run_opf(mpc, mpopt);
    r2 = opf.dmc.export(opf.dm, opf.dm.source);
    va = r2.bus(:, VA);
    vm = r2.bus(:, VM);
    pg = r2.gen(:, PG);
    qg = r2.gen(:, QG);
    t_ok(opf.success, [t 'success 2']);
    t_is(va, eva, 9, [t 'va']);
    t_is(vm, evm, 9, [t 'vm']);
    t_is(pg, epg, 9, [t 'pg']);
    t_is(qg, eqg, 9, [t 'qg']);
end

t_end;

if nargout
    pf_ = pf;
    if nargout > 1
        opf_ = opf;
    end
end
