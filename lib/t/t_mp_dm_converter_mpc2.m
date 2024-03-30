function tssk = t_mp_dm_converter_mpc2(quiet)
% t_mp_dm_converter_mpc2 - Tests for mp.dm_converter_mpc2.

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

% define_constants;
if quiet
    verbose = 0;
else
    verbose = 1;
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

casefile = 't_case_ext';
mpc0 = loadcase(casefile);
mpc0.bus(end, MU_VMIN) = 0;         %% add columns
mpc0.gen(end, MU_QMIN) = 0;         %% add columns
mpc0.branch(end, MU_ANGMAX) = 0;    %% add columns
dmc = mp.dm_converter_mpc2().build();

t_begin(27, quiet);

t = 'dmc constructor : ';
dmc = mp.dm_converter_mpc2();
t_ok(isa(dmc, 'mp.dm_converter'), [t 'class']);
t_is(length(dmc.element_classes), 5, 12, [t '# of element_classes']);
t_ok(isempty(dmc.elements), [t 'elements empty']);

t = 'dmc.build() : ';
dmc.build();
t_is(length(dmc.elements), 5, 12, [t '# of elements']);
t_str_match(dmc.elements{1}.name, 'bus', [t 'element{1} is bus']);
t_str_match(dmc.elements{2}.name, 'gen', [t 'element{2} is gen']);
t_str_match(dmc.elements{3}.name, 'load', [t 'element{3} is load']);
t_str_match(dmc.elements{4}.name, 'branch', [t 'element{4} is branch']);
t_str_match(dmc.elements{5}.name, 'shunt', [t 'element{5} is shunt']);

t = 'dm constructor : ';
dm = mp.data_model_cpf();
t_ok(isa(dm, 'mp.data_model'), [t 'class']);
t_is(length(dm.element_classes), 5, 12, [t '# of element_classes']);
t_ok(isempty(dm.elements), [t 'elements empty']);

t = 'dm.build(mpc, dmc) : ';
dm.build(mpc0, dmc);
t_is(length(dm.elements), 5, 12, [t '# of elements']);
t_str_match(dm.elements{1}.name, 'bus', [t 'element{1} is bus']);
t_str_match(dm.elements{2}.name, 'gen', [t 'element{2} is gen']);
t_str_match(dm.elements{3}.name, 'load', [t 'element{3} is load']);
t_str_match(dm.elements{4}.name, 'branch', [t 'element{4} is branch']);
t_str_match(dm.elements{5}.name, 'shunt', [t 'element{5} is shunt']);

t = 'mpc = dmc.export(dm, mpc0) : ';
mpc = dmc.export(dm, mpc0);
t_ok(isequal(mpc, mpc0), [t 'mpc == mpc0']);

t = 'mpc = dmc.export(updated_dm, mpc0) : ';
mpc1 = mpc0;
mpc1.bus([2;3;6;9], VM) = [1.02;1.03;1.06;1.09];
dm.elements.bus.tab.vm([2;3;6;9]) = [1.02;1.03;1.06;1.09];
mpc1.bus([3;5;7], PD) = [30;50;70];
mpc1.bus([3;5;7], QD) = [18;30;42];
load0 = dm.elements.load.tab;
load = [load0; load0(1:2, :)];
k = [4;1;5];
load.uid(k) = [4;1;5];
load.bus(k) = [30;5;6];
load.source_uid(k) = [3;5;7];
load.pd(k) = [30;50;70];
load.qd(k) = [18;30;42];
dm.elements.load.tab = load;
dm.elements.load.rebuild(dm);
mpc = dmc.export(dm, mpc0);
t_ok(~isequal(mpc, mpc0), [t 'mpc ~= mpc0']);
t_ok( isequal(mpc, mpc1), [t 'mpc == updated_mpc']);

t = 'mpc = dmc.export(dm) : ';
mpc1 = struct( ...
    'version', mpc0.version, ...
    'baseMVA', mpc0.baseMVA, ...
    'bus', mpc0.bus, ...
    'gen', mpc0.gen, ...
    'branch', mpc0.branch, ...
    'gencost', mpc0.gencost, ...
    'bus_name', {mpc0.bus_name} );
dmc = mp.dm_converter_mpc2().build();
dm = mp.data_model_opf().build(mpc1, dmc);
mpc = dmc.export(dm);
mpc1.gen(:, MBASE) = 0;         % not included dme_gen
mpc1.branch(7, BR_STATUS) = 0;  % connected to offline bus
mpc1.gen(3, GEN_STATUS) = 0;    % connected to offline bus
t_ok(isequal(mpc, mpc1), [t 'mpc == mpc0']);

t = 'mpc = dmce.export(dme, mpc, ''qcost'') : ';
%% define named indices into data matrices
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
dme = dm.elements.gen;
dmce = dmc.elements.gen;
dme.tab.cost_qg = dme.tab.cost_pg;
dme.tab.cost_qg.poly_coef = dme.tab.cost_qg.poly_coef ./ 5;
ng = size(mpc.gen, 1);
epc = mpc.gencost(1:ng, :);
eqc = epc; eqc([1;3], COST) = eqc([1;3], COST) / 5;
mpc = dmce.export(dme, mpc, {'cost_qg'});
t_is(size(mpc.gencost, 1), 2*ng, 12, [t 'size(gencost, 1)']);
t_is(mpc.gencost(1:ng,:), epc, 12, [t 'pcost']);
t_is(mpc.gencost(ng+1:2*ng,:), eqc, 12, [t 'qcost']);

t = 'mpc = dmce.export(dme, mpc, ''qcost'', ridx) : ';
ridx = [2;3];
dme.tab.cost_pg.poly_coef = dme.tab.cost_pg.poly_coef * 2;
dme.tab.cost_pg.pwl_qty = dme.tab.cost_pg.pwl_qty + 10;
dme.tab.cost_pg.pwl_cost = dme.tab.cost_pg.pwl_cost + 10;
dme.tab.cost_qg.poly_coef = dme.tab.cost_qg.poly_coef * 3;
dme.tab.cost_qg.pwl_qty = dme.tab.cost_qg.pwl_qty + 1;
dme.tab.cost_qg.pwl_cost = dme.tab.cost_qg.pwl_cost + 1;
epc(2, COST:end) = epc(2, COST:end) + 10;
epc(3, COST:end) = epc(3, COST:end) * 2;
eqc(2, COST:end) = eqc(2, COST:end) + 1;
eqc(3, COST:end) = eqc(3, COST:end) * 3;
mpc = dmce.export(dme, mpc, {'cost_pg', 'cost_qg'}, ridx);
t_is(mpc.gencost(1:ng,:), epc, 12, [t 'pcost']);
t_is(mpc.gencost(ng+1:2*ng,:), eqc, 12, [t 'qcost']);

t_end;
