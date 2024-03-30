function obj = t_nm_element(quiet, out_ac)
% t_nm_element - Tests for mp.nm_element.

%   MATPOWER
%   Copyright (c) 2019-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 1
    quiet = 0;
end

t_begin(453, quiet);

define_constants;
if quiet
    verbose = 0;
else
    verbose = 1;
end

casefile = 't_case9_gizmo';
mpopt = mpoption('out.all', 0, 'verbose', 0);
dmc = mp.dm_converter_mpc2().build();

%%-----  DC formulation  -----
t = 'mp.net_model_dc() : ';
dc = mp.net_model_dc();
t_str_match(dc.name, 'network', [t 'name']);
t_str_match(class(dc), 'mp.net_model_dc', [t 'class']);
t_str_match(dc.find_form_class(), 'mp.form_dc', [t 'formulation class']);
t_str_match(dc.form_name, 'DC', [t 'formulation name']);
t_str_match(dc.form_tag, 'dc', [t 'formulation tag']);
p = dc.model_params();
t_is(length(p), 3, 12, [t '# of model params']);
t_str_match(p{1}, 'B', [t 'B']);
t_str_match(p{2}, 'K', [t 'K']);
t_str_match(p{3}, 'p', [t 'p']);
t_is(dc.nk, 0, 12, [t 'nk']);
t_is(dc.np, 0, 12, [t 'np']);
t_is(dc.nz, 0, 12, [t 'nz']);
t_is(dc.nv, 0, 12, [t 'nv']);
t_is(length(fieldnames(dc.set_types)), 5, 12, [t '# of set types']);
t_str_match(dc.set_types.node, 'NODES', [t 'set_types.node']);
t_str_match(dc.set_types.state, 'STATES', [t 'set_types.state']);
t_str_match(dc.set_types.va, 'VOLTAGE VARS (va)', [t 'set_types.va']);
t_str_match(dc.set_types.z, 'NON-VOLTAGE VARS (z)', [t 'set_types.z']);
t_is(length(dc.elements), 0, 12, [t '# of element types']);

t = 'dc.build(dm, dmc) : ';
dm = mp.data_model().build(rundcpf(loadcase(casefile), mpopt), dmc);
mpc = dm.source;
t_ok(mpc.success, [t 'solved power flow']);
dc.build(dm);
t_is(dc.nk, 1, 12, [t 'nk']);
t_is(dc.np, 24, 12, [t 'np']);
t_is(dc.nz, 3, 12, [t 'nz']);
t_is(dc.nv, 9, 12, [t 'nv']);
t_is(length(dc.elements), 4, 12, [t '# of element types']);

nme = dc.elements;
t = 'mp.nme_bus_dc : '; k = 1;
t_str_match(nme{k}.name, 'bus', [t 'name']);
t_str_match(class(nme{k}), 'mp.nme_bus_dc', [t 'class']);
t_str_match(nme{k}.find_form_class(), 'mp.form_dc', [t 'formulation class']);
t_str_match(nme{k}.form_name, 'DC', [t 'formulation name']);
t_str_match(nme{k}.form_tag, 'dc', [t 'formulation tag']);
t_is(nme{k}.nk, 9, 12, [t 'nk']);
t_is(nme{k}.np, 0, 12, [t 'np']);
t_is(nme{k}.nz, 0, 12, [t 'nz']);
t_ok(isempty(nme{k}.B), [t 'B']);
t_ok(isempty(nme{k}.K), [t 'K']);
t_ok(isempty(nme{k}.p), [t 'p']);
t_ok(isempty(nme{k}.C), [t 'C']);
t_ok(isempty(nme{k}.D), [t 'D']);
% nme{k}
% C = nme{k}.C
% D = nme{k}.D

t = 'mp.nme_gen_dc : '; k = 2;
t_str_match(nme{k}.name, 'gen', [t 'name']);
t_str_match(class(nme{k}), 'mp.nme_gen_dc', [t 'class']);
t_str_match(nme{k}.find_form_class(), 'mp.form_dc', [t 'formulation class']);
t_str_match(nme{k}.form_name, 'DC', [t 'formulation name']);
t_str_match(nme{k}.form_tag, 'dc', [t 'formulation tag']);
t_is(nme{k}.nk, 3, 12, [t 'nk']);
t_is(nme{k}.np, 1, 12, [t 'np']);
t_is(nme{k}.nz, 1, 12, [t 'nz']);
t_ok(isempty(nme{k}.B), [t 'B']);
t_is(nme{k}.K, -speye(3), 12, [t 'K']);
t_ok(isempty(nme{k}.p), [t 'p']);
t_is(nme{k}.C, sparse([1;3;2], [1;2;3], 1, 9, 3), 12, [t 'C']);
t_is(nme{k}.D, speye(3), 12, [t 'D']);
% nme{k}
% C = nme{k}.C
% D = nme{k}.D

t = 'mp.nme_load_dc : '; k = 3;
t_str_match(nme{k}.name, 'load', [t 'name']);
t_str_match(class(nme{k}), 'mp.nme_load_dc', [t 'class']);
t_str_match(nme{k}.find_form_class(), 'mp.form_dc', [t 'formulation class']);
t_str_match(nme{k}.form_name, 'DC', [t 'formulation name']);
t_str_match(nme{k}.form_tag, 'dc', [t 'formulation tag']);
t_is(nme{k}.nk, 3, 12, [t 'nk']);
t_is(nme{k}.np, 1, 12, [t 'np']);
t_is(nme{k}.nz, 0, 12, [t 'nz']);
t_ok(isempty(nme{k}.B), [t 'B']);
t_ok(isempty(nme{k}.K), [t 'K']);
t_is(nme{k}.p, [0.9 1 1.25]', 12, [t 'p']);
eC = sparse([5 7 9], 1:3, 1, 9, 3);
t_is(nme{k}.C, eC, 12, [t 'C']);
t_ok(isempty(nme{k}.D), [t 'D']);
% nme{k}
% C = nme{k}.C
% D = nme{k}.D

t = 'mp.nme_branch_dc : '; k = 4;
t_str_match(nme{k}.name, 'branch', [t 'name']);
t_str_match(class(nme{k}), 'mp.nme_branch_dc', [t 'class']);
t_str_match(nme{k}.find_form_class(), 'mp.form_dc', [t 'formulation class']);
t_str_match(nme{k}.form_name, 'DC', [t 'formulation name']);
t_str_match(nme{k}.form_tag, 'dc', [t 'formulation tag']);
t_is(nme{k}.nk, 9, 12, [t 'nk']);
t_is(nme{k}.np, 2, 12, [t 'np']);
t_is(nme{k}.nz, 0, 12, [t 'nz']);
nl = size(mpc.branch, 1);
stat = mpc.branch(:, BR_STATUS);    %% ones at in-service branches
tm = ones(nl, 1);                   %% default tap ratio = 1
i = find(mpc.branch(:, TAP));       %% indices of non-zero tap ratios
tm(i) = mpc.branch(i, TAP);         %% assign non-zero tap ratios
b = stat ./ mpc.branch(:, BR_X);    %% series susceptance
b = b ./ tm;
Pfinj = b .* (-mpc.branch(:, SHIFT) * pi/180);
Bdc = sparse( ...
    [1:nl 1:nl nl+1:2*nl nl+1:2*nl]', ...
    [1:nl nl+1:2*nl 1:nl nl+1:2*nl]', ...
    [b; -b; -b; b], ...
    2*nl, 2*nl );
pdc = [Pfinj; -Pfinj];
t_is(nme{k}.B, Bdc, 12, [t 'B']);
t_ok(isempty(nme{k}.K), [t 'K']);
t_is(nme{k}.p, pdc, 12, [t 'p']);
t_is(nme{k}.C, sparse([1;4;5;3;6;7;8;8;9;4;5;6;6;7;8;2;9;4], [1:18]', 1, 9, 18), 12, [t 'C']);
t_ok(isempty(nme{k}.D), [t 'D']);
% nme{k}
% C = nme{k}.C
% D = nme{k}.D

t = 'dc.C';
C = dc.C;
t_is(C, sparse([1;3;2;5;7;9;1;4;5;3;6;7;8;8;9;4;5;6;6;7;8;2;9;4], [1:24]', 1, 9, 24), 12, t);

t = 'dc.D';
D = dc.D;
t_is(D, speye(3), 12, t);

t = 'dc.port_inj_power(x) : ';
v = mpc.bus(:, VA) * pi/180;
z = mpc.gen(:, PG) / mpc.baseMVA;
x = [v;z];
P  = dc.port_inj_power(x, 1);
eP = [-0.67 -0.85 -1.63 0.9 1.0 1.25 0.67 0.289673913 -0.610326086 0.85 0.239673913 -0.760326086 -1.63 0.869673913 -0.380326086 -0.67 -0.289673913 0.610326086 -0.85 -0.239673913 0.760326086 1.63 -0.869673913 0.380326086]';
t_is(P, eP, 8, t);
t_is(C * P, 0, 12, [t ' : C * P == 0']);

t = 'dc.port_inj_power(A''*x, 0)';
nv = dc.nv;
nz = dc.nz;
np = dc.np;
A = [   C sparse(nv, nz);
        sparse(nz, np) D    ];
P1 = dc.port_inj_power(A'*x, 0);
t_is(P1, eP, 8, t);

t = 'dc.port_inj_power(x, 1, [3;1])';
P2  = dc.port_inj_power(x, 1, [3;1]);
t_is(P2, eP([3;1]), 8, t);

t = 'dc.elements.gen : ';
gen = dc.elements.gen;
t_str_match(gen.name, 'gen', [t 'name']);
t_str_match(class(gen), 'mp.nme_gen_dc', [t 'class']);

t = 'gen.port_inj_power(x, 1)';
pg = gen.port_inj_power(x, 1);
epg = -[0.67; 0.85; 1.63];
t_is(pg, epg, 12, t);

t = 'gen.port_inj_power(x, 1, [3;1])';
pg31 = gen.port_inj_power(x, 1, [3;1]);
t_is(pg31, epg([3;1]), 12, t);

t = 'gen.port_inj_power(x, 1, 2)';
pg2 = gen.port_inj_power(x, 1, 2);
t_is(pg2, epg(2), 12, t);

%%-----  AC formulation  -----
t = 'mp.net_model_acp() : ';
ac = mp.net_model_acp();
t_str_match(ac.name, 'network', [t 'name']);
t_str_match(class(ac), 'mp.net_model_acp', [t 'class']);
t_str_match(ac.find_form_class(), 'mp.form_acp', [t 'formulation class']);
t_str_match(ac.form_name, 'AC-polar', [t 'formulation name']);
t_str_match(ac.form_tag, 'acp', [t 'formulation tag']);
p = ac.model_params();
t_is(length(p), 6, 12, [t '# of model params']);
t_str_match(p{1}, 'Y', [t 'Y']);
t_str_match(p{2}, 'L', [t 'L']);
t_str_match(p{3}, 'M', [t 'M']);
t_str_match(p{4}, 'N', [t 'N']);
t_str_match(p{5}, 'i', [t 'i']);
t_str_match(p{6}, 's', [t 's']);
t_is(ac.nk, 0, 12, [t 'nk']);
t_is(ac.np, 0, 12, [t 'np']);
t_is(ac.nz, 0, 12, [t 'nz']);
t_is(ac.nv, 0, 12, [t 'nv']);
t_is(length(fieldnames(ac.set_types)), 7, 12, [t '# of set types']);
t_str_match(ac.set_types.node, 'NODES', [t 'set_types.node']);
t_str_match(ac.set_types.state, 'STATES', [t 'set_types.state']);
t_str_match(ac.set_types.va, 'VOLTAGE ANG VARS (va)', [t 'set_types.va']);
t_str_match(ac.set_types.vm, 'VOLTAGE MAG VARS (vm)', [t 'set_types.vm']);
t_str_match(ac.set_types.zr, 'NON-VOLTAGE VARS REAL (zr)', [t 'set_types.zr']);
t_str_match(ac.set_types.zi, 'NON-VOLTAGE VARS IMAG (zi)', [t 'set_types.zi']);
t_is(length(ac.elements), 0, 12, [t '# of element types']);

t = 'ac.build(dm) : ';
dm = mp.data_model().build(runpf(loadcase(casefile), mpopt), dmc);
mpc = dm.source;
t_ok(mpc.success, [t 'solved power flow']);
ac.build(dm);
t_is(ac.nk, 1, 12, [t 'nk']);
t_is(ac.np, 24, 12, [t 'np']);
t_is(ac.nz, 3, 12, [t 'nz']);
t_is(ac.nv, 18, 12, [t 'nv']);
t_is(length(ac.elements), 4, 12, [t '# of element types']);

nme = ac.elements;
t = 'mp.nme_bus_acp : '; k = 1;
t_str_match(nme{k}.name, 'bus', [t 'name']);
t_str_match(class(nme{k}), 'mp.nme_bus_acp', [t 'class']);
t_str_match(nme{k}.find_form_class(), 'mp.form_acp', [t 'formulation class']);
t_str_match(nme{k}.form_name, 'AC-polar', [t 'formulation name']);
t_str_match(nme{k}.form_tag, 'acp', [t 'formulation tag']);
t_is(nme{k}.nk, 9, 12, [t 'nk']);
t_is(nme{k}.np, 0, 12, [t 'np']);
t_is(nme{k}.nz, 0, 12, [t 'nz']);
t_ok(isempty(nme{k}.Y), [t 'Y']);
t_ok(isempty(nme{k}.L), [t 'L']);
t_ok(isempty(nme{k}.M), [t 'M']);
t_ok(isempty(nme{k}.N), [t 'N']);
t_ok(isempty(nme{k}.i), [t 'i']);
t_ok(isempty(nme{k}.s), [t 's']);
t_ok(isempty(nme{k}.C), [t 'C']);
t_ok(isempty(nme{k}.D), [t 'D']);
% nme{k}
% C = nme{k}.C
% D = nme{k}.D

t = 'mp.nme_gen_acp : '; k = 2;
t_str_match(nme{k}.name, 'gen', [t 'name']);
t_str_match(class(nme{k}), 'mp.nme_gen_acp', [t 'class']);
t_str_match(nme{k}.find_form_class(), 'mp.form_acp', [t 'formulation class']);
t_str_match(nme{k}.form_name, 'AC-polar', [t 'formulation name']);
t_str_match(nme{k}.form_tag, 'acp', [t 'formulation tag']);
t_is(nme{k}.nk, 3, 12, [t 'nk']);
t_is(nme{k}.np, 1, 12, [t 'np']);
t_is(nme{k}.nz, 1, 12, [t 'nz']);
t_ok(isempty(nme{k}.Y), [t 'Y']);
t_ok(isempty(nme{k}.L), [t 'L']);
t_ok(isempty(nme{k}.M), [t 'M']);
t_is(nme{k}.N, -speye(3), 12, [t 'N']);
t_ok(isempty(nme{k}.i), [t 'i']);
t_ok(isempty(nme{k}.s), [t 's']);
t_is(nme{k}.C, sparse([1;3;2], [1;2;3], 1, 9, 3), 12, [t 'C']);
t_is(nme{k}.D, speye(3), 12, [t 'D']);
% nme{k}
% C = nme{k}.C
% D = nme{k}.D

t = 'mp.nme_load_acp : '; k = 3;
t_str_match(nme{k}.name, 'load', [t 'name']);
t_str_match(class(nme{k}), 'mp.nme_load_acp', [t 'class']);
t_str_match(nme{k}.find_form_class(), 'mp.form_acp', [t 'formulation class']);
t_str_match(nme{k}.form_name, 'AC-polar', [t 'formulation name']);
t_str_match(nme{k}.form_tag, 'acp', [t 'formulation tag']);
t_is(nme{k}.nk, 3, 12, [t 'nk']);
t_is(nme{k}.np, 1, 12, [t 'np']);
t_is(nme{k}.nz, 0, 12, [t 'nz']);
t_ok(isempty(nme{k}.Y), [t 'Y']);
t_ok(isempty(nme{k}.L), [t 'L']);
t_ok(isempty(nme{k}.M), [t 'M']);
t_ok(isempty(nme{k}.N), [t 'N']);
t_ok(isempty(nme{k}.i), [t 'i']);
t_is(nme{k}.s, [0.9 1 1.25]' + 1j*[0.3 0.35 0.5]', 12, [t 's']);
eC = sparse([5 7 9], 1:3, 1, 9, 3);
t_is(nme{k}.C, eC, 12, [t 'C']);
t_ok(isempty(nme{k}.D), [t 'D']);
% nme{k}
% C = nme{k}.C
% D = nme{k}.D

t = 'mp.nme_branch_acp : '; k = 4;
t_str_match(nme{k}.name, 'branch', [t 'name']);
t_str_match(class(nme{k}), 'mp.nme_branch_acp', [t 'class']);
t_str_match(nme{k}.find_form_class(), 'mp.form_acp', [t 'formulation class']);
t_str_match(nme{k}.form_name, 'AC-polar', [t 'formulation name']);
t_str_match(nme{k}.form_tag, 'acp', [t 'formulation tag']);
t_is(nme{k}.nk, 9, 12, [t 'nk']);
t_is(nme{k}.np, 2, 12, [t 'np']);
t_is(nme{k}.nz, 0, 12, [t 'nz']);
nl = size(mpc.branch, 1);
stat = mpc.branch(:, BR_STATUS);    %% ones at in-service branches
tm = ones(nl, 1);                   %% default tap ratio = 1
i = find(mpc.branch(:, TAP));       %% indices of non-zero tap ratios
tm(i) = mpc.branch(i, TAP);         %% assign non-zero tap ratios
T = tm .* exp(1j*pi/180 * mpc.branch(:, SHIFT)); %% add phase shifters
Ys = stat ./ (mpc.branch(:, BR_R) + 1j * mpc.branch(:, BR_X));  %% series admittance
Bc = stat .* mpc.branch(:, BR_B);   %% line charging susceptance
Ytt = Ys + 1j*Bc/2;
Yff = Ytt ./ (T .* conj(T));
Yft = - Ys ./ conj(T);
Ytf = - Ys ./ T;
eY = sparse( ...
    [1:nl 1:nl nl+1:2*nl nl+1:2*nl]', ...
    [1:nl nl+1:2*nl 1:nl nl+1:2*nl]', ...
    [Yff; Yft; Ytf; Ytt], 2*nl, 2*nl );
t_is(nme{k}.Y, eY, 12, [t 'Y']);
t_ok(isempty(nme{k}.L), [t 'L']);
t_ok(isempty(nme{k}.M), [t 'M']);
t_ok(isempty(nme{k}.N), [t 'N']);
t_ok(isempty(nme{k}.i), [t 'i']);
t_ok(isempty(nme{k}.s), [t 's']);
t_is(nme{k}.C, sparse([1;4;5;3;6;7;8;8;9;4;5;6;6;7;8;2;9;4], [1:18]', 1, 9, 18), 12, [t 'C']);
t_ok(isempty(nme{k}.D), [t 'D']);
% nme{k}
% C = nme{k}.C
% D = nme{k}.D

t = 'ac.C';
C = ac.C;
t_is(C, sparse([1;3;2;5;7;9;1;4;5;3;6;7;8;8;9;4;5;6;6;7;8;2;9;4], [1:24]', 1, 9, 24), 12, t);

t = 'ac.D';
D = ac.D;
t_is(D, speye(3), 12, t);

t = 'S = ac.port_inj_power(x_)';
v_ = mpc.bus(:, VM) .* exp(1j * mpc.bus(:, VA) * pi/180);
z_ = (mpc.gen(:, PG) + 1j * mpc.gen(:, QG)) / mpc.baseMVA;
x_ = [v_;z_];
S  = ac.port_inj_power(x_);
eS = [-0.7195470 -0.85 -1.63 0.9 1.0 1.25 0.7195470 0.3072828 -0.5944531 0.85 0.2410613 -0.7598935 -1.63 0.8650443 -0.4096011 -0.7195470 -0.3055468 0.6089386 -0.85 -0.2401064 0.7649556 1.63 -0.8403988 0.4122642]' ...
  + 1j * [-0.2406895 0.0364902 -0.1446011 0.3 0.35 0.5 0.2406895 -0.0058585 -0.1631204 -0.0364902 0.0453679 -0.1059923 0.0227618 -0.0253242 -0.3571801 -0.2075304 -0.1368795 -0.1242746 0.0789067 -0.2440076 0.0025623 0.1446011 -0.1428198 0.2133889]';
t_is(S, eS, 6, t);
t_is(C * S, 0, 10, [t ' : C * S == 0']);

t = 'S = ac.port_inj_power(A''*x_, 0)';
nv = ac.nv;
nz = ac.nz;
np = ac.np;
A = [   C sparse(nv/2, nz);
        sparse(nz, np) D    ];
S1 = ac.port_inj_power(A'*x_, 0);
t_is(S1, eS, 6, t);

t = 'ac.port_inj_power(x_, 1, [3;1])';
S2  = ac.port_inj_power(x_, 1, [3;1]);
t_is(S2, eS([3;1]), 6, t);

t = '[S, Sva, Svm, Szr, Szi] = ac.port_inj_power(x_, 1) : ';
[S, Sva, Svm, Szr, Szi] = ac.port_inj_power(x_, 1);
t_is(S, eS, 6, [t 'S']);
t_is(size(Sva), [24 9], 12, [t 'size(Sva)']);
t_is(size(Svm), [24 9], 12, [t 'size(Svm)']);
t_is(size(Szr), [24 3], 12, [t 'size(Szr)']);
t_is(size(Szi), [24 3], 12, [t 'size(Szi)']);

t = '[S, Sva, Svm, Szr, Szi] = ac.port_inj_power(x_, 1, [3;2;1]) : ';
[S1, Sva1, Svm1, Szr1, Szi1] = ac.port_inj_power(x_, 1, [3;2;1]);
t_is(S1, eS([3;2;1]), 6, [t 'S']);
t_is(Sva1, Sva([3;2;1], :), 12, [t 'Sva']);
t_is(Svm1, Svm([3;2;1], :), 12, [t 'Svm']);
t_is(Szr1, Szr([3;2;1], :), 12, [t 'Szr']);
t_is(Szi1, Szi([3;2;1], :), 12, [t 'Szi']);

t = 'ac.elements.gen : ';
gen = ac.elements.gen;
t_str_match(gen.name, 'gen', [t 'name']);
t_str_match(class(gen), 'mp.nme_gen_acp', [t 'class']);

t = 'S = gen.port_inj_power(x_, 1)';
Sg = gen.port_inj_power(x_, 1);
eSg = -[0.7195470; 0.85; 1.63] + 1j * [-0.2406895; 0.0364902; -0.1446011];
t_is(Sg, eSg, 6, t);

t = 'S = gen.port_inj_power(x_, 1, [3;1])';
Sg31 = gen.port_inj_power(x_, 1, [3;1]);
t_is(Sg31, eSg([3;1]), 6, t);

t = 'S = gen.port_inj_power(x_, 1, 2)';
Sg2 = gen.port_inj_power(x_, 1, 2);
t_is(Sg2, eSg(2), 6, t);

t = '[S, Sva, Svm, Szr, Szi] = gen.port_inj_power(x_, 1)';
[Sg, Sgva, Sgvm, Sgzr, Sgzi] = gen.port_inj_power(x_, 1);
t_is(Sg, S(1:3), 12, t);
t_is(Sgva, Sva(1:3, :), 12, [t 'Sva']);
t_is(Sgvm, Svm(1:3, :), 12, [t 'Svm']);
t_is(Sgzr, Szr(1:3, :), 12, [t 'Szr']);
t_is(Sgzi, Szi(1:3, :), 12, [t 'Szi']);

t = '[S, Sva, Svm, Szr, Szi] = gen.port_inj_power(x_, 1, [3;1])';
[Sg, Sgva, Sgvm, Sgzr, Sgzi] = gen.port_inj_power(x_, 1, [3;1]);
t_is(Sg, S([3;1]), 12, t);
t_is(Sgva, Sva([3;1], :), 12, [t 'Sva']);
t_is(Sgvm, Svm([3;1], :), 12, [t 'Svm']);
t_is(Sgzr, Szr([3;1], :), 12, [t 'Szr']);
t_is(Sgzi, Szi([3;1], :), 12, [t 'Szi']);

t = '[S, Sva, Svm, Szr, Szi] = gen.port_inj_power(x_, 1, 2)';
[Sg, Sgva, Sgvm, Sgzr, Sgzi] = gen.port_inj_power(x_, 1, 2);
t_is(Sg, S(2), 12, t);
t_is(Sgva, Sva(2, :), 12, [t 'Sva']);
t_is(Sgvm, Svm(2, :), 12, [t 'Svm']);
t_is(Sgzr, Szr(2, :), 12, [t 'Szr']);
t_is(Sgzi, Szi(2, :), 12, [t 'Szi']);

t = 'H = ac.port_inj_power_hess(x_, lam, 1) : ';
lam = [1:np]' / 100;
H = ac.port_inj_power_hess(x_, lam, 1);
t_is(size(H), [24 24], 12, t);

t = 'H = ac.port_inj_power_hess(x_, lam(1:3), 1, [3;2;1]) : ';
H = ac.port_inj_power_hess(x_, lam(1:3), 1, [3;2;1]);
t_is(size(H), [24 24], 12, t);

t = 'I = ac.port_inj_current(x_)';
I  = ac.port_inj_current(x_);
eI = [-0.7195470 -0.8440196 -1.6311322 0.8988176 1.0183565 1.2619557 0.7195470 0.3113025 -0.5961880 0.8440196 0.2416340 -0.7720819 -1.6311322 0.8647643 -0.3982054 -0.7195470 -0.3026295 0.6023856 -0.8440196 -0.2462745 0.7663679 1.6311322 -0.8637503 0.4082444]' ...
  + 1j * [0.2406895 -0.1070623 -0.1312139 -0.3714244 -0.3440707 -0.6196286 -0.2406895 -0.0071427 0.2095040 0.1070623 -0.0371169 0.0991665 -0.1312139 0.0829015 0.4043547 0.2406895 0.1619204 0.1441792 -0.1070623 0.2449042 0.0483124 0.1312139 0.2152738 -0.2335468]';
t_is(I, eI, 6, t);
t_is(C * I, 0, 10, [t ' : C * I == 0']);

t = 'I = ac.port_inj_current(A''*x_, 0)';
I1 = ac.port_inj_current(A'*x_, 0);
t_is(I1, eI, 6, t);

t = '[I, Iva, Ivm, Izr, Izi] = ac.port_inj_current(x_, 1) : ';
[I, Iva, Ivm, Izr, Izi] = ac.port_inj_current(x_, 1);
t_is(I, eI, 6, [t 'I']);
t_is(size(Iva), [24 9], 12, [t 'size(Iva)']);
t_is(size(Ivm), [24 9], 12, [t 'size(Ivm)']);
t_is(size(Izr), [24 3], 12, [t 'size(Izr)']);
t_is(size(Izi), [24 3], 12, [t 'size(Izi)']);

t = '[I, Iva, Ivm, Izr, Izi] = ac.port_inj_current(x_, 1, [3;2;1]) : ';
[I1, Iva1, Ivm1, Izr1, Izi1] = ac.port_inj_current(x_, 1, [3;2;1]);
t_is(I1, eI([3;2;1]), 6, [t 'I']);
t_is(Iva1, Iva([3;2;1], :), 12, [t 'Iva']);
t_is(Ivm1, Ivm([3;2;1], :), 12, [t 'Ivm']);
t_is(Izr1, Izr([3;2;1], :), 12, [t 'Izr']);
t_is(Izi1, Izi([3;2;1], :), 12, [t 'Izi']);

t = 'I = gen.port_inj_current(x_, 1)';
Ig = gen.port_inj_current(x_, 1);
eIg = -[0.7195470; 0.8440197; 1.6311322] + 1j * [0.2406896; -0.1070623; -0.1312139];
t_is(Ig, eIg, 6, t);

t = 'I = gen.port_inj_current(x_, 1, [3;1])';
Ig31 = gen.port_inj_current(x_, 1, [3;1]);
t_is(Ig31, eIg([3;1]), 6, t);

t = 'I = gen.port_inj_current(x_, 1, 2)';
Ig2 = gen.port_inj_current(x_, 1, 2);
t_is(Ig2, eIg(2), 6, t);

t = '[I, Iva, Ivm, Izr, Izi] = gen.port_inj_current(x_, 1)';
[Ig, Igva, Igvm, Igzr, Igzi] = gen.port_inj_current(x_, 1);
t_is(Ig, I(1:3), 12, t);
t_is(Igva, Iva(1:3, :), 12, [t 'Iva']);
t_is(Igvm, Ivm(1:3, :), 12, [t 'Ivm']);
t_is(Igzr, Izr(1:3, :), 12, [t 'Izr']);
t_is(Igzi, Izi(1:3, :), 12, [t 'Izi']);

t = '[I, Iva, Ivm, Izr, Izi] = gen.port_inj_current(x_, 1, [3;1])';
[Ig, Igva, Igvm, Igzr, Igzi] = gen.port_inj_current(x_, 1, [3;1]);
t_is(Ig, I([3;1]), 12, t);
t_is(Igva, Iva([3;1], :), 12, [t 'Iva']);
t_is(Igvm, Ivm([3;1], :), 12, [t 'Ivm']);
t_is(Igzr, Izr([3;1], :), 12, [t 'Izr']);
t_is(Igzi, Izi([3;1], :), 12, [t 'Izi']);

t = '[I, Iva, Ivm, Izr, Izi] = gen.port_inj_current(x_, 1, 2)';
[Ig, Igva, Igvm, Igzr, Igzi] = gen.port_inj_current(x_, 1, 2);
t_is(Ig, I(2), 12, t);
t_is(Igva, Iva(2, :), 12, [t 'Iva']);
t_is(Igvm, Ivm(2, :), 12, [t 'Ivm']);
t_is(Igzr, Izr(2, :), 12, [t 'Izr']);
t_is(Igzi, Izi(2, :), 12, [t 'Izi']);

%% AC Newton power flow
t = 'mp.task_pf().run(mpc, mpopt) : ';
mpc = loadcase(casefile);
pf = mp.task_pf();
success = pf.run(mpc, mpopt);
v_ = pf.nm.soln.v;
success = pf.mm.soln.eflag > 0;
i = pf.mm.soln.output.iterations;
ev_ = [1 0.985795245 0.996534978 0.986136280 0.973075428 1.002808832 0.985586887 0.993996115 0.954862527]' + ...
 1j * [0 0.167951584 0.083174736 -0.041445908 -0.068338709 0.033715184 0.010692065 0.066005818 -0.072633402]';
t_is(v_, ev_, 8, [t 'x']);
t_is(success, 1, 12, [t 'success']);
t_is(i, 4, 12, [t 'i']);

%%-----  AC gizmo test formulation  -----
t = 'mp.net_model_acp_test() : ';
ac = mp.net_model_acp_test();
t_str_match(ac.name, 'network', [t 'name']);
t_str_match(class(ac), 'mp.net_model_acp_test', [t 'class']);
t_str_match(ac.find_form_class(), 'mp.form_acp', [t 'formulation class']);
t_str_match(ac.form_name, 'AC-polar', [t 'formulation name']);
t_str_match(ac.form_tag, 'acp', [t 'formulation tag']);
p = ac.model_params();
t_is(length(p), 6, 12, [t '# of model params']);
t_str_match(p{1}, 'Y', [t 'Y']);
t_str_match(p{2}, 'L', [t 'L']);
t_str_match(p{3}, 'M', [t 'M']);
t_str_match(p{4}, 'N', [t 'N']);
t_str_match(p{5}, 'i', [t 'i']);
t_str_match(p{6}, 's', [t 's']);
t_is(ac.nk, 0, 12, [t 'nk']);
t_is(ac.np, 0, 12, [t 'np']);
t_is(ac.nz, 0, 12, [t 'nz']);
t_is(ac.nv, 0, 12, [t 'nv']);
t_is(length(fieldnames(ac.set_types)), 7, 12, [t '# of set types']);
t_str_match(ac.set_types.node, 'NODES', [t 'set_types.node']);
t_str_match(ac.set_types.state, 'STATES', [t 'set_types.state']);
t_str_match(ac.set_types.va, 'VOLTAGE ANG VARS (va)', [t 'set_types.va']);
t_str_match(ac.set_types.vm, 'VOLTAGE MAG VARS (vm)', [t 'set_types.vm']);
t_str_match(ac.set_types.zr, 'NON-VOLTAGE VARS REAL (zr)', [t 'set_types.zr']);
t_str_match(ac.set_types.zi, 'NON-VOLTAGE VARS IMAG (zi)', [t 'set_types.zi']);
t_is(length(ac.elements), 0, 12, [t '# of element types']);

%% AC Newton power flow
t = 'mp.task_pf().run(mpc, mpopt) : ';
mpc = loadcase(casefile);
ref = find(mpc.bus(:, BUS_TYPE) == REF);
mpc.gen(ref, PG) = 1.7165997325858 * mpc.baseMVA;
mpc.gen(:, QG) = [0.2570733353840 0.0079004398259 -0.1749046999314].' * mpc.baseMVA;
% mpopt = mpoption(mpopt, 'verbose', 2);
pf = mp.task_pf();
warn_id = 'update_z:multiple_nodes';
s1 = warning('query', warn_id);
warning('off', warn_id);
% mpopt.exp.dmc_element_classes = @mp.dmce_gizmo_mpc2;
% mpopt.exp.dm_element_classes = @mp.dme_gizmo;
% mpopt.exp.network_model_class = @mp.net_model_acp_test;
% success = pf.run(mpc, mpopt);
success = pf.run(mpc, mpopt, {mp.xt_gizmo()});
warning(s1.state, warn_id);
v_ = pf.nm.soln.v;
success = pf.mm.soln.eflag > 0;
i = pf.mm.soln.output.iterations;
ev_ = [1 0.997294645 0.995662368 0.990141911 0.968051969 0.991260416 0.959080764 0.986861859 0.962323832]' + ...
 1j * [0 -0.073507764 -0.093040043 -0.093177920 -0.149175997 -0.140499450 -0.207407401 -0.144679179 -0.173550729]';
t_is(v_, ev_, 8, [t 'x']);
t_is(success, 1, 12, [t 'success']);
t_is(i, 4, 12, [t 'i']);

t = 'ac.build(dm) : ';
mpc.bus(:, VA) = angle(v_) * 180/pi;
mpc.bus(:, VM) = abs(v_);
dmc = mp.dm_converter_mpc2().modify_element_classes(@mp.dmce_gizmo_mpc2).build();
dm = mp.data_model().modify_element_classes(@mp.dme_gizmo).build(mpc, dmc);
ac = mp.net_model_acp_test().build(dm);
t_is(ac.nk, 1, 12, [t 'nk']);
t_is(ac.np, 30, 12, [t 'np']);
t_is(ac.nz, 7, 12, [t 'nz']);
t_is(ac.nv, 18, 12, [t 'nv']);
t_is(length(ac.elements), 5, 12, [t '# of element types']);

nme = ac.elements;
t = 'mp.nme_bus_acp : '; k = 1;
t_str_match(nme{k}.name, 'bus', [t 'name']);
t_str_match(class(nme{k}), 'mp.nme_bus_acp', [t 'class']);
t_str_match(nme{k}.find_form_class(), 'mp.form_acp', [t 'formulation class']);
t_str_match(nme{k}.form_name, 'AC-polar', [t 'formulation name']);
t_str_match(nme{k}.form_tag, 'acp', [t 'formulation tag']);
t_is(nme{k}.nk, 9, 12, [t 'nk']);
t_is(nme{k}.np, 0, 12, [t 'np']);
t_is(nme{k}.nz, 0, 12, [t 'nz']);
t_ok(isempty(nme{k}.Y), [t 'Y']);
t_ok(isempty(nme{k}.L), [t 'L']);
t_ok(isempty(nme{k}.M), [t 'M']);
t_ok(isempty(nme{k}.N), [t 'N']);
t_ok(isempty(nme{k}.i), [t 'i']);
t_ok(isempty(nme{k}.s), [t 's']);
t_ok(isempty(nme{k}.C), [t 'C']);
t_ok(isempty(nme{k}.D), [t 'D']);
% nme{k}
% C = nme{k}.C
% D = nme{k}.D

t = 'mp.nme_gen_acp : '; k = 2;
t_str_match(nme{k}.name, 'gen', [t 'name']);
t_str_match(class(nme{k}), 'mp.nme_gen_acp', [t 'class']);
t_str_match(nme{k}.find_form_class(), 'mp.form_acp', [t 'formulation class']);
t_str_match(nme{k}.form_name, 'AC-polar', [t 'formulation name']);
t_str_match(nme{k}.form_tag, 'acp', [t 'formulation tag']);
t_is(nme{k}.nk, 3, 12, [t 'nk']);
t_is(nme{k}.np, 1, 12, [t 'np']);
t_is(nme{k}.nz, 1, 12, [t 'nz']);
t_ok(isempty(nme{k}.Y), [t 'Y']);
t_ok(isempty(nme{k}.L), [t 'L']);
t_ok(isempty(nme{k}.M), [t 'M']);
t_is(nme{k}.N, -speye(3), 12, [t 'N']);
t_ok(isempty(nme{k}.i), [t 'i']);
t_ok(isempty(nme{k}.s), [t 's']);
t_is(nme{k}.C, sparse([1;3;2], [1;2;3], 1, 9, 3), 12, [t 'C']);
t_is(nme{k}.D, sparse(1:3, 1:3, 1, 7, 3), 12, [t 'D']);
% nme{k}
% C = nme{k}.C
% D = nme{k}.D

t = 'mp.nme_load_acp : '; k = 3;
t_str_match(nme{k}.name, 'load', [t 'name']);
t_str_match(class(nme{k}), 'mp.nme_load_acp', [t 'class']);
t_str_match(nme{k}.find_form_class(), 'mp.form_acp', [t 'formulation class']);
t_str_match(nme{k}.form_name, 'AC-polar', [t 'formulation name']);
t_str_match(nme{k}.form_tag, 'acp', [t 'formulation tag']);
t_is(nme{k}.nk, 3, 12, [t 'nk']);
t_is(nme{k}.np, 1, 12, [t 'np']);
t_is(nme{k}.nz, 0, 12, [t 'nz']);
t_ok(isempty(nme{k}.Y), [t 'Y']);
t_ok(isempty(nme{k}.L), [t 'L']);
t_ok(isempty(nme{k}.M), [t 'M']);
t_ok(isempty(nme{k}.N), [t 'N']);
t_ok(isempty(nme{k}.i), [t 'i']);
t_is(nme{k}.s, [0.9 1 1.25]' + 1j*[0.3 0.35 0.5]', 12, [t 's']);
eC = sparse([5 7 9], 1:3, 1, 9, 3);
t_is(nme{k}.C, eC, 12, [t 'C']);
t_ok(isempty(nme{k}.D), [t 'D']);
% nme{k}
% C = nme{k}.C
% D = nme{k}.D

t = 'mp.nme_branch_acp : '; k = 4;
t_str_match(nme{k}.name, 'branch', [t 'name']);
t_str_match(class(nme{k}), 'mp.nme_branch_acp', [t 'class']);
t_str_match(nme{k}.find_form_class(), 'mp.form_acp', [t 'formulation class']);
t_str_match(nme{k}.form_name, 'AC-polar', [t 'formulation name']);
t_str_match(nme{k}.form_tag, 'acp', [t 'formulation tag']);
t_is(nme{k}.nk, 9, 12, [t 'nk']);
t_is(nme{k}.np, 2, 12, [t 'np']);
t_is(nme{k}.nz, 0, 12, [t 'nz']);
nl = size(mpc.branch, 1);
stat = mpc.branch(:, BR_STATUS);    %% ones at in-service branches
tm = ones(nl, 1);                   %% default tap ratio = 1
i = find(mpc.branch(:, TAP));       %% indices of non-zero tap ratios
tm(i) = mpc.branch(i, TAP);         %% assign non-zero tap ratios
T = tm .* exp(1j*pi/180 * mpc.branch(:, SHIFT));    %% add phase shifters
Ys = stat ./ (mpc.branch(:, BR_R) + 1j * mpc.branch(:, BR_X));  %% series admittance
Bc = stat .* mpc.branch(:, BR_B);   %% line charging susceptance
Ytt = Ys + 1j*Bc/2;
Yff = Ytt ./ (T .* conj(T));
Yft = - Ys ./ conj(T);
Ytf = - Ys ./ T;
eY = sparse( ...
    [1:nl 1:nl nl+1:2*nl nl+1:2*nl]', ...
    [1:nl nl+1:2*nl 1:nl nl+1:2*nl]', ...
    [Yff; Yft; Ytf; Ytt], 2*nl, 2*nl );
t_is(nme{k}.Y, eY, 12, [t 'Y']);
t_ok(isempty(nme{k}.L), [t 'L']);
t_ok(isempty(nme{k}.M), [t 'M']);
t_ok(isempty(nme{k}.N), [t 'N']);
t_ok(isempty(nme{k}.i), [t 'i']);
t_ok(isempty(nme{k}.s), [t 's']);
t_is(nme{k}.C, sparse([1;4;5;3;6;7;8;8;9;4;5;6;6;7;8;2;9;4], [1:18]', 1, 9, 18), 12, [t 'C']);
t_ok(isempty(nme{k}.D), [t 'D']);
% nme{k}
% C = nme{k}.C
% D = nme{k}.D

t = 'mp.nme_gizmo_acp : '; k = 5;
t_str_match(nme{k}.name, 'gizmo', [t 'name']);
t_str_match(class(nme{k}), 'mp.nme_gizmo_acp', [t 'class']);
t_str_match(nme{k}.find_form_class(), 'mp.form_acp', [t 'formulation class']);
t_str_match(nme{k}.form_name, 'AC-polar', [t 'formulation name']);
t_str_match(nme{k}.form_tag, 'acp', [t 'formulation tag']);
t_is(nme{k}.nk, 2, 12, [t 'nk']);
t_is(nme{k}.np, 3, 12, [t 'np']);
t_is(nme{k}.nz, 2, 12, [t 'nz']);
nk = size(mpc.gizmo, 1);
y1 = mpc.gizmo(:,  4) + 1j * mpc.gizmo(:,  5);
y2 = mpc.gizmo(:,  6) + 1j * mpc.gizmo(:,  7);
ll = mpc.gizmo(:,  8) + 1j * mpc.gizmo(:,  9);
ii = mpc.gizmo(:, 10) + 1j * mpc.gizmo(:, 11);
m1 = mpc.gizmo(:, 12) + 1j * mpc.gizmo(:, 13);
m2 = mpc.gizmo(:, 14) + 1j * mpc.gizmo(:, 15);
nn = mpc.gizmo(:, 16) + 1j * mpc.gizmo(:, 17);
ss = mpc.gizmo(:, 18) + 1j * mpc.gizmo(:, 19);
zz = zeros(nk, 1);
j1 = (1:nk); j2 = nk+j1; j3 = nk+j2;
eY = sparse( ...
    [j1 j1 j1 j2 j2 j2 j3 j3 j3]', ...
    [j1 j2 j3 j1 j2 j3 j1 j2 j3]', ...
    [y1; zz; -y1; zz; y2; zz; -y1; zz; y1], 3*nk, 3*nk );
eL = sparse( ...
    [j1 j1 j2 j2 j3 j3 ]', ...
    [j1 j2 j1 j2 j1 j2 ]', ...
    [zz; ll; zz; -ll; zz; zz], 3*nk, 2*nk );
ei = [-ii; ii; zz];
eM = sparse( ...
    [j1 j1 j1 j2 j2 j2 j3 j3 j3]', ...
    [j1 j2 j3 j1 j2 j3 j1 j2 j3]', ...
    [m1; -m1; zz; -m1; m1; zz; zz; zz; m2], 3*nk, 3*nk );
eN = sparse( ...
    [j1 j1 j2 j2 j3 j3 ]', ...
    [j1 j2 j1 j2 j1 j2 ]', ...
    [zz; zz; nn; zz; -nn; zz], 3*nk, 2*nk );
es = [zz; -ss; ss];
t_is(nme{k}.Y, eY, 12, [t 'Y']);
t_is(nme{k}.L, eL, 12, [t 'L']);
t_is(nme{k}.M, eM, 12, [t 'M']);
t_is(nme{k}.N, eN, 12, [t 'N']);
t_is(nme{k}.i, ei, 12, [t 'i']);
t_is(nme{k}.s, es, 12, [t 's']);
t_is(nme{k}.C, sparse([1;5;2;7;3;9], [1:6]', 1, 9, 6), 12, [t 'C']);
t_is(nme{k}.D, sparse(4:7, 1:4, 1, 7, 4), 12, [t 'D']);
% nme{k}
% C = nme{k}.C
% D = nme{k}.D

t = 'ac.C';
C = ac.C;
t_is(C, sparse([1;3;2;5;7;9;1;4;5;3;6;7;8;8;9;4;5;6;6;7;8;2;9;4;1;5;2;7;3;9], [1:30]', 1, 9, 30), 12, t);

t = 'ac.D';
D = ac.D;
t_is(D, speye(7), 12, t);

t = 'S = ac.port_inj_power(x_)';
v_ = mpc.bus(:, VM) .* exp(1j * mpc.bus(:, VA) * pi/180);
z_ = [ (mpc.gen(:, PG) + 1j * mpc.gen(:, QG)) / mpc.baseMVA;
      mpc.gizmo(:, 20) + 1j * mpc.gizmo(:, 21);
      mpc.gizmo(:, 22) + 1j * mpc.gizmo(:, 23)  ];
x_ = [v_;z_];
S  = ac.port_inj_power(x_);
eS = [-1.7165997 -0.85 -1.63 0.9 1 1.25 1.6176722 0.6367558 -0.0934284 0.8133635 0.7191480 -0.9250656 -1.1479321 0.2152553 -0.9710411 -1.6176722 -0.6297178 0.0942154 -0.8133635 -0.7127809 0.9326768 1.1479321 -0.2135503 0.9809164 0.0989275 -0.1768537 0.4820678 0.6378466 0.0366364 -0.0654085]' ...
  + 1j * [-0.2570733 -0.0079004 0.1749047 0.3 0.35 0.5 0.1711473 -0.0147713 -0.2748400 -0.0005591 0.0335496 -0.1518847 0.0000208 -0.0705294 -0.1207424 -0.0187288 -0.1010674 -0.0728762 0.0393266 -0.1849791 0.0705086 0.0827665 -0.2193974 0.0335002 0.0859259 0.0759074 -0.2576712 -0.0131360 0.0084596 -0.1598600]';
t_is(S, eS, 6, t);
t_is(C * S, 0, 10, [t ' : C * S == 0']);

t = 'S = ac.port_inj_power(A''*x_, 0)';
nv = ac.nv;
nz = ac.nz;
np = ac.np;
A = [   C sparse(nv/2, nz);
        sparse(nz, np) D    ];
S1 = ac.port_inj_power(A'*x_, 0);
t_is(S1, eS, 6, t);

t = 'S = ac.port_inj_power(x_, 1, [3;1])';
S2  = ac.port_inj_power(x_, 1, [3;1]);
t_is(S2, eS([3;1]), 6, t);

t = '[S, Sva, Svm, Szr, Szi] = ac.port_inj_power(x_, 1) : ';
[S, Sva, Svm, Szr, Szi] = ac.port_inj_power(x_, 1);
t_is(S, eS, 6, [t 'S']);
t_is(size(Sva), [30 9], 12, [t 'size(Sva)']);
t_is(size(Svm), [30 9], 12, [t 'size(Svm)']);
t_is(size(Szr), [30 7], 12, [t 'size(Szr)']);
t_is(size(Szi), [30 7], 12, [t 'size(Szi)']);

t = '[S, Sva, Svm, Szr, Szi] = ac.port_inj_power(x_, 1, [3;2;1]) : ';
[S1, Sva1, Svm1, Szr1, Szi1] = ac.port_inj_power(x_, 1, [3;2;1]);
t_is(S1, eS([3;2;1]), 6, [t 'S']);
t_is(Sva1, Sva([3;2;1], :), 12, [t 'Sva']);
t_is(Svm1, Svm([3;2;1], :), 12, [t 'Svm']);
t_is(Szr1, Szr([3;2;1], :), 12, [t 'Szr']);
t_is(Szi1, Szi([3;2;1], :), 12, [t 'Szi']);

t = 'ac.elements.gen : ';
gen = ac.elements.gen;
t_str_match(gen.name, 'gen', [t 'name']);
t_str_match(class(gen), 'mp.nme_gen_acp', [t 'class']);

t = 'S = gen.port_inj_power(x_, 1)';
Sg = gen.port_inj_power(x_, 1);
eSg = [-1.7165997; -0.85; -1.63] + 1j * [-0.2570733; -0.0079004; 0.1749047];
t_is(Sg, eSg, 6, t);

t = 'S = gen.port_inj_power(x_, 1, [3;1])';
Sg31 = gen.port_inj_power(x_, 1, [3;1]);
t_is(Sg31, eSg([3;1]), 6, t);

t = 'S = gen.port_inj_power(x_, 1, 2)';
Sg2 = gen.port_inj_power(x_, 1, 2);
t_is(Sg2, eSg(2), 6, t);

t = '[S, Sva, Svm, Szr, Szi] = gen.port_inj_power(x_, 1)';
[Sg, Sgva, Sgvm, Sgzr, Sgzi] = gen.port_inj_power(x_, 1);
t_is(Sg, S(1:3), 12, t);
t_is(Sgva, Sva(1:3, :), 12, [t 'Sva']);
t_is(Sgvm, Svm(1:3, :), 12, [t 'Svm']);
t_is(Sgzr, Szr(1:3, :), 12, [t 'Szr']);
t_is(Sgzi, Szi(1:3, :), 12, [t 'Szi']);

t = '[S, Sva, Svm, Szr, Szi] = gen.port_inj_power(x_, 1, [3;1])';
[Sg, Sgva, Sgvm, Sgzr, Sgzi] = gen.port_inj_power(x_, 1, [3;1]);
t_is(Sg, S([3;1]), 12, t);
t_is(Sgva, Sva([3;1], :), 12, [t 'Sva']);
t_is(Sgvm, Svm([3;1], :), 12, [t 'Svm']);
t_is(Sgzr, Szr([3;1], :), 12, [t 'Szr']);
t_is(Sgzi, Szi([3;1], :), 12, [t 'Szi']);

t = '[S, Sva, Svm, Szr, Szi] = gen.port_inj_power(x_, 1, 2)';
[Sg, Sgva, Sgvm, Sgzr, Sgzi] = gen.port_inj_power(x_, 1, 2);
t_is(Sg, S(2), 12, t);
t_is(Sgva, Sva(2, :), 12, [t 'Sva']);
t_is(Sgvm, Svm(2, :), 12, [t 'Svm']);
t_is(Sgzr, Szr(2, :), 12, [t 'Szr']);
t_is(Sgzi, Szi(2, :), 12, [t 'Szi']);

t = 'H = ac.port_inj_power_hess(x_, lam, 1) : ';
lam = [1:np]' / 100;
H = ac.port_inj_power_hess(x_, lam, 1);
t_is(size(H), [32 32], 12, t);

t = 'H = ac.port_inj_power_hess(x_, lam(1:3), 1, [3;2;1]) : ';
H = ac.port_inj_power_hess(x_, lam(1:3), 1, [3;2;1]);
t_is(size(H), [32 32], 12, t);

t = 'I = ac.port_inj_current(x_)';
I  = ac.port_inj_current(x_);
eI = [-1.7165997 -0.8455779 -1.6384471 0.8614893 0.9206882 1.1672710 1.6176722 0.6388419 -0.0515376 0.8098875 0.7064979 -0.8887234 -1.1387426 0.2237887 -0.9553581 -1.6176722 -0.6196959 0.1033895 -0.8098875 -0.6701417 0.9149538 1.1387426 -0.1750995 0.9788303 0.0989275 -0.1902557 0.4997044 0.6381768 0.0356904 -0.0368132]' ...
  + 1j * [0.2570733 0.0869502 -0.0546138 -0.4426554 -0.5640375 -0.7300876 -0.1711473 -0.0452001 0.2918523 -0.0751186 -0.1339831 0.3505571 0.1669245 0.0386598 0.2977641 0.1711473 0.1998974 0.0588645 0.0751186 0.3377938 -0.2055844 -0.1669245 0.2595655 -0.1259471 -0.0859259 -0.0490943 0.2215384 -0.1243133 -0.0118315 0.1727578]';
t_is(I, eI, 6, t);
t_is(C * I, 0, 10, [t ' : C * I == 0']);

t = 'I = ac.port_inj_current(A''*x_, 0)';
I1 = ac.port_inj_current(A'*x_, 0);
t_is(I1, eI, 6, t);

t = '[I, Iva, Ivm, Izr, Izi] = ac.port_inj_current(x_, 1) : ';
[I, Iva, Ivm, Izr, Izi] = ac.port_inj_current(x_, 1);
t_is(I, eI, 6, [t 'I']);
t_is(size(Iva), [30 9], 12, [t 'size(Iva)']);
t_is(size(Ivm), [30 9], 12, [t 'size(Ivm)']);
t_is(size(Izr), [30 7], 12, [t 'size(Izr)']);
t_is(size(Izi), [30 7], 12, [t 'size(Izi)']);

t = '[I, Iva, Ivm, Izr, Izi] = ac.port_inj_current(x_, 1, [3;2;1]) : ';
[I1, Iva1, Ivm1, Izr1, Izi1] = ac.port_inj_current(x_, 1, [3;2;1]);
t_is(I1, eI([3;2;1]), 6, [t 'I']);
t_is(Iva1, Iva([3;2;1], :), 12, [t 'Iva']);
t_is(Ivm1, Ivm([3;2;1], :), 12, [t 'Ivm']);
t_is(Izr1, Izr([3;2;1], :), 12, [t 'Izr']);
t_is(Izi1, Izi([3;2;1], :), 12, [t 'Izi']);

t = 'I = gen.port_inj_current(x_, 1)';
Ig = gen.port_inj_current(x_, 1);
eIg = -[1.7165997; 0.8455780; 1.6384471] + 1j * [0.2570733; 0.0869502; -0.0546139];
t_is(Ig, eIg, 6, t);

t = 'I = gen.port_inj_current(x_, 1, [3;1])';
Ig31 = gen.port_inj_current(x_, 1, [3;1]);
t_is(Ig31, eIg([3;1]), 6, t);

t = 'I = gen.port_inj_current(x_, 1, 2)';
Ig2 = gen.port_inj_current(x_, 1, 2);
t_is(Ig2, eIg(2), 6, t);

t = '[I, Iva, Ivm, Izr, Izi] = gen.port_inj_current(x_, 1)';
[Ig, Igva, Igvm, Igzr, Igzi] = gen.port_inj_current(x_, 1);
t_is(Ig, I(1:3), 12, t);
t_is(Igva, Iva(1:3, :), 12, [t 'Iva']);
t_is(Igvm, Ivm(1:3, :), 12, [t 'Ivm']);
t_is(Igzr, Izr(1:3, :), 12, [t 'Izr']);
t_is(Igzi, Izi(1:3, :), 12, [t 'Izi']);

t = '[I, Iva, Ivm, Izr, Izi] = gen.port_inj_current(x_, 1, [3;1])';
[Ig, Igva, Igvm, Igzr, Igzi] = gen.port_inj_current(x_, 1, [3;1]);
t_is(Ig, I([3;1]), 12, t);
t_is(Igva, Iva([3;1], :), 12, [t 'Iva']);
t_is(Igvm, Ivm([3;1], :), 12, [t 'Ivm']);
t_is(Igzr, Izr([3;1], :), 12, [t 'Izr']);
t_is(Igzi, Izi([3;1], :), 12, [t 'Izi']);

t = '[I, Iva, Ivm, Izr, Izi] = gen.port_inj_current(x_, 1, 2)';
[Ig, Igva, Igvm, Igzr, Igzi] = gen.port_inj_current(x_, 1, 2);
t_is(Ig, I(2), 12, t);
t_is(Igva, Iva(2, :), 12, [t 'Iva']);
t_is(Igvm, Ivm(2, :), 12, [t 'Ivm']);
t_is(Igzr, Izr(2, :), 12, [t 'Izr']);
t_is(Igzi, Izi(2, :), 12, [t 'Izi']);

t = 'ac.set_type_idx_map(''node'', ..., 1) : ';
n = ac.set_type_idx_map('node', (1:ac.getN('node')), dm, 1);
e = struct( 'name', 'bus', ...
            'idx', [], ...
            'i', [1:9], ...
            'e', [1:9], ...
            'ID', [1 2 30 4 5 6 7 8 9], ...
            'j', [1:9] );
t_ok(isequal(n, e), t);

t = 'ac.set_type_idx_map(''port'', ..., 1) : ';
p = ac.set_type_idx_map('port', (1:ac.getN('port')), dm, 1);
e = struct( 'name', {'branch', 'branch', 'gen', 'gizmo', 'gizmo', 'gizmo', 'load'}, ...
            'idx', {{1}, {2}, [], {1}, {2}, {3}, []}, ...
            'i', {[1:9], [1:9], [1:3], [1 2], [1 2], [1 2], [1:3]}, ...
            'e', {[1:9], [1:9], [1:3], [1 2], [1 2], [1 2], [1:3]}, ...
            'ID', {[1:9], [1:9], [1:3], [1 2], [1 2], [1 2], [1:3]}, ...
            'j', {[7:15], [16:24], [1:3], [25 26], [27 28], [29 30], [4:6]}   );
t_ok(isequal(p, e), t);

t = 'ac.set_type_idx_map(''state'', ..., 1) : ';
s = ac.set_type_idx_map('state', (1:ac.getN('state')), dm, 1);
e = struct( 'name', {'gen', 'gizmo', 'gizmo'}, ...
            'idx', {[], {1}, {2}}, ...
            'i', {[1:3], [1 2], [1 2]}, ...
            'e', {[1:3], [1 2], [1 2]}, ...
            'ID', {[1:3], [1 2], [1 2]}, ...
            'j', {[1:3], [4 5], [6 7]}   );
t_ok(isequal(s, e), t);

t = 'ac.set_type_label(''node'', ...) : ';
i = [1 3 9];
e = {'bus 1', 'bus 30', 'bus 9'};
for k = 1:length(e)
    t_ok(isequal(ac.set_type_label('node', i(k), dm), e{k}), [t e{k}]);
end
t_ok(isequal(ac.set_type_label('node', i, dm), e), [t e{k}]);
% ac.set_type_label('node', (1:ac.getN('node'))', dm)
% for k = 1:ac.getN('node')
%     fprintf('node %d : %s\n', k, ac.set_type_label('node', k, dm));
% end

t = 'ac.set_type_label(''port'', ...) : ';
i = [2 6 7 22 29];
e = {'gen 2', 'load 3', 'branch(1) 1', 'branch(2) 7', 'gizmo(3) 1'};
for k = 1:length(e)
    t_ok(isequal(ac.set_type_label('port', i(k), dm), e{k}), [t e{k}]);
end
t_ok(isequal(ac.set_type_label('port', i, dm), e), [t e{k}]);
% ac.set_type_label('port', (1:ac.getN('port'))', dm)
% for k = 1:ac.getN('port')
%     fprintf('port %d : %s\n', k, ac.set_type_label('port', k, dm));
% end

t = 'ac.set_type_label(''state'', ...) : ';
i = [2 3 5 6];
e = {'gen 2', 'gen 3', 'gizmo(1) 2', 'gizmo(2) 1'};
for k = 1:length(e)
    t_ok(isequal(ac.set_type_label('state', i(k), dm), e{k}), [t e{k}]);
end
t_ok(isequal(ac.set_type_label('state', i, dm), e), [t e{k}]);
% ac.set_type_label('state', (1:ac.getN('state'))', dm)
% for k = 1:ac.getN('state')
%     fprintf('state %d : %s\n', k, ac.set_type_label('state', k, dm));
% end

% disp(dc)
% disp(ac)

t_end;

if nargout
    if nargin < 2
        out_ac = 0;
    end
    if out_ac
        obj = ac;
    else
        obj = dc;
    end
end
