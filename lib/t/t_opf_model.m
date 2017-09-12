function t_opf_model(quiet)
%T_OPF_MODEL Tests for OPF_MODEL.

%   MATPOWER
%   Copyright (c) 2012-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

num_tests = 591;

t_begin(num_tests, quiet);

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%%-----  opt_model  -----
t = 'constructor';
om = opf_model;
t_ok(isa(om, 'opf_model'), t);
t_ok(isa(om, 'opt_model'), t);

%%-----  add_var  -----
t = 'add_var';
vN = 0;
vNS = 0;
t_ok(om.getN('var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(om.get('var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'om.add_var(''Va'', 4)';
nVa = 4;
om.add_var('Va', nVa);
vNS = vNS + 1; vN = vN + nVa;
t_ok(om.getN('var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(om.get('var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'om.add_var(''Pg'', 3, Pg0, Pgmin, Pgmax)';
nPg = 3;
om.add_var('Pg', nPg, [2;4;6], [1;2;3], [10;20;30]);
vNS = vNS + 1; vN = vN + nPg;
t_ok(om.getN('var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(om.get('var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'om.add_var(''Vm1'', 5, V0, Vmin, Vmax, ''I'')';
V0 = [1;1;1;1;1];
Vmin = zeros(5, 1);
Vmax = 1 + 0.01*(1:5)';
vt = 'I';
om.add_var('Vm1', 5, V0, Vmin, Vmax, vt);
vNS = vNS + 1; vN = vN + 5;
t_ok(om.getN('var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(om.get('var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'om.add_var(''Vm2'', 5, V0, Vmin, Vmax, ''CIBIC'')';
nVm2 = 5;
vt = 'CIBIC';
om.add_var('Vm2', nVm2, V0, Vmin, Vmax, vt);
vNS = vNS + 1; vN = vN + nVm2;
t_ok(om.getN('var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(om.get('var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'om.init_indexed_name(''var'', ''x'', dims)';
om.init_indexed_name('var', 'x', {2,2});
t_ok(om.getN('var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(om.get('var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'om.add_var(''x'', {1,1}, 2)';
om.add_var('x', {1,1}, 2);
vNS = vNS + 1; vN = vN + 2;
t_ok(om.getN('var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(om.get('var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'om.add_var(''x'', {1,2}, 2, x0(1,2))';
om.add_var('x', {1,2}, 2, [-1;-2]);
vNS = vNS + 1; vN = vN + 2;
t_ok(om.getN('var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(om.get('var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'om.add_var(''x'', {2,1}, 3)';
om.add_var('x', {2,1}, 3);
vNS = vNS + 1; vN = vN + 3;
t_ok(om.getN('var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(om.get('var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'om.add_var(''x'', {2,2}, 2, x0(2,2), xmin(2,2), xmax(2,2))';
om.add_var('x', {2,2}, 2, [1;0],[0;-1],[2;1]);
vNS = vNS + 1; vN = vN + 2;
t_ok(om.getN('var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(om.get('var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'om.init_indexed_name(''var'', ''y'', {2,3,4})';
om.init_indexed_name('var', 'y', {2,3,4});
vt0 = {'C', 'I', 'B'};
for i = 1:2
    for j = 1:3
        for k = 1:4
            n = i+j+k;
            if i == 1
                vt = vt0{j};
            else
                vt = char(vt0{j} * ones(1, n));
                vt(j+1) = vt0{1+rem(j,3)};
            end
%             fprintf('%d %d %d : %s\n', i, j, k, vt);
            t = sprintf('om.add_var(''y'', {%d,%d,%d}, y0, ymin, ymax, vt)', i,j,k);
            om.add_var('y', {i,j,k}, n, 10*(n:-1:1)', -1*(n:-1:1)', 100+(n:-1:1)', vt);
            vNS = vNS + 1; vN = vN + n;
            t_ok(om.getN('var') == vN, sprintf('%s : var.N  = %d', t, vN));
            t_ok(om.get('var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));
        end
    end
end

%%-----  getN  -----
t = 'om.getN(''var'', ''Pg'') == 3';
t_ok(om.getN('var', 'Pg') == 3, t);

t = 'size(om.getN(''var'', ''x'')) == [2 2]';
t_is(size(om.getN('var', 'x')), [2,2], 14, t);

t = 'om.getN(''var'', ''x'')(1,2) == 2';
N = om.getN('var', 'x');
t_is(N(1,2), 2, 14, t);

t = 'om.getN(''var'', ''x'', {2,1}) == 3';
t_is(om.getN('var', 'x', {2,1}), 3, 14, t);

t = 'om.getN(''var'', ''y'', {2,1,3}) == 6';
t_is(om.getN('var', 'y', {2,1,3}), 6, 14, t);

t = 'om.getN(''var'')';
t_is(om.getN('var'), vN, 14, t);

%%-----  get_idx  -----
t = 'get_idx : var';
vv = om.get_idx();
t_is([vv.i1.Pg vv.iN.Pg vv.N.Pg], [5 7 3], 14, [t ' : Pg']);
t_is(size(vv.i1.x), [2, 2], 14, [t ' : size(vv.i1.x)']);
t_is([vv.i1.x(2,1) vv.iN.x(2,1) vv.N.x(2,1)], [22 24 3], 14, [t ' : x(2,1)']);
t_is(size(vv.i1.y), [2, 3, 4], 14, [t ' : size(vv.i1.y)']);
t_is([vv.i1.y(2,2,4) vv.iN.y(2,2,4) vv.N.y(2,2,4)], [133 140 8], 14, [t ' : y(2,2,4)']);

t = 'get_idx(''var'')';
vv = om.get_idx('var');
t_is([vv.i1.Pg vv.iN.Pg vv.N.Pg], [5 7 3], 14, [t ' : Pg']);
t_is(size(vv.i1.x), [2, 2], 14, [t ' : size(vv.i1.x)']);
t_is([vv.i1.x(2,1) vv.iN.x(2,1) vv.N.x(2,1)], [22 24 3], 14, [t ' : x(2,1)']);
t_is(size(vv.i1.y), [2, 3, 4], 14, [t ' : size(vv.i1.y)']);
t_is([vv.i1.y(2,2,4) vv.iN.y(2,2,4) vv.N.y(2,2,4)], [133 140 8], 14, [t ' : y(2,2,4)']);

%%-----  params_var  -----
t = 'om.params_var(''Va'')';
[v0, vl, vu] = om.params_var('Va');
t_ok(~any(v0), [t ' : v0']);
t_ok(all(isinf(vl) & vl < 0), [t ' : vl']);
t_ok(all(isinf(vu) & vu > 0), [t ' : vu']);

t = 'om.params_var(''Pg'')';
[v0, vl, vu] = om.params_var('Pg');
t_is(v0, [2;4;6], 14, [t ' : v0']);
t_is(vl, [1;2;3], 14, [t ' : vl']);
t_is(vu, [10;20;30], 14, [t ' : vu']);

t = 'om.params_var(''Vm1'')';
[v0, vl, vu, vt] = om.params_var('Vm1');
t_is(double(vt), double('I'), 14, [t ' : vt']);

t = 'om.params_var(''Vm2'')';
[v0, vl, vu, vt] = om.params_var('Vm2');
t_is(double(vt), double('CIBIC'), 14, [t ' : vt']);

t = 'om.params_var(''x'')';
[v0, vl, vu, vt] = om.params_var('x');
t_is(size(v0), [2,2], 14, [t ' : size(v0)']);
t_is(v0{2,2}, [1;0], 14, [t ' : v0{2,2}']);
t_is(vl{2,2}, [0;-1], 14, [t ' : vl{2,2}']);
t_is(vu{2,2}, [2;1], 14, [t ' : vu{2,2}']);
t_is(double(vt{2,2}), double('C'), 14, [t ' : vt{2,2}']);

for i = 1:2
    for j = 1:3
        for k = 1:4
            n = i+j+k;
            if i == 1
                vt = vt0{j};
            else
                vt = char(vt0{j} * ones(1, n));
                vt(j+1) = vt0{1+rem(j,3)};
            end
            t = sprintf('om.params_var(''y'', {%d,%d,%d})', i, j, k);
            [v0, vl, vu, gvt] = om.params_var('y', {i,j,k});
            t_is(v0, 10*(n:-1:1)', 14, [t ' : v0']);
            t_is(vl, -1*(n:-1:1)', 14, [t ' : vl']);
            t_is(vu, 100+(n:-1:1)', 14, [t ' : vu']);
            t_is(gvt, vt, 14, [t ' : vt']);
        end
    end
end

t = 'om.params_var()';
[v0, vl, vu, vt] = om.params_var();
t_ok(length(v0) == om.getN('var'), [t ' : length(v0)']);
t_ok(length(vl) == om.getN('var'), [t ' : length(vl)']);
t_ok(length(vu) == om.getN('var'), [t ' : length(vu)']);
t_is(v0(vv.i1.x(2,2):vv.iN.x(2,2)), [1;0], 14, [t ' : v0(vv.i1.x(2,2):vv.iN.x(2,2))']);
t_is(vl(vv.i1.x(2,2):vv.iN.x(2,2)), [0;-1], 14, [t ' : vl(vv.i1.x(2,2):vv.iN.x(2,2))']);
t_is(vu(vv.i1.x(2,2):vv.iN.x(2,2)), [2;1], 14, [t ' : vu(vv.i1.x(2,2):vv.iN.x(2,2))']);
t_is(vt(vv.i1.x(2,2):vv.iN.x(2,2)), 'C', 14, [t ' : vt(vv.i1.x(2,2):vv.iN.x(2,2))']);
n = 8;
t_is(v0(vv.i1.y(2,2,4):vv.iN.y(2,2,4)), 10*(n:-1:1)', 14, [t ' : v0(vv.i1.y(2,2,4):vv.iN.y(2,2,4))']);
t_is(vl(vv.i1.y(2,2,4):vv.iN.y(2,2,4)), -1*(n:-1:1)', 14, [t ' : vl(vv.i1.y(2,2,4):vv.iN.y(2,2,4))']);
t_is(vu(vv.i1.y(2,2,4):vv.iN.y(2,2,4)), 100+(n:-1:1)', 14, [t ' : vu(vv.i1.y(2,2,4):vv.iN.y(2,2,4))']);
t_is(vt(vv.i1.y(2,2,4):vv.iN.y(2,2,4)), 'IIBIIIII', 14, [t ' : vt(vv.i1.y(2,2,4):vv.iN.y(2,2,4))']);
vt0 = 'CCCCCCCIIIIICIBICCCCCCCCCCCCCCCCCCCCCCCCCCCCIIIIIIIIIIIIIIIIIIIIIIBBBBBBBBBBBBBBBBBBBBBBBBBBCICCCICCCCICCCCCICCCCCIIBIIIIBIIIIIBIIIIIIBIIIIIBBBCBBBBBCBBBBBBCBBBBBBBCBBBBB';
t_is(vt, vt0, 14, [t ' : vt']);

%%-----  varsets_len  -----
t = 'om.varsets_len(vs) : ';
vs = om.varsets_cell2struct({'Pg'});
t_is(om.varsets_len(vs), 3, 14, [t '{''Pg''}']);

vs = om.varsets_cell2struct({'Pg', 'Va'});
t_is(om.varsets_len(vs), 7, 14, [t '{''Pg'', ''Va''}']);

vs = struct('name', 'x', 'idx', {{1,1},{2,1}});
t_is(om.varsets_len(vs), 5, 14, [t '''x'', {{1,1},{2,1}}']);

vs = om.varsets_cell2struct({'x'});
t_is(om.varsets_len(vs), 9, 14, [t '{''x''}']);

vs = om.varsets_cell2struct({'x', 'y', 'Pg'});
t_is(om.varsets_len(vs), 156, 14, [t '{''x'', ''y'', ''Pg''}']);

vs = om.varsets_cell2struct({});
t_is(om.varsets_len(vs), om.var.N, 14, [t '<all>']);

%%-----  varsets_idx  -----
t = 'om.varsets_idx(vs) : ';
vv = om.get_idx('var');
vs = om.varsets_cell2struct({'Pg'});
t_is(om.varsets_idx(vs), [vv.i1.Pg:vv.iN.Pg], 14, [t '{''Pg''}']);

vs = om.varsets_cell2struct({'Pg', 'Va'});
t_is(om.varsets_idx(vs), [vv.i1.Pg:vv.iN.Pg vv.i1.Va:vv.iN.Va], 14, [t '{''Pg'', ''Va''}']);

vs = struct('name', 'x', 'idx', {{1,1},{2,1}});
t_is(om.varsets_idx(vs), [vv.i1.x(1,1):vv.iN.x(1,1) vv.i1.x(2,1):vv.iN.x(2,1)], 14, [t '''x'', {{1,1},{2,1}}']);

vs = om.varsets_cell2struct({'x'});
t_is(om.varsets_idx(vs), [vv.i1.x(1,1):vv.iN.x(1,1) vv.i1.x(1,2):vv.iN.x(1,2) vv.i1.x(2,1):vv.iN.x(2,1) vv.i1.x(2,2):vv.iN.x(2,2)], 14, [t '{''x''}']);

vs = om.varsets_cell2struct({});
t_is(om.varsets_idx(vs), 1:om.var.N, 14, [t '<all>']);

%%-----  varsets_x  -----
t = 'varsets_x(x, vs) : ';
x = (1:om.var.N)';
vs = om.varsets_cell2struct({'Pg'});
xx = om.varsets_x(x, vs);
t_is(length(xx), 1, 14, [t '{''Pg''} : length']);
t_is(xx{1}, [vv.i1.Pg:vv.iN.Pg]', 14, [t '{''Pg''} : 1']);

vs = om.varsets_cell2struct({'Pg', 'Va'});
xx = om.varsets_x(x, vs);
t_is(length(xx), 2, 14, [t '{''Pg'', ''Va''} : length']);
t_is(xx{1}, [vv.i1.Pg:vv.iN.Pg]', 14, [t '{''Pg'', ''Va''} : 1']);
t_is(xx{2}, [vv.i1.Va:vv.iN.Va]', 14, [t '{''Pg'', ''Va''} : 2']);

vs = struct('name', 'x', 'idx', {{1,1},{2,1}});
xx = om.varsets_x(x, vs);
t_is(length(xx), 2, 14, [t '''x'', {{1,1},{2,1}} : length']);
t_is(xx{1}, [vv.i1.x(1,1):vv.iN.x(1,1)]', 14, [t '''x'', {{1,1},{2,1}} : 1']);
t_is(xx{2}, [vv.i1.x(2,1):vv.iN.x(2,1)]', 14, [t '''x'', {{1,1},{2,1}} : 2']);

vs = om.varsets_cell2struct({'x'});
xx = om.varsets_x(x, vs);
t_is(length(xx), 4, 14, [t '{''x''} : length']);
t_is(xx{1}, [vv.i1.x(1,1):vv.iN.x(1,1)]', 14, [t '{''x''} : 1']);
t_is(xx{2}, [vv.i1.x(1,2):vv.iN.x(1,2)]', 14, [t '{''x''} : 2']);
t_is(xx{3}, [vv.i1.x(2,1):vv.iN.x(2,1)]', 14, [t '{''x''} : 3']);
t_is(xx{4}, [vv.i1.x(2,2):vv.iN.x(2,2)]', 14, [t '{''x''} : 4']);

vs = om.varsets_cell2struct({'x', 'y', 'Pg'});
xx = om.varsets_x(x, vs);
t_is(length(xx), 29, 14, [t '{''x'', ''y'', ''Pg''} : length']);
t_is(xx{ 1}, [vv.i1.x(1,1):vv.iN.x(1,1)]', 14, [t '{''x'', ''y'', ''Pg''} :  1']);
t_is(xx{ 2}, [vv.i1.x(1,2):vv.iN.x(1,2)]', 14, [t '{''x'', ''y'', ''Pg''} :  2']);
t_is(xx{ 3}, [vv.i1.x(2,1):vv.iN.x(2,1)]', 14, [t '{''x'', ''y'', ''Pg''} :  3']);
t_is(xx{ 4}, [vv.i1.x(2,2):vv.iN.x(2,2)]', 14, [t '{''x'', ''y'', ''Pg''} :  4']);
t_is(xx{ 5}, [vv.i1.y(1,1,1):vv.iN.y(1,1,1)]', 14, [t '{''x'', ''y'', ''Pg''} :  5']);
t_is(xx{ 6}, [vv.i1.y(1,1,2):vv.iN.y(1,1,2)]', 14, [t '{''x'', ''y'', ''Pg''} :  6']);
t_is(xx{ 7}, [vv.i1.y(1,1,3):vv.iN.y(1,1,3)]', 14, [t '{''x'', ''y'', ''Pg''} :  7']);
t_is(xx{ 8}, [vv.i1.y(1,1,4):vv.iN.y(1,1,4)]', 14, [t '{''x'', ''y'', ''Pg''} :  8']);
t_is(xx{ 9}, [vv.i1.y(1,2,1):vv.iN.y(1,2,1)]', 14, [t '{''x'', ''y'', ''Pg''} :  9']);
t_is(xx{10}, [vv.i1.y(1,2,2):vv.iN.y(1,2,2)]', 14, [t '{''x'', ''y'', ''Pg''} : 10']);
t_is(xx{11}, [vv.i1.y(1,2,3):vv.iN.y(1,2,3)]', 14, [t '{''x'', ''y'', ''Pg''} : 11']);
t_is(xx{12}, [vv.i1.y(1,2,4):vv.iN.y(1,2,4)]', 14, [t '{''x'', ''y'', ''Pg''} : 12']);
t_is(xx{13}, [vv.i1.y(1,3,1):vv.iN.y(1,3,1)]', 14, [t '{''x'', ''y'', ''Pg''} : 13']);
t_is(xx{14}, [vv.i1.y(1,3,2):vv.iN.y(1,3,2)]', 14, [t '{''x'', ''y'', ''Pg''} : 14']);
t_is(xx{15}, [vv.i1.y(1,3,3):vv.iN.y(1,3,3)]', 14, [t '{''x'', ''y'', ''Pg''} : 15']);
t_is(xx{16}, [vv.i1.y(1,3,4):vv.iN.y(1,3,4)]', 14, [t '{''x'', ''y'', ''Pg''} : 16']);
t_is(xx{17}, [vv.i1.y(2,1,1):vv.iN.y(2,1,1)]', 14, [t '{''x'', ''y'', ''Pg''} : 17']);
t_is(xx{18}, [vv.i1.y(2,1,2):vv.iN.y(2,1,2)]', 14, [t '{''x'', ''y'', ''Pg''} : 18']);
t_is(xx{19}, [vv.i1.y(2,1,3):vv.iN.y(2,1,3)]', 14, [t '{''x'', ''y'', ''Pg''} : 19']);
t_is(xx{20}, [vv.i1.y(2,1,4):vv.iN.y(2,1,4)]', 14, [t '{''x'', ''y'', ''Pg''} : 20']);
t_is(xx{21}, [vv.i1.y(2,2,1):vv.iN.y(2,2,1)]', 14, [t '{''x'', ''y'', ''Pg''} : 21']);
t_is(xx{22}, [vv.i1.y(2,2,2):vv.iN.y(2,2,2)]', 14, [t '{''x'', ''y'', ''Pg''} : 22']);
t_is(xx{23}, [vv.i1.y(2,2,3):vv.iN.y(2,2,3)]', 14, [t '{''x'', ''y'', ''Pg''} : 23']);
t_is(xx{24}, [vv.i1.y(2,2,4):vv.iN.y(2,2,4)]', 14, [t '{''x'', ''y'', ''Pg''} : 24']);
t_is(xx{25}, [vv.i1.y(2,3,1):vv.iN.y(2,3,1)]', 14, [t '{''x'', ''y'', ''Pg''} : 25']);
t_is(xx{26}, [vv.i1.y(2,3,2):vv.iN.y(2,3,2)]', 14, [t '{''x'', ''y'', ''Pg''} : 26']);
t_is(xx{27}, [vv.i1.y(2,3,3):vv.iN.y(2,3,3)]', 14, [t '{''x'', ''y'', ''Pg''} : 27']);
t_is(xx{28}, [vv.i1.y(2,3,4):vv.iN.y(2,3,4)]', 14, [t '{''x'', ''y'', ''Pg''} : 28']);
t_is(xx{29}, [vv.i1.Pg:vv.iN.Pg]', 14, [t '{''x'', ''y'', ''Pg''} : 29']);

vs = om.varsets_cell2struct({});
xx = om.varsets_x(x, vs);
t_is(length(xx), om.var.N, 14, [t '<all> : length']);
t_is(xx, [1:om.var.N]', 14, [t '<all>']);

t = 'varsets_x(x, vs, ''vector'') : ';
vs = om.varsets_cell2struct({'Pg'});
xx = om.varsets_x(x, vs, 'vector');
t_is(length(xx), vv.N.Pg, 14, [t '{''Pg''} : length']);
t_is(xx, [vv.i1.Pg:vv.iN.Pg]', 14, [t '{''Pg''}']);

vs = om.varsets_cell2struct({'Pg', 'Va'});
xx = om.varsets_x(x, vs, 'vector');
t_is(length(xx), vv.N.Va + vv.N.Pg, 14, [t '{''Pg'', ''Va''} : length']);
t_is(xx, [vv.i1.Pg:vv.iN.Pg vv.i1.Va:vv.iN.Va]', 14, [t '{''Pg'', ''Va''}']);

vs = struct('name', 'x', 'idx', {{1,1},{2,1}});
xx = om.varsets_x(x, vs, 'vector');
t_is(length(xx), vv.N.x(1,1) + vv.N.x(2,1), 14, [t '''x'', {{1,1},{2,1}} : length']);
t_is(xx, [vv.i1.x(1,1):vv.iN.x(1,1) vv.i1.x(2,1):vv.iN.x(2,1)]', 14, [t '''x'', {{1,1},{2,1}}']);

vs = om.varsets_cell2struct({'x'});
xx = om.varsets_x(x, vs, 'vector');
t_is(length(xx), sum(vv.N.x(:)), 14, [t '{''x''} : length']);
t_is(xx, [vv.i1.x(1,1):vv.iN.x(1,1) vv.i1.x(1,2):vv.iN.x(1,2) vv.i1.x(2,1):vv.iN.x(2,1) vv.i1.x(2,2):vv.iN.x(2,2)]', 14, [t '{''x''}']);

vs = om.varsets_cell2struct({'x', 'y', 'Pg'});
xx = om.varsets_x(x, vs, 'vector');
t_is(length(xx), sum(vv.N.x(:))+sum(vv.N.y(:))+vv.N.Pg, 14, [t '{''x'', ''y'', ''Pg''} : length']);
t_is(xx, [  vv.i1.x(1,1):vv.iN.x(1,1) vv.i1.x(1,2):vv.iN.x(1,2) ...
            vv.i1.x(2,1):vv.iN.x(2,1) vv.i1.x(2,2):vv.iN.x(2,2) ...
            vv.i1.y(1,1,1):vv.iN.y(1,1,1) vv.i1.y(1,1,2):vv.iN.y(1,1,2) ...
            vv.i1.y(1,1,3):vv.iN.y(1,1,3) vv.i1.y(1,1,4):vv.iN.y(1,1,4) ...
            vv.i1.y(1,2,1):vv.iN.y(1,2,1) vv.i1.y(1,2,2):vv.iN.y(1,2,2) ...
            vv.i1.y(1,2,3):vv.iN.y(1,2,3) vv.i1.y(1,2,4):vv.iN.y(1,2,4) ...
            vv.i1.y(1,3,1):vv.iN.y(1,3,1) vv.i1.y(1,3,2):vv.iN.y(1,3,2) ...
            vv.i1.y(1,3,3):vv.iN.y(1,3,3) vv.i1.y(1,3,4):vv.iN.y(1,3,4) ...
            vv.i1.y(2,1,1):vv.iN.y(2,1,1) vv.i1.y(2,1,2):vv.iN.y(2,1,2) ...
            vv.i1.y(2,1,3):vv.iN.y(2,1,3) vv.i1.y(2,1,4):vv.iN.y(2,1,4) ...
            vv.i1.y(2,2,1):vv.iN.y(2,2,1) vv.i1.y(2,2,2):vv.iN.y(2,2,2) ...
            vv.i1.y(2,2,3):vv.iN.y(2,2,3) vv.i1.y(2,2,4):vv.iN.y(2,2,4) ...
            vv.i1.y(2,3,1):vv.iN.y(2,3,1) vv.i1.y(2,3,2):vv.iN.y(2,3,2) ...
            vv.i1.y(2,3,3):vv.iN.y(2,3,3) vv.i1.y(2,3,4):vv.iN.y(2,3,4) ...
            vv.i1.Pg:vv.iN.Pg]', 14, [t '{''x'', ''y'', ''Pg''}']);

vs = om.varsets_cell2struct({});
xx = om.varsets_x(x, vs, 'vector');
t_is(length(xx), om.var.N, 14, [t '<all> : length']);
t_is(xx, [1:om.var.N]', 14, [t '<all>']);

%%-----  add_lin_constraint  -----
t = 'add_lin_constraint';
lN = 0;
lNS = 0;
t_ok(om.getN('lin') == lN, sprintf('%s : lin.N  = %d', t, lN));
t_ok(om.get('lin', 'NS') == lNS, sprintf('%s : lin.NS = %d', t, lNS));

t = 'om.add_lin_constraint(''Pmis'', A, l, u, {''Va'', ''Pg''})';
A = sparse([1:3 1:3 1:3]', [1:3 4:6 7 7 7]', [1 1 1 -1 -1 -1 2 3 4]', 3, 7);
l = -(1:3)'; u = (1:3)';
om.add_lin_constraint('Pmis', A, l, u, {'Va', 'Pg'});
lNS = lNS + 1; lN = lN + 3;
t_ok(om.getN('lin') == lN, sprintf('%s : lin.N  = %d', t, lN));
t_ok(om.get('lin', 'NS') == lNS, sprintf('%s : lin.NS = %d', t, lNS));

t = 'om.add_lin_constraint(''Qmis'', A, l, u)';
A = sparse([1:3 1:3 1:3]', [1:3 4:6 7 7 7]', [1 1 1 -1 -1 -1 2 3 4]', 3, vN);
om.add_lin_constraint('Qmis', A, l, u);
lNS = lNS + 1; lN = lN + 3;
t_ok(om.getN('lin') == lN, sprintf('%s : lin.N  = %d', t, lN));
t_ok(om.get('lin', 'NS') == lNS, sprintf('%s : lin.NS = %d', t, lNS));

t = 'om.init_indexed_name(''lin'', ''mylin'', {2, 2})';
om.init_indexed_name('lin', 'mylin', {2, 2});
t_ok(om.getN('lin') == lN, sprintf('%s : lin.N  = %d', t, lN));
t_ok(om.get('lin', 'NS') == lNS, sprintf('%s : lin.NS = %d', t, lNS));

for i = 1:2
    for j = 1:2
        t = sprintf('om.add_lin_constraint(''mylin'', {%d,%d}, A, l, u, vs)', i,j);
        A = sparse([1:(i+j) 1:(i+j)]', [1:(i+j) 5*ones(1,i+j)]', ...
            [ones(i+j,1);-ones(i+j,1)], i+j, 3+2+(i==2 && j==1));
        l = -ones(i+j, 1); u = [];
        vs = struct('name', {'Pg', 'x'}, 'idx', {{}, {i,j}});
        om.add_lin_constraint('mylin', {i, j}, A, l, u, vs);
        lNS = lNS + 1; lN = lN + i+j;
        t_ok(om.getN('lin') == lN, sprintf('%s : lin.N  = %d', t, lN));
        t_ok(om.get('lin', 'NS') == lNS, sprintf('%s : lin.NS = %d', t, lNS));
    end
end

%%-----  add_nln_constraint (equality)  -----
t = 'add_nln_constraint (equality)';
neN = 0;
neNS = 0;
t_ok(om.getN('nle') == neN, sprintf('%s : nle.N  = %d', t, neN));
t_ok(om.get('nle', 'NS') == neNS, sprintf('%s : nle.NS = %d', t, neNS));

t = 'om.add_nln_constraint(''Pmise'', N, 1, fcn, hess, {''Pg'', ''Va''})';
N = 4;
fcn = @(x)my_fcn(x, N, 2);
hess = @(x, lam)my_hess(x, lam, 10);
om.add_nln_constraint('Pmise', N, 1, fcn, hess, {'Pg', 'Va'});
neNS = neNS + 1; neN = neN + N;
t_ok(om.getN('nle') == neN, sprintf('%s : nle.N  = %d', t, neN));
t_ok(om.get('nle', 'NS') == neNS, sprintf('%s : nle.NS = %d', t, neNS));

t = 'om.add_nln_constraint(''Qmise'', N, 1, fcn, hess)';
N = 3;
fcn = @(x)my_fcn(x, N, 2);
hess = @(x, lam)my_hess(x, lam, 10);
om.add_nln_constraint('Qmise', N, 1, fcn, hess);
neNS = neNS + 1; neN = neN + N;
t_ok(om.getN('nle') == neN, sprintf('%s : nle.N  = %d', t, neN));
t_ok(om.get('nle', 'NS') == neNS, sprintf('%s : nle.NS = %d', t, neNS));

t = 'om.add_nln_constraint({''P'',''Q'',''R''}, [3;2;1], 1, fcn, hess, {''Pg'', ''Va''})';
N = [3;2;1];
fcn = @(x)my_fcn(x, sum(N), 2);
hess = @(x, lam)my_hess(x, lam, 10);
om.add_nln_constraint({'P', 'Q', 'R'}, N, 1, fcn, hess, {'Pg', 'Va'});
neNS = neNS + length(N); neN = neN + sum(N);
t_ok(om.getN('nle') == neN, sprintf('%s : nle.N  = %d', t, neN));
t_ok(om.get('nle', 'NS') == neNS, sprintf('%s : nle.NS = %d', t, neNS));

t = 'om.init_indexed_name(''nle'', ''mynle'', {2, 2})';
om.init_indexed_name('nle', 'mynle', {2, 2});
t_ok(om.getN('nle') == neN, sprintf('%s : nle.N  = %d', t, neN));
t_ok(om.get('nle', 'NS') == neNS, sprintf('%s : nle.NS = %d', t, neNS));

for i = 1:2
    for j = 1:2
        t = sprintf('om.add_nln_constraint(''mynle'', {%d,%d}, N, 1, fcn, hess, vs)', i,j);
        N = i+j;
        fcn = @(x)my_fcn(x, N, i);
        hess = @(x, lam)my_hess(x, lam, j);
        vs = struct('name', {'Pg', 'x'}, 'idx', {{}, {i,j}});
        om.add_nln_constraint('mynle', {i, j}, N, 1, fcn, hess, vs);
        neNS = neNS + 1; neN = neN + N;
        t_ok(om.getN('nle') == neN, sprintf('%s : nle.N  = %d', t, neN));
        t_ok(om.get('nle', 'NS') == neNS, sprintf('%s : nle.NS = %d', t, neNS));
    end
end

%%-----  add_nln_constraint (inequality)  -----
t = 'add_nln_constraint (inequality)';
niN = 0;
niNS = 0;
t_ok(om.getN('nli') == niN, sprintf('%s : nli.N  = %d', t, niN));
t_ok(om.get('nli', 'NS') == niNS, sprintf('%s : nli.NS = %d', t, niNS));

t = 'om.add_nln_constraint(''Pmisi'', N, 0, fcn, hess, {''Pg'', ''Va''})';
N = 3;
fcn = @(x)my_fcn(x, N, -2);
hess = @(x, lam)my_hess(x, lam, -10);
om.add_nln_constraint('Pmisi', N, 0, fcn, hess, {'Pg', 'Va'});
niNS = niNS + 1; niN = niN + N;
t_ok(om.getN('nli') == niN, sprintf('%s : nli.N  = %d', t, niN));
t_ok(om.get('nli', 'NS') == niNS, sprintf('%s : nli.NS = %d', t, niNS));

t = 'om.add_nln_constraint(''Qmisi'', N, 0, fcn, hess)';
N = 2;
fcn = @(x)my_fcn(x, N, -2);
hess = @(x, lam)my_hess(x, lam, -10);
om.add_nln_constraint('Qmisi', N, 0, fcn, hess);
niNS = niNS + 1; niN = niN + N;
t_ok(om.getN('nli') == niN, sprintf('%s : nli.N  = %d', t, niN));
t_ok(om.get('nli', 'NS') == niNS, sprintf('%s : nli.NS = %d', t, niNS));

t = 'om.init_indexed_name(''nli'', ''mynli'', {2, 2})';
om.init_indexed_name('nli', 'mynli', {2, 2});
t_ok(om.getN('nli') == niN, sprintf('%s : nli.N  = %d', t, niN));
t_ok(om.get('nli', 'NS') == niNS, sprintf('%s : nli.NS = %d', t, niNS));

for i = 1:2
    for j = 1:2
        t = sprintf('om.add_nln_constraint(''mynli'', {%d,%d}, N, 0, fcn, hess, vs)', i,j);
        N = i+j-1;
        fcn = @(x)my_fcn(x, N, i);
        hess = @(x, lam)my_hess(x, lam, j);
        vs = struct('name', {'Pg', 'x'}, 'idx', {{}, {i,j}});
        om.add_nln_constraint('mynli', {i, j}, N, 0, fcn, hess, vs);
        niNS = niNS + 1; niN = niN + N;
        t_ok(om.getN('nli') == niN, sprintf('%s : nli.N  = %d', t, niN));
        t_ok(om.get('nli', 'NS') == niNS, sprintf('%s : nli.NS = %d', t, niNS));
    end
end

%%-----  get_idx  -----
t = 'get_idx : lin';
[vv, ll] = om.get_idx();
t_is([ll.i1.Qmis ll.iN.Qmis ll.N.Qmis], [4 6 3], 14, [t ' : Qmis']);
t_is(size(ll.i1.mylin), [2, 2], 14, [t ' : size(ll.i1.mylin)']);
t_is([ll.i1.mylin(2,1) ll.iN.mylin(2,1) ll.N.mylin(2,1)], [12 14 3], 14, [t ' : mylin(2,1)']);

t = 'get_idx(''lin'')';
ll = om.get_idx('lin');
t_is([ll.i1.Qmis ll.iN.Qmis ll.N.Qmis], [4 6 3], 14, [t ' : Qmis']);
t_is(size(ll.i1.mylin), [2, 2], 14, [t ' : size(ll.i1.mylin)']);
t_is([ll.i1.mylin(2,1) ll.iN.mylin(2,1) ll.N.mylin(2,1)], [12 14 3], 14, [t ' : mylin(2,1)']);

%%-----  params_lin_constraint  -----
t = 'params_lin_constraint';
[A, l, u] = om.params_lin_constraint();
t_ok(issparse(A), [t ' : issparse(A)']);
t_is(size(A), [lN, vN], 14, [t ' : size(A)']);
t_is(length(l), lN, 14, [t ' : length(l)']);
t_is(length(u), lN, 14, [t ' : length(u)']);
AA = sparse([1:3 1:3 1:3]', [1:3 4:6 7 7 7]', [1 1 1 -1 -1 -1 2 3 4]', 3, vN);
t_is(full(A(ll.i1.Qmis:ll.iN.Qmis, :)), full(AA), 14, [t ' : A(<Qmis>,:)']);
t_is(full(A(ll.i1.mylin(2,1):ll.iN.mylin(2,1), vv.i1.Pg:vv.iN.Pg)), eye(3,3), 14, [t ' : A(<mylin(2,1)>,<Pg>)']);
t_is(full(A(ll.i1.mylin(2,1):ll.iN.mylin(2,1), vv.i1.x(2,1):vv.iN.x(2,1))), [0 -1 0;0 -1 0;0 -1 0], 14, [t ' : A(<mylin(2,1)>,<x(2,1)>)']);

%om

%%-----  eval_nln_constraint  -----
t = 'eval_nln_constraint';
x = (1:om.var.N)';
[g, dg] = om.eval_nln_constraint(x, 1);
t_is(length(g), neN, 14, [t ' : length(g)']);
t_ok(issparse(dg), [t ' : issparse(dg)']);
t_is(size(dg), [neN, vN], 14, [t ' : size(dg)']);
t_is(g, [7 8 9 3 3 4 5 7 8 9 3 4 5 6 7 6 7 8 7 8 9 7 8 9 27]', 14, [t ' : g']);
e = [[  1 2 3 4 7 6 7;
        0 0 0 0 0 2 0;
        0 0 0 0 0 0 2;
        2 0 0 0 0 0 0 ] zeros(4, vN-7) ];
t_is(full(dg(1:4, :)), e, 14, [t ' : dg(1:4, :)   [Pmise]']);
e = [[3 2:vN]; [0 2 0 zeros(1, vN-3)]; [0 0 2 zeros(1, vN-3)]];
t_is(full(dg(5:7, :)), e, 14, [t ' : dg(5:7, :)   [Qmise]']);
e = [[  1 2 3 4 7 6 7;
        0 0 0 0 0 2 0;
        0 0 0 0 0 0 2;
        2 0 0 0 0 0 0;
        0 2 0 0 0 0 0;
        0 0 2 0 0 0 0 ] zeros(6, vN-7) ];
t_is(full(dg(8:13, :)), e, 14, [t ' : dg(8:13, :)  [mynle(1,1)]']);
e = [[0 0 0 0 6 6 7; 0 0 0 0 0 1 0] zeros(2, 10) [18 19; 0 0] zeros(2, vN-19)];
t_is(full(dg(14:15, :)), e, 14, [t ' : dg(14:15, :) [mynle(1,1)]']);
e = [[0 0 0 0 6 6 7; 0 0 0 0 0 1 0; 0 0 0 0 0 0 1] zeros(3, 12) [20 21; 0 0; 0 0] zeros(3, vN-21)];
t_is(full(dg(16:18, :)), e, 14, [t ' : dg(16:18, :) [mynle(1,2)]']);
e = [[0 0 0 0 7 6 7; 0 0 0 0 0 2 0; 0 0 0 0 0 0 2] zeros(3, 14) [22 23 24; 0 0 0; 0 0 0] zeros(3, vN-24)];
t_is(full(dg(19:21, :)), e, 14, [t ' : dg(19:21, :) [mynle(2,1)]']);
e = [[0 0 0 0 7 6 7; 0 0 0 0 0 2 0; 0 0 0 0 0 0 2; 0 0 0 0 0 0 0] zeros(4, 17) [25 26; 0 0; 0 0; 2 0] zeros(4, vN-26)];
t_is(full(dg(22:25, :)), e, 14, [t ' : dg(22:25, :) [mynle(2,2)]']);
% g
% full(dg)
% full(dg)'

[h, dh] = om.eval_nln_constraint(x, 0);
t_is(length(h), niN, 14, [t ' : length(h)']);
t_ok(issparse(dh), [t ' : issparse(dh)']);
t_is(size(dh), [niN, vN], 14, [t ' : size(dh)']);
t_is(h, [3 4 5 -1 0 6 6 7 7 8 7 8 9]', 14, [t ' : h']);
e = [[  1 2 3 4 3 6 7;
        0 0 0 0 0 -2 0;
        0 0 0 0 0 0 -2 ] zeros(3, vN-7) ];
t_is(full(dh(1:3, :)), e, 14, [t ' : dh(1:3, :)   [Pmisi]']);
e = [[-1 2:vN]; [0 -2 zeros(1, vN-2)]];
t_is(full(dh(4:5, :)), e, 14, [t ' : dh(5:7, :)   [Qmisi]']);
e = [[0 0 0 0 6 6 7] zeros(1, 10) [18 19] zeros(1, vN-19)];
t_is(full(dh(6, :)), e, 14, [t ' : dh(6, :)     [mynli(1,1)]']);
e = [[0 0 0 0 6 6 7; 0 0 0 0 0 1 0] zeros(2, 12) [20 21; 0 0] zeros(2, vN-21)];
t_is(full(dh(7:8, :)), e, 14, [t ' : dh(7:8, :)   [mynli(1,2)]']);
e = [[0 0 0 0 7 6 7; 0 0 0 0 0 2 0] zeros(2, 14) [22 23 24; 0 0 0] zeros(2, vN-24)];
t_is(full(dh(9:10, :)), e, 14, [t ' : dh(9:10, :)  [mynli(2,1)]']);
e = [[0 0 0 0 7 6 7; 0 0 0 0 0 2 0; 0 0 0 0 0 0 2] zeros(3, 17) [25 26; 0 0; 0 0] zeros(3, vN-26)];
t_is(full(dh(11:13, :)), e, 14, [t ' : dh(11:13, :) [mynli(2,2)]']);
% h
%full(dh)'

t = 'eval_nln_constraint_hess';
lam = (1:neN)'/100;
d2G = om.eval_nln_constraint_hess(x, lam, 1);
t_ok(issparse(d2G), [t ' : issparse(d2G)']);
t_is(size(d2G), [vN, vN], 14, [t ' : size(d2G)']);
% t_is(full(d2G(27:end, :)), zeros(vN-26, vN), 14, [t ' : d2G(27:end, :)']);
% t_is(full(d2G(:, 27:end)), zeros(vN, vN-26), 14, [t ' : d2G(:, 27:end)']);
e = sparse([1:3 5:7 25], [1:3 5:7 25], [33.2 24.18 26.2 56.8 62.86 60.76 27.25], vN, vN);
t_is(d2G, e, 13, [t ' : d2G']);

% d2G

lam = -(1:niN)'/100;
d2H = om.eval_nln_constraint_hess(x, lam, 0);
t_ok(issparse(d2H), [t ' : issparse(d2H)']);
t_is(size(d2H), [vN, vN], 14, [t ' : size(d2H)']);
% t_is(full(d2H(27:end, :)), zeros(vN-26, vN), 14, [t ' : d2H(27:end, :)']);
% t_is(full(d2H(:, 27:end)), zeros(vN, vN-26), 14, [t ' : d2H(:, 27:end)']);
e = sparse([1:2 5:7], [1:2 5:7], [-9.04 -8.05 20.66 18.68 5.84], vN, vN);
t_is(d2H, e, 13, [t ' : d2H']);

%d2H

%%-----  add_quad_cost  -----
t = 'add_quad_cost';
qcN = 0;
qcNS = 0;
t_ok(om.getN('qdc') == qcN, sprintf('%s : qdc.N  = %d', t, qcN));
t_ok(om.get('qdc', 'NS') == qcNS, sprintf('%s : qdc.NS = %d', t, qcNS));

t = 'om.add_quad_cost(''qc1'', <mat>Q, c, k, {''Pg'', ''Va''})';
n = nVa + nPg;
Q1 = sparse(1:n, 1:n, 1:n, n, n) + sparse(1:n, n:-1:1, 1:n, n, n);
c1 = 10*(1:n)';
k1 = n;
om.add_quad_cost('qc1', Q1, c1, k1, {'Pg', 'Va'});
qcNS = qcNS + 1; qcN = qcN + 1;
t_ok(om.getN('qdc') == qcN, sprintf('%s : qdc.N  = %d', t, qcN));
t_ok(om.get('qdc', 'NS') == qcNS, sprintf('%s : qdc.NS = %d', t, qcNS));

t = 'om.add_quad_cost(''qc2'', <mat>Q, c)';
n = om.getN('var');
Q2 = sparse(1, 1:n, 1:n, n, n) + sparse(1:n, 1, n:-1:1, n, n);
c2 = 10*(n:-1:1)';
om.add_quad_cost('qc2', Q2, c2);
qcNS = qcNS + 1; qcN = qcN + 1;
t_ok(om.getN('qdc') == qcN, sprintf('%s : qdc.N  = %d', t, qcN));
t_ok(om.get('qdc', 'NS') == qcNS, sprintf('%s : qdc.NS = %d', t, qcNS));

t = 'om.add_quad_cost(''qc3'', <vec>Q, c, k, {''Vm2'', ''Pg''})';
n = nVm2 + nPg;
Q3 = 2*(1:n)';
c3 = -1*(1:n)';
k3 = -n;
om.add_quad_cost('qc3', Q3, c3, k3, {'Vm2', 'Pg'});
qcNS = qcNS + 1; qcN = qcN + n;
t_ok(om.getN('qdc') == qcN, sprintf('%s : qdc.N  = %d', t, qcN));
t_ok(om.get('qdc', 'NS') == qcNS, sprintf('%s : qdc.NS = %d', t, qcNS));

t = 'om.add_quad_cost(''qc4'', <vec>Q, [], 0, vs)';
n = om.getN('var', 'x', {2,1}) + om.getN('var', 'y', {1,1,1});
Q4 = 1./(1:n)';
vs = struct('name', {'x', 'y'}, 'idx', {{2,1}, {1,1,1}});
om.add_quad_cost('qc4', Q4, [], 0, vs);
qcNS = qcNS + 1; qcN = qcN + n;
t_ok(om.getN('qdc') == qcN, sprintf('%s : qdc.N  = %d', t, qcN));
t_ok(om.get('qdc', 'NS') == qcNS, sprintf('%s : qdc.NS = %d', t, qcNS));

t = 'om.add_quad_cost(''qc5'', [], c, k, {''Pg'', ''Va''})';
n = nVa + nPg;
c5 = 100*(1:n)';
k5 = (1:n)';
om.add_quad_cost('qc5', [], c5, k5, {'Pg', 'Va'});
qcNS = qcNS + 1; qcN = qcN + n;
t_ok(om.getN('qdc') == qcN, sprintf('%s : qdc.N  = %d', t, qcN));
t_ok(om.get('qdc', 'NS') == qcNS, sprintf('%s : qdc.NS = %d', t, qcNS));

t = 'om.add_quad_cost(''qc6'', [], c)';
n = om.getN('var');
c6 = -(1:n)';
om.add_quad_cost('qc6', [], c6);
qcNS = qcNS + 1; qcN = qcN + n;
t_ok(om.getN('qdc') == qcN, sprintf('%s : qdc.N  = %d', t, qcN));
t_ok(om.get('qdc', 'NS') == qcNS, sprintf('%s : qdc.NS = %d', t, qcNS));

t = 'om.init_indexed_name(''qdc'', ''qc'', {2,2})';
om.init_indexed_name('qdc', 'qc', {2,2});
t_ok(om.getN('qdc') == qcN, sprintf('%s : qdc.N  = %d', t, qcN));
t_ok(om.get('qdc', 'NS') == qcNS, sprintf('%s : qdc.NS = %d', t, qcNS));

for i = 1:2
    for j = 1:2
        t = 'om.add_quad_cost(''qc'', {i, j}, cp, vs)';
        n = nPg + om.getN('var', 'x', {i,j});
        QQ = sparse(1:n, 1:n, 1:n, n, n) + sparse(1:n, n*ones(n,1), 1:n, n, n);
        cc = -2*(1:n)';
        kk = 1000;
        vs = struct('name', {'Pg', 'x'}, 'idx', {{}, {i,j}});
        om.add_quad_cost('qc', {i, j}, QQ, cc, kk, vs);
        qcNS = qcNS + 1; qcN = qcN + 1;
        t_ok(om.getN('qdc') == qcN, sprintf('%s : qdc.N  = %d', t, qcN));
        t_ok(om.get('qdc', 'NS') == qcNS, sprintf('%s : qdc.NS = %d', t, qcNS));
    end
end

%%-----  params_quad_cost  -----
t = 'om.params_quad_cost(''qc1'')';
[Q, c, k, vs] = om.params_quad_cost('qc1');
t_is(Q, Q1, 14, [t, ' : Q']);
t_is(c, c1, 14, [t, ' : c']);
t_is(k, k1, 14, [t, ' : k']);
vs1 = struct('name', {'Pg', 'Va'}, 'idx', {{}, {}});
t_ok(isequal(vs, vs1), [t, ' : vs']);

t = 'om.params_quad_cost(''qc2'')';
[Q, c, k, vs] = om.params_quad_cost('qc2');
t_is(Q, Q2, 14, [t, ' : Q']);
t_is(c, c2, 14, [t, ' : c']);
t_is(k, 0, 14, [t, ' : k']);
t_ok(isequal(vs, {}), [t, ' : vs']);

t = 'om.params_quad_cost(''qc3'')';
[Q, c, k] = om.params_quad_cost('qc3');
t_is(Q, Q3, 14, [t, ' : Q']);
t_is(c, c3, 14, [t, ' : c']);
t_is(k, k3, 14, [t, ' : k']);

t = 'om.params_quad_cost(''qc4'')';
[Q, c] = om.params_quad_cost('qc4');
t_is(Q, Q4, 14, [t, ' : Q']);
t_ok(isempty(c), [t, ' : c']);
% t_is(k, 0, 14, [t, ' : k']);

t = 'om.params_quad_cost(''qc5'')';
[Q, c, k] = om.params_quad_cost('qc5');
t_ok(isempty(Q), [t, ' : Q']);
t_is(c, c5, 14, [t, ' : c']);
t_is(k, k5, 14, [t, ' : k']);

t = 'om.params_quad_cost(''qc6'')';
[Q, c, k] = om.params_quad_cost('qc6');
t_ok(isempty(Q), [t, ' : Q']);
t_is(c, c6, 14, [t, ' : c']);
t_is(k, 0, 14, [t, ' : k']);

for i = 1:2
    for j = 1:2
        t = sprintf('om.params_quad_cost(''qc'', {%d, %d})', i, j);
        n = nPg + om.getN('var', 'x', {i,j});
        QQ = sparse(1:n, 1:n, 1:n, n, n) + sparse(1:n, n*ones(n,1), 1:n, n, n);
        cc = -2*(1:n)';
        kk = 1000;
        [Q, c, k] = om.params_quad_cost('qc', {i,j});
        t_is(Q, QQ, 14, [t, ' : Q']);
        t_is(c, cc, 14, [t, ' : c']);
        t_is(k, kk, 14, [t, ' : k']);
    end
end

t = 'om.params_quad_cost()';
[Q, c, k] = om.params_quad_cost();
% [ii, jj, ss] = find(Q)
ii = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 1 2 7 1 3 6 1 4 5 1 4 5 1 3 6 1 2 7 1 1 1 1 1 1 13 1 14 1 15 1 16 1 17 1 18 1 5 6 7 18 19 1 20 1 5 6 7 20 21 1 22 1 23 1 5 6 7 22 23 24 1 25 1 5 6 7 25 26 1 27 1 28 1 29 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]';
jj = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 9 10 11 12 13 13 14 14 15 15 16 16 17 17 18 18 19 19 19 19 19 19 20 20 21 21 21 21 21 21 22 22 23 23 24 24 24 24 24 24 24 25 25 26 26 26 26 26 26 27 27 28 28 29 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170]';
ss = [179 169 168 167 166 165 164 163 162 161 160 159 158 157 156 155 154 153 152 151 150 149 148 147 146 145 144 143 142 141 140 139 138 137 136 135 134 133 132 131 130 129 128 127 126 125 124 123 122 121 120 119 118 117 116 115 114 113 112 111 110 109 108 107 106 105 104 103 102 101 100 99 98 97 96 95 94 93 92 91 90 89 88 87 86 85 84 83 82 81 80 79 78 77 76 75 74 73 72 71 70 69 68 67 66 65 64 63 62 61 60 59 58 57 56 55 54 53 52 51 50 49 48 47 46 45 44 43 42 41 40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 2 5 3 3 6 2 4 7 1 5 7 17 6 6 24 7 5 31 8 9 10 11 12 13 2 14 4 15 6 16 8 17 10 18 4 19 1 2 3 4 10 20 4 21 1 2 3 4 10 22 5 23 5.5 24 1 2 3 4 5 12.3333333333333333 25 4 26 1 2 3 4 10 27 0.25 28 0.2 29 0.1666666666666667 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170];
QQ = sparse(ii, jj, ss, om.var.N, om.var.N);
cc = [ 2139 2238 2337 2436 1751 1841 1931 1622 1611 1600 1589 1578 1566 1554 1542 1530 1518 1504 1491 1482 1469 1460 1447 1434 1427 1414 1413 1402 1391 1380 1369 1358 1347 1336 1325 1314 1303 1292 1281 1270 1259 1248 1237 1226 1215 1204 1193 1182 1171 1160 1149 1138 1127 1116 1105 1094 1083 1072 1061 1050 1039 1028 1017 1006 995 984 973 962 951 940 929 918 907 896 885 874 863 852 841 830 819 808 797 786 775 764 753 742 731 720 709 698 687 676 665 654 643 632 621 610 599 588 577 566 555 544 533 522 511 500 489 478 467 456 445 434 423 412 401 390 379 368 357 346 335 324 313 302 291 280 269 258 247 236 225 214 203 192 181 170 159 148 137 126 115 104 93 82 71 60 49 38 27 16 5 -6 -17 -28 -39 -50 -61 -72 -83 -94 -105 -116 -127 -138 -149 -160]';
t_is(Q, QQ, 14, [t, ' : Q']);
t_is(c, cc, 14, [t, ' : c']);
t_is(k, 4027, 14, [t, ' : k']);

%%-----  eval_quad_cost  -----
t = 'om.eval_quad_cost(x, ''qc1'')';
x = (1:om.var.N)';
[Q, c, k, vs] = om.params_quad_cost('qc1');
xx = om.varsets_x(x, vs, 'vector');
ef = 1/2 * xx'*Q*xx + c'*xx + k;
edf = Q*xx + c;
f = om.eval_quad_cost(x, 'qc1');
t_is(f, ef, 14, [t, ' : f']);
[f, df] = om.eval_quad_cost(x, 'qc1');
t_is(f, ef, 14, [t, ' : f']);
t_is(df, edf, 14, [t, ' : df']);
[f, df, d2f] = om.eval_quad_cost(x, 'qc1');
t_is(f, ef, 14, [t, ' : f']);
t_is(df, edf, 14, [t, ' : df']);
t_is(d2f, Q, 14, [t, ' : d2f']);

t = 'om.eval_quad_cost(x, ''qc2'')';
[Q, c, k, vs] = om.params_quad_cost('qc2');
xx = om.varsets_x(x, vs, 'vector');
ef = 1/2 * xx'*Q*xx + c'*xx + k;
edf = Q*xx + c;
[f, df, d2f] = om.eval_quad_cost(x, 'qc2');
t_is(f, ef, 14, [t, ' : f']);
t_is(df, edf, 14, [t, ' : df']);
t_is(d2f, Q, 14, [t, ' : d2f']);

t = 'om.eval_quad_cost(x, ''qc3'')';
[Q, c, k, vs] = om.params_quad_cost('qc3');
xx = om.varsets_x(x, vs, 'vector');
ef = 1/2 * Q.*xx.^2 + c.*xx + k;
edf = Q.*xx + c;
[f, df, d2f] = om.eval_quad_cost(x, 'qc3');
t_is(f, ef, 14, [t, ' : f']);
t_is(df, edf, 14, [t, ' : df']);
t_is(d2f, Q, 14, [t, ' : d2f']);

t = 'om.eval_quad_cost(x, ''qc4'')';
[Q, c, k, vs] = om.params_quad_cost('qc4');
xx = om.varsets_x(x, vs, 'vector');
ef = 1/2 * Q.*xx.^2 + k;
edf = Q.*xx;
[f, df, d2f] = om.eval_quad_cost(x, 'qc4');
t_is(f, ef, 14, [t, ' : f']);
t_is(df, edf, 14, [t, ' : df']);
t_is(d2f, Q, 14, [t, ' : d2f']);

t = 'om.eval_quad_cost(x, ''qc5'')';
[Q, c, k, vs] = om.params_quad_cost('qc5');
xx = om.varsets_x(x, vs, 'vector');
ef = c.*xx + k;
edf = c;
[f, df, d2f] = om.eval_quad_cost(x, 'qc5');
t_is(f, ef, 14, [t, ' : f']);
t_is(df, edf, 14, [t, ' : df']);
t_is(d2f, sparse(length(xx), 1), 14, [t, ' : d2f']);

t = 'om.eval_quad_cost(x, ''qc6'')';
[Q, c, k, vs] = om.params_quad_cost('qc6');
xx = x;
ef = c.*xx + k;
edf = c;
[f, df, d2f] = om.eval_quad_cost(x, 'qc6');
t_is(f, ef, 14, [t, ' : f']);
t_is(df, edf, 14, [t, ' : df']);
t_is(d2f, sparse(length(x), 1), 14, [t, ' : d2f']);

for i = 1:2
    for j = 1:2
        t = sprintf('om.eval_quad_cost(x, ''qc'', {%d, %d})', i, j);
        [Q, c, k, vs] = om.params_quad_cost('qc', {i,j});
        xx = om.varsets_x(x, vs, 'vector');
        ef = 1/2 * xx'*Q*xx + c'*xx + k;
        edf = Q*xx + c;
        [f, df, d2f] = om.eval_quad_cost(x, 'qc', {i,j});
        t_is(f, ef, 14, [t, ' : f']);
        t_is(df, edf, 14, [t, ' : df']);
        t_is(d2f, Q, 14, [t, ' : d2f']);
    end
end

t = 'om.eval_quad_cost(x)';
[Q, c, k] = om.params_quad_cost();
xx = x;
ef = 1/2 * xx'*Q*xx + c'*xx + k;
edf = Q*xx + c;
[f, df, d2f] = om.eval_quad_cost(x);
t_is(f, ef, 14, [t, ' : f']);
t_is(df, edf, 14, [t, ' : df']);
t_is(d2f, Q, 14, [t, ' : d2f']);


%%-----  add_legacy_cost  -----
t = 'add_legacy_cost';
cN = 0;
cNS = 0;
t_ok(om.getN('cost') == cN, sprintf('%s : cost.N  = %d', t, cN));
t_ok(om.get('cost', 'NS') == cNS, sprintf('%s : cost.NS = %d', t, cNS));

t = 'om.add_legacy_cost(''ucost'', cp, {''Va'', ''Pg''})';
cp = struct('N', sparse([1:2 1:2 1:2]', [1:4 5 7]', [1 1 -1 -1 2 2]', 2,7), ...
            'Cw', [2;3]);
om.add_legacy_cost('ucost', cp, {'Va', 'Pg'});
cNS = cNS + 1; cN = cN + 2;
t_ok(om.getN('cost') == cN, sprintf('%s : cost.N  = %d', t, cN));
t_ok(om.get('cost', 'NS') == cNS, sprintf('%s : cost.NS = %d', t, cNS));

t = 'om.add_legacy_cost(''vcost'', cp)';
cp = struct('N', sparse([1:2 1:2 1:2]', [1:4 5 7]', [1 1 -1 -1 2 2]', 2, vN), ...
            'Cw', [2;3]);
om.add_legacy_cost('vcost', cp);
cNS = cNS + 1; cN = cN + 2;
t_ok(om.getN('cost') == cN, sprintf('%s : cost.N  = %d', t, cN));
t_ok(om.get('cost', 'NS') == cNS, sprintf('%s : cost.NS = %d', t, cNS));

t = 'om.init_indexed_name(''cost'', ''wc'', {2,2})';
om.init_indexed_name('cost', 'wc', {2,2});
t_ok(om.getN('cost') == cN, sprintf('%s : cost.N  = %d', t, cN));
t_ok(om.get('cost', 'NS') == cNS, sprintf('%s : cost.NS = %d', t, cNS));

for i = 1:2
    for j = 1:2
        t = 'om.add_legacy_cost(''wc'', {i, j}, cp, vs)';
        cp.N = sparse([1:(i+j) 1:(i+j)]', [1:(i+j) 5*ones(1,i+j)]', ...
            [ones(i+j,1);-ones(i+j,1)], i+j, 3+2+(i==2 && j==1));
        cp.Cw = (i+j:-1:1)';
        if i == 2
            cp.H = sparse((1:i+j)', (1:i+j)', (1:i+j)', i+j, i+j);
        end
        vs = struct('name', {'Pg', 'x'}, 'idx', {{}, {i,j}});
        om.add_legacy_cost('wc', {i, j}, cp, vs);
        cNS = cNS + 1; cN = cN + i+j;
        t_ok(om.getN('cost') == cN, sprintf('%s : cost.N  = %d', t, cN));
        t_ok(om.get('cost', 'NS') == cNS, sprintf('%s : cost.NS = %d', t, cNS));
    end
end

t = 'params_legacy_cost';
cp = om.params_legacy_cost();
t_ok(isfield(cp, 'N'), t);

%%-----  get_idx  -----
t = 'get_idx : cost';
[vv, ll, nne, nni, cc] = om.get_idx();
t_is([cc.i1.vcost cc.iN.vcost cc.N.vcost], [3 4 2], 14, [t ' : vcost']);
t_is(size(cc.i1.wc), [2, 2], 14, [t ' : size(cc.i1.wc)']);
t_is([cc.i1.wc(2,1) cc.iN.wc(2,1) cc.N.wc(2,1)], [10 12 3], 14, [t ' : wc(2,1)']);

t = 'get_idx(''var'', ''cost'', ''lin'')';
[vv, cc, ll] = om.get_idx('var', 'cost', 'lin');
t_is([vv.i1.Pg vv.iN.Pg vv.N.Pg], [5 7 3], 14, [t ' : Pg']);
t_is(size(vv.i1.x), [2, 2], 14, [t ' : size(vv.i1.x)']);
t_is([vv.i1.x(2,1) vv.iN.x(2,1) vv.N.x(2,1)], [22 24 3], 14, [t ' : x(2,1)']);
t_is(size(vv.i1.y), [2, 3, 4], 14, [t ' : size(vv.i1.y)']);
t_is([vv.i1.y(2,2,4) vv.iN.y(2,2,4) vv.N.y(2,2,4)], [133 140 8], 14, [t ' : y(2,2,4)']);
t_is([cc.i1.vcost cc.iN.vcost cc.N.vcost], [3 4 2], 14, [t ' : vcost']);
t_is(size(cc.i1.wc), [2, 2], 14, [t ' : size(cc.i1.wc)']);
t_is([cc.i1.wc(2,1) cc.iN.wc(2,1) cc.N.wc(2,1)], [10 12 3], 14, [t ' : wc(2,1)']);
t_is([ll.i1.Qmis ll.iN.Qmis ll.N.Qmis], [4 6 3], 14, [t ' : Qmis']);
t_is(size(ll.i1.mylin), [2, 2], 14, [t ' : size(ll.i1.mylin)']);
t_is([ll.i1.mylin(2,1) ll.iN.mylin(2,1) ll.N.mylin(2,1)], [12 14 3], 14, [t ' : mylin(2,1)']);

%%-----  params_legacy_cost  -----
t = 'om.params_legacy_cost(''ucost'')';
cp = om.params_legacy_cost('ucost');
n = nVa + nPg;
N = sparse([1:2 1:2 1:2]', [1:4 5 7]', [1 1 -1 -1 2 2]', 2, n);
t_is(full(cp.N), full(N), 14, [t, ' : N']);
t_is(cp.Cw, [2;3], 14, [t, ' : Cw']);
t_is(full(cp.H), zeros(2,2), 14, [t, ' : H']);
t_is(cp.dd, ones(2,1),  14, [t, ' : dd']);
t_is(cp.rh, zeros(2,1), 14, [t, ' : rh']);
t_is(cp.kk, zeros(2,1), 14, [t, ' : kk']);
t_is(cp.mm, ones(2,1),  14, [t, ' : mm']);

t = 'om.params_legacy_cost(''vcost'')';
cp = om.params_legacy_cost('vcost');
N = sparse([1:2 1:2 1:2]', [1:4 5 7]', [1 1 -1 -1 2 2]', 2, vN);
t_is(full(cp.N), full(N), 14, [t, ' : N']);
t_is(cp.Cw, [2;3], 14, [t, ' : Cw']);
t_is(full(cp.H), zeros(2,2), 14, [t, ' : H']);
t_is(cp.dd, ones(2,1),  14, [t, ' : dd']);
t_is(cp.rh, zeros(2,1), 14, [t, ' : rh']);
t_is(cp.kk, zeros(2,1), 14, [t, ' : kk']);
t_is(cp.mm, ones(2,1),  14, [t, ' : mm']);

t = 'om.params_legacy_cost(''wc'') : error';
try
    cp = om.params_legacy_cost('wc')
    t_ok(0, t);
catch
    t_ok(strfind(lasterr, '@opt_model/params_legacy_cost: legacy cost set ''wc'' requires an IDX arg'), t);
end

t = 'om.params_legacy_cost(''wc'', {1,2})';
cp = om.params_legacy_cost('wc', {1,2});
n = nPg + om.getN('var', 'x', {1,2});
N = sparse([1:3 1:3]', [1:3 nPg+2*ones(1,3)]', [ones(3,1);-ones(3,1)], 3, n);
t_is(full(cp.N), full(N), 14, [t, ' : N']);
t_is(cp.Cw, [3;2;1], 14, [t, ' : Cw']);
t_is(full(cp.H), zeros(3,3), 14, [t, ' : H']);
t_is(cp.dd, ones(3,1),  14, [t, ' : dd']);
t_is(cp.rh, zeros(3,1), 14, [t, ' : rh']);
t_is(cp.kk, zeros(3,1), 14, [t, ' : kk']);
t_is(cp.mm, ones(3,1),  14, [t, ' : mm']);

t = 'om.params_legacy_cost(''wc'', {2,1})';
cp = om.params_legacy_cost('wc', {2,1});
n = nPg + om.getN('var', 'x', {2,1});
N = sparse([1:3 1:3]', [1:3 nPg+2*ones(1,3)]', [ones(3,1);-ones(3,1)], 3, n);
t_is(full(cp.N), full(N), 14, [t, ' : N']);
t_is(cp.Cw, [3;2;1], 14, [t, ' : Cw']);
H = sparse(1:3, 1:3, 1:3, 3, 3);
t_is(full(cp.H), full(H), 14, [t, ' : H']);
t_is(cp.dd, ones(3,1),  14, [t, ' : dd']);
t_is(cp.rh, zeros(3,1), 14, [t, ' : rh']);
t_is(cp.kk, zeros(3,1), 14, [t, ' : kk']);
t_is(cp.mm, ones(3,1),  14, [t, ' : mm']);

t = 'om.params_legacy_cost()';
cp = om.params_legacy_cost();
t_ok(issparse(cp.N), [t ' : issparse(cp.N)']);
t_is(size(cp.N), [cN, vN], 14, [t ' : size(cp.N)']);
t_is(size(cp.H), [cN, cN], 14, [t ' : size(cp.H)']);
t_is(length(cp.Cw), cN, 14, [t ' : length(cp.Cw)']);
t_is(length(cp.dd), cN, 14, [t ' : length(cp.dd)']);
t_is(length(cp.rh), cN, 14, [t ' : length(cp.rh)']);
t_is(length(cp.kk), cN, 14, [t ' : length(cp.kk)']);
t_is(length(cp.mm), cN, 14, [t ' : length(cp.mm)']);
N = sparse([1:2 1:2]', [1:4]', [1 1 -1 -1]', 2, 4);
Cw = [2;3];
H = zeros(2,2);
t_is(full(cp.N(cc.i1.vcost:cc.iN.vcost, vv.i1.Va:vv.iN.Va)), full(N), 14, [t ' : N(<vcost>,<Va>)']);
N = sparse([1:2]', [1 3]', [2 2]', 2, 3);
t_is(full(cp.N(cc.i1.vcost:cc.iN.vcost, vv.i1.Pg:vv.iN.Pg)), full(N), 14, [t ' : N(<vcost>,<Pg>)']);
t_is(full(cp.Cw(cc.i1.vcost:cc.iN.vcost)), Cw, 14, [t ' : Cw(<vcost>)']);
t_is(full(cp.H(cc.i1.vcost:cc.iN.vcost, cc.i1.vcost:cc.iN.vcost)), full(H), 14, [t ' : H(<vcost>,<vcost>)']);
N = sparse([1:3]', [1:3]', [1 1 1]', 3, 3);
t_is(full(cp.N(cc.i1.wc(1,2):cc.iN.wc(1,2), vv.i1.Pg:vv.iN.Pg)), full(N), 14, [t ' : N(<wc(1,2)>,<Pg>)']);
N = sparse([1:3]', [2 2 2]', [-1 -1 -1]', 3, 2);
Cw = [3;2;1];
H = zeros(3,3);
t_is(full(cp.N(cc.i1.wc(1,2):cc.iN.wc(1,2), vv.i1.x(1,2):vv.iN.x(1,2))), full(N), 14, [t ' : N(<wc(1,2)>,<x(1,2)>)']);
t_is(full(cp.Cw(cc.i1.wc(1,2):cc.iN.wc(1,2))), Cw, 14, [t ' : Cw(<wc(1,2)>)']);
t_is(full(cp.H(cc.i1.wc(1,2):cc.iN.wc(1,2), cc.i1.wc(1,2):cc.iN.wc(1,2))), full(H), 14, [t ' : H(<wc(1,2)>,<wc(1,2)>)']);

N = sparse([1:3]', [1:3]', [1 1 1]', 3, 3);
t_is(full(cp.N(cc.i1.wc(2,1):cc.iN.wc(2,1), vv.i1.Pg:vv.iN.Pg)), full(N), 14, [t ' : N(<wc(2,1)>,<Pg>)']);
N = sparse([1:3]', [2 2 2]', [-1 -1 -1]', 3, 3);
Cw = [3;2;1];
H = full(sparse((1:3)', (1:3)', (1:3)', 3, 3));
t_is(full(cp.N(cc.i1.wc(2,1):cc.iN.wc(2,1), vv.i1.x(2,1):vv.iN.x(2,1))), full(N), 14, [t ' : N(<wc(2,1)>,<x(2,1)>)']);
t_is(full(cp.Cw(cc.i1.wc(2,1):cc.iN.wc(2,1))), Cw, 14, [t ' : Cw(<wc(2,1)>)']);
t_is(full(cp.H(cc.i1.wc(2,1):cc.iN.wc(2,1), cc.i1.wc(2,1):cc.iN.wc(2,1))), full(H), 14, [t ' : H(<wc(2,1)>,<wc(2,1)>)']);

t_is(cp.dd, ones(cN,1),  14, [t, ' : dd']);
t_is(cp.rh, zeros(cN,1), 14, [t, ' : rh']);
t_is(cp.kk, zeros(cN,1), 14, [t, ' : kk']);
t_is(cp.mm, ones(cN,1),  14, [t, ' : mm']);

%%-----  eval_legacy_cost  -----
t = 'om.eval_legacy_cost(x)';
x = [1:7 rand(1,10) 8:(vN-10)]';
f = om.eval_legacy_cost(x);
t_is(f, 343, 14, t);

t = 'om.eval_legacy_cost(''ucost'')';
f = om.eval_legacy_cost(x, 'ucost');
t_is(f, 52, 14, t);

t = 'om.eval_legacy_cost(''wc'', {2,1})';
f = om.eval_legacy_cost(x, 'wc', {2,1});
t_is(f, 91, 14, t);

t = 'om.eval_legacy_cost(''wc'')';
f = om.eval_legacy_cost(x, 'wc');
t_is(f, 239, 14, t);

%%-----  add_nln_cost  -----
t = 'add_nln_cost';
nlcN = 0;
nlcNS = 0;
t_ok(om.getN('nlc') == nlcN, sprintf('%s : nlc.N  = %d', t, nlcN));
t_ok(om.get('nlc', 'NS') == nlcNS, sprintf('%s : nlc.NS = %d', t, nlcNS));

t = 'om.add_nln_cost(''ucost'', 1, fcn, {''Va'', ''Pg''})';
cp = struct('N', sparse([1:2 1:2 1:2]', [1:4 5 7]', [1 1 -1 -1 2 2]', 2,7), ...
            'Cw', [2;3]);
fcn = @(x)my_legacy_cost_fcn(x, cp, om, {'Va', 'Pg'});
om.add_nln_cost('ucost', 1, fcn, {'Va', 'Pg'});
nlcNS = nlcNS + 1; nlcN = nlcN + 1;
t_ok(om.getN('nlc') == nlcN, sprintf('%s : nlc.N  = %d', t, nlcN));
t_ok(om.get('nlc', 'NS') == nlcNS, sprintf('%s : nlc.NS = %d', t, nlcNS));

t = 'om.add_nln_cost(''vcost'', cp)';
cp = struct('N', sparse([1:2 1:2 1:2]', [1:4 5 7]', [1 1 -1 -1 2 2]', 2, vN), ...
            'Cw', [2;3]);
fcn = @(x)my_legacy_cost_fcn(x, cp, om);
om.add_nln_cost('vcost', 1, fcn);
nlcNS = nlcNS + 1; nlcN = nlcN + 1;
t_ok(om.getN('nlc') == nlcN, sprintf('%s : nlc.N  = %d', t, nlcN));
t_ok(om.get('nlc', 'NS') == nlcNS, sprintf('%s : nlc.NS = %d', t, nlcNS));

t = 'om.init_indexed_name(''nlc'', ''wc'', {2,2})';
om.init_indexed_name('nlc', 'wc', {2,2});
t_ok(om.getN('nlc') == nlcN, sprintf('%s : nlc.N  = %d', t, nlcN));
t_ok(om.get('nlc', 'NS') == nlcNS, sprintf('%s : nlc.NS = %d', t, nlcNS));

for i = 1:2
    for j = 1:2
        t = 'om.add_nln_cost(''wc'', {i, j}, cp, vs)';
        cp.N = sparse([1:(i+j) 1:(i+j)]', [1:(i+j) 5*ones(1,i+j)]', ...
            [ones(i+j,1);-ones(i+j,1)], i+j, 3+2+(i==2 && j==1));
        cp.Cw = (i+j:-1:1)';
        if i == 2
            cp.H = sparse((1:i+j)', (1:i+j)', (1:i+j)', i+j, i+j);
        end
        vs = struct('name', {'Pg', 'x'}, 'idx', {{}, {i,j}});
        fcn = @(x)my_legacy_cost_fcn(x, cp, om, vs);
        om.add_nln_cost('wc', {i, j}, 1, fcn, vs);
        nlcNS = nlcNS + 1; nlcN = nlcN + 1;
        t_ok(om.getN('nlc') == nlcN, sprintf('%s : nlc.N  = %d', t, nlcN));
        t_ok(om.get('nlc', 'NS') == nlcNS, sprintf('%s : nlc.NS = %d', t, nlcNS));
    end
end

%%-----  params_nln_cost  -----
t = 'om.params_nln_cost(''ucost'')';
[N, fcn] = om.params_nln_cost('ucost');
t_is(N, 1, 14, [t, ' : N']);
t_ok(isa(fcn, 'function_handle'), [t, ' : fcn']);

t = 'om.params_nln_cost(''vcost'')';
[N, fcn] = om.params_nln_cost('vcost');
t_is(N, 1, 14, [t, ' : N']);
t_ok(isa(fcn, 'function_handle'), [t, ' : fcn']);

t = 'om.params_nln_cost(''wc'') : error';
try
    [N, fcn] = om.params_nln_cost('wc')
    t_ok(0, t);
catch
    t_ok(strfind(lasterr, '@opt_model/params_nln_cost: general nonlinear cost set ''wc'' requires an IDX arg'), t);
end

t = 'om.params_nln_cost(''wc'', {1,2})';
[N, fcn] = om.params_nln_cost('wc', {1,2});
t_is(N, 1, 14, [t, ' : N']);
t_ok(isa(fcn, 'function_handle'), [t, ' : fcn']);

t = 'om.params_nln_cost(''wc'', {2,1})';
[N, fcn] = om.params_nln_cost('wc', {2,1});
t_is(N, 1, 14, [t, ' : N']);
t_ok(isa(fcn, 'function_handle'), [t, ' : fcn']);

%%-----  eval_nln_cost  -----
t = 'om.eval_nln_cost(x) : ';
x = [1:7 rand(1,10) 8:(vN-10)]';
f = om.eval_nln_cost(x);
ef = 343;
ii = [1 2 3 4 5 6 7 19 21 23 25 26];
jj = [1 1 1 1 1 1 1 1 1 1 1 1];
ss = [4 6 -4 -6 1 -26 -29 -3 -6 34 -3 52];
edf = full(sparse(ii,jj,ss,vN,1));
ii = [5 23 26 6 23 26 7 23 26 5 6 7 23 25 26 5 6 7 25 26];
jj = [5 5 5 6 6 6 7 7 7 23 23 23 23 25 25 26 26 26 26 26];
ss = [2 -1 -1 4 -2 -2 6 -3 -3 -1 -2 -3 6 4 -4 -1 -2 -3 -4 10];
ed2f = full(sparse(ii,jj,ss,vN,vN));
t_is(f, ef, 14, [t 'f']);
[f, df] = om.eval_nln_cost(x);
t_is(f, ef, 14, [t 'f']);
t_is(df, edf, 14, [t 'df']);
[f, df, d2f] = om.eval_nln_cost(x);
t_is(f, ef, 14, [t 'f']);
t_is(df, edf, 14, [t 'df']);
t_is(d2f, ed2f, 14, [t 'd2f']);

t = 'om.eval_nln_cost(''ucost'') : ';
f = om.eval_nln_cost(x, 'ucost');
ef = 52;
edf = [2; 3; -2; -3; 4; 0; 6];
t_is(f, ef, 14, [t 'f']);
[f, df] = om.eval_nln_cost(x, 'ucost');
t_is(f, ef, 14, [t 'f']);
t_is(df, edf, 14, [t 'df']);
[f, df, d2f] = om.eval_nln_cost(x, 'ucost');
t_is(f, ef, 14, [t 'f']);
t_is(df, edf, 14, [t 'df']);
t_is(full(d2f), zeros(7,7), 14, [t 'd2f']);

t = 'om.eval_nln_cost(''wc'', {2,1}) : ';
f = om.eval_nln_cost(x, 'wc', {2,1});
ef = 91;
edf = [-5; -12; -17; 0; 34; 0];
ii = [1 5 2 5 3 5 1 2 3 5];
jj = [1 1 2 2 3 3 5 5 5 5];
ss = [1 -1 2 -2 3 -3 -1 -2 -3 6];
ed2f = full(sparse(ii,jj,ss,6,6));
t_is(f, ef, 14, [t 'f']);
[f, df] = om.eval_nln_cost(x, 'wc', {2,1});
t_is(f, ef, 14, [t 'f']);
t_is(df, edf, 14, [t 'df']);
[f, df, d2f] = om.eval_nln_cost(x, 'wc', {2,1});
t_is(f, ef, 14, [t 'f']);
t_is(df, edf, 14, [t 'df']);
t_is(d2f, ed2f, 14, [t 'd2f']);

t = 'om.eval_nln_cost(''wc'') : ';
f = om.eval_nln_cost(x, 'wc');
t_is(f, 239, 14, [t 'f']);

% om
% om = struct(om);

t_end

function [g, dg] = my_fcn(x, p1, p2)
if iscell(x)
    xx = [];
    for k = 1:length(x)
        xx = [xx; x{k}];
    end
else
    xx = x;
end
M = p1;
N = length(xx);
if M > N
    error('M <= length(x)');
end
g = xx(1:M) + p2;
dg = sparse(1:M, 1:M, p2, M, N) + sparse(1, 1:N, xx, M, N);

function d2G = my_hess(x, lam, p3)
if iscell(x)
    xx = [];
    for k = 1:length(x)
        xx = [xx; x{k}];
    end
else
    xx = x;
end
N = length(xx);
M = length(lam);
MM = min(M, N);
d2G = sparse(1:MM, 1:MM, xx(1:MM) + lam(1:MM) + p3, N, N);
%full(d2G(1:MM,1:MM))

function [varargout] = my_legacy_cost_fcn(x, cp, om, vs)
[nw, nx] = size(cp.N);
if ~isfield(cp, 'H') || isempty(cp.H)
    cp.H = sparse(nw, nw);
end
if ~isfield(cp, 'dd') || isempty(cp.dd)
    cp.dd = ones(nw, 1);
end
if ~isfield(cp, 'rh') || isempty(cp.rh)
    cp.rh = zeros(nw, 1);
end
if ~isfield(cp, 'kk') || isempty(cp.kk)
    cp.kk = zeros(nw, 1);
end
if ~isfield(cp, 'mm') || isempty(cp.mm)
    cp.mm = ones(nw, 1);
end
if nargin < 4
    vs = {};
end
if iscell(x)
    xx = vertcat(x{:});
else
    xx = x;
end
vs = om.varsets_cell2struct(vs);
[varargout{1:nargout}] = opf_legacy_user_cost_fcn(xx, cp);
