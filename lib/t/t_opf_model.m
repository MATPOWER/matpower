function t_opf_model(quiet)
%T_OPF_MODEL Tests for OPF_MODEL.

%   MATPOWER
%   Copyright (c) 2012-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

num_tests = 305;

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

%%-----  add_vars  -----
t = 'add_vars';
vN = 0;
vNS = 0;
t_ok(getN(om, 'var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(get(om, 'var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'add_vars(om, ''Va'', 4)';
om = add_vars(om, 'Va', 4);
vNS = vNS + 1; vN = vN + 4;
t_ok(getN(om, 'var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(get(om, 'var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'add_vars(om, ''Pg'', 3, Pg0, Pgmin, Pgmax)';
om = add_vars(om, 'Pg', 3, [2;4;6], [1;2;3], [10;20;30]);
vNS = vNS + 1; vN = vN + 3;
t_ok(getN(om, 'var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(get(om, 'var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'add_vars(om, ''Vm1'', 5, V0, Vmin, Vmax, ''I'')';
V0 = [1;1;1;1;1];
Vmin = zeros(5, 1);
Vmax = 1 + 0.01*(1:5)';
vt = 'I';
om = add_vars(om, 'Vm1', 5, V0, Vmin, Vmax, vt);
vNS = vNS + 1; vN = vN + 5;
t_ok(getN(om, 'var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(get(om, 'var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'add_vars(om, ''Vm2'', 5, V0, Vmin, Vmax, ''CIBIC'')';
vt = 'CIBIC';
om = add_vars(om, 'Vm2', 5, V0, Vmin, Vmax, vt);
vNS = vNS + 1; vN = vN + 5;
t_ok(getN(om, 'var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(get(om, 'var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'add_vars(om, ''x'', dims)';
om = add_vars(om, 'x', {2,2});
t_ok(getN(om, 'var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(get(om, 'var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'add_vars(om, ''x'', {1,1}, 2)';
om = add_vars(om, 'x', {1,1}, 2);
vNS = vNS + 1; vN = vN + 2;
t_ok(getN(om, 'var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(get(om, 'var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'add_vars(om, ''x'', {1,2}, 2, x0(1,2))';
om = add_vars(om, 'x', {1,2}, 2, [-1;-2]);
vNS = vNS + 1; vN = vN + 2;
t_ok(getN(om, 'var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(get(om, 'var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'add_vars(om, ''x'', {2,1}, 3)';
om = add_vars(om, 'x', {2,1}, 3);
vNS = vNS + 1; vN = vN + 3;
t_ok(getN(om, 'var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(get(om, 'var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'add_vars(om, ''x'', {2,2}, 2, x0(2,2), xmin(2,2), xmax(2,2))';
om = add_vars(om, 'x', {2,2}, 2, [1;0],[0;-1],[2;1]);
vNS = vNS + 1; vN = vN + 2;
t_ok(getN(om, 'var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(get(om, 'var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'add_vars(om, ''y'', {2,3,4})';
om = add_vars(om, 'y', {2,3,4});
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
            t = sprintf('add_vars(om, ''y'', {%d,%d,%d}, y0, ymin, ymax, vt)', i,j,k);
            om = add_vars(om, 'y', {i,j,k}, n, 10*(n:-1:1)', -1*(n:-1:1)', 100+(n:-1:1)', vt);
            vNS = vNS + 1; vN = vN + n;
            t_ok(getN(om, 'var') == vN, sprintf('%s : var.N  = %d', t, vN));
            t_ok(get(om, 'var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));
        end
    end
end

%%-----  getN  -----
t = 'getN(om, ''var'', ''Pg'') == 3';
t_ok(getN(om, 'var', 'Pg') == 3, t);

t = 'size(getN(om, ''var'', ''x'')) == [2 2]';
t_is(size(getN(om, 'var', 'x')), [2,2], 1e-14, t);

t = 'getN(om, ''var'', ''x'')(1,2) == 2';
N = getN(om, 'var', 'x');
t_is(N(1,2), 2, 1e-14, t);

t = 'getN(om, ''var'', ''x'', {2,1}) == 3';
t_is(getN(om, 'var', 'x', {2,1}), 3, 1e-14, t);

t = 'getN(om, ''var'', ''y'', {2,1,3}) == 6';
t_is(getN(om, 'var', 'y', {2,1,3}), 6, 1e-14, t);

t = 'getN(om, ''var'')';
t_is(getN(om, 'var'), vN, 1e-14, t);

%%-----  get_idx  -----
t = 'get_idx : var';
vv = get_idx(om);
t_is([vv.i1.Pg vv.iN.Pg vv.N.Pg], [5 7 3], 1e-14, [t ' : Pg']);
t_is(size(vv.i1.x), [2, 2], 1e-14, [t ' : size(vv.i1.x)']);
t_is([vv.i1.x(2,1) vv.iN.x(2,1) vv.N.x(2,1)], [22 24 3], 1e-14, [t ' : x(2,1)']);
t_is(size(vv.i1.y), [2, 3, 4], 1e-14, [t ' : size(vv.i1.y)']);
t_is([vv.i1.y(2,2,4) vv.iN.y(2,2,4) vv.N.y(2,2,4)], [133 140 8], 1e-14, [t ' : y(2,2,4)']);

%%-----  getv  -----
t = 'getv(om, ''Va'')';
[v0, vl, vu] = getv(om, 'Va');
t_ok(~any(v0), [t ' : v0']);
t_ok(all(isinf(vl) & vl < 0), [t ' : vl']);
t_ok(all(isinf(vu) & vu > 0), [t ' : vu']);

t = 'getv(om, ''Pg'')';
[v0, vl, vu] = getv(om, 'Pg');
t_is(v0, [2;4;6], 1e-14, [t ' : v0']);
t_is(vl, [1;2;3], 1e-14, [t ' : vl']);
t_is(vu, [10;20;30], 1e-14, [t ' : vu']);

t = 'getv(om, ''Vm1'')';
[v0, vl, vu, vt] = getv(om, 'Vm1');
t_is(double(vt), double('I'), 1e-14, [t ' : vt']);

t = 'getv(om, ''Vm2'')';
[v0, vl, vu, vt] = getv(om, 'Vm2');
t_is(double(vt), double('CIBIC'), 1e-14, [t ' : vt']);

t = 'getv(om, ''x'')';
[v0, vl, vu, vt] = getv(om, 'x');
t_is(size(v0), [2,2], 1e-14, [t ' : size(v0)']);
t_is(v0{2,2}, [1;0], 1e-14, [t ' : v0{2,2}']);
t_is(vl{2,2}, [0;-1], 1e-14, [t ' : vl{2,2}']);
t_is(vu{2,2}, [2;1], 1e-14, [t ' : vu{2,2}']);
t_is(double(vt{2,2}), double('C'), 1e-14, [t ' : vt{2,2}']);

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
            t = sprintf('getv(om, ''y'', {%d,%d,%d})', i, j, k);
            [v0, vl, vu, gvt] = getv(om, 'y', {i,j,k});
            t_is(v0, 10*(n:-1:1)', 1e-14, [t ' : v0']);
            t_is(vl, -1*(n:-1:1)', 1e-14, [t ' : vl']);
            t_is(vu, 100+(n:-1:1)', 1e-14, [t ' : vu']);
            t_is(gvt, vt, 1e-14, [t ' : vt']);
        end
    end
end

t = 'getv(om)';
[v0, vl, vu, vt] = getv(om);
t_ok(length(v0) == getN(om, 'var'), [t ' : length(v0)']);
t_ok(length(vl) == getN(om, 'var'), [t ' : length(vl)']);
t_ok(length(vu) == getN(om, 'var'), [t ' : length(vu)']);
t_is(v0(vv.i1.x(2,2):vv.iN.x(2,2)), [1;0], 1e-14, [t ' : v0(vv.i1.x(2,2):vv.iN.x(2,2))']);
t_is(vl(vv.i1.x(2,2):vv.iN.x(2,2)), [0;-1], 1e-14, [t ' : vl(vv.i1.x(2,2):vv.iN.x(2,2))']);
t_is(vu(vv.i1.x(2,2):vv.iN.x(2,2)), [2;1], 1e-14, [t ' : vu(vv.i1.x(2,2):vv.iN.x(2,2))']);
t_is(vt(vv.i1.x(2,2):vv.iN.x(2,2)), 'C', 1e-14, [t ' : vt(vv.i1.x(2,2):vv.iN.x(2,2))']);
n = 8;
t_is(v0(vv.i1.y(2,2,4):vv.iN.y(2,2,4)), 10*(n:-1:1)', 1e-14, [t ' : v0(vv.i1.y(2,2,4):vv.iN.y(2,2,4))']);
t_is(vl(vv.i1.y(2,2,4):vv.iN.y(2,2,4)), -1*(n:-1:1)', 1e-14, [t ' : vl(vv.i1.y(2,2,4):vv.iN.y(2,2,4))']);
t_is(vu(vv.i1.y(2,2,4):vv.iN.y(2,2,4)), 100+(n:-1:1)', 1e-14, [t ' : vu(vv.i1.y(2,2,4):vv.iN.y(2,2,4))']);
t_is(vt(vv.i1.y(2,2,4):vv.iN.y(2,2,4)), 'IIBIIIII', 1e-14, [t ' : vt(vv.i1.y(2,2,4):vv.iN.y(2,2,4))']);
vt0 = 'CCCCCCCIIIIICIBICCCCCCCCCCCCCCCCCCCCCCCCCCCCIIIIIIIIIIIIIIIIIIIIIIBBBBBBBBBBBBBBBBBBBBBBBBBBCICCCICCCCICCCCCICCCCCIIBIIIIBIIIIIBIIIIIIBIIIIIBBBCBBBBBCBBBBBBCBBBBBBBCBBBBB';
t_is(vt, vt0, 1e-14, [t ' : vt']);

%%-----  add_constraints  -----
t = 'add_constraints';
lN = 0;
lNS = 0;
t_ok(getN(om, 'lin') == lN, sprintf('%s : lin.N  = %d', t, lN));
t_ok(get(om, 'lin', 'NS') == lNS, sprintf('%s : lin.NS = %d', t, lNS));

t = 'add_constraints(om, ''Pmis'', A, l, u, {''Va'', ''Pg''})';
A = sparse([1:3 1:3 1:3]', [1:3 4:6 7 7 7]', [1 1 1 -1 -1 -1 2 3 4]', 3, 7);
l = -(1:3)'; u = (1:3)';
om = add_constraints(om, 'Pmis', A, l, u, {'Va', 'Pg'});
lNS = lNS + 1; lN = lN + 3;
t_ok(getN(om, 'lin') == lN, sprintf('%s : lin.N  = %d', t, lN));
t_ok(get(om, 'lin', 'NS') == lNS, sprintf('%s : lin.NS = %d', t, lNS));

t = 'add_constraints(om, ''Qmis'', A, l, u)';
A = sparse([1:3 1:3 1:3]', [1:3 4:6 7 7 7]', [1 1 1 -1 -1 -1 2 3 4]', 3, vN);
om = add_constraints(om, 'Qmis', A, l, u);
lNS = lNS + 1; lN = lN + 3;
t_ok(getN(om, 'lin') == lN, sprintf('%s : lin.N  = %d', t, lN));
t_ok(get(om, 'lin', 'NS') == lNS, sprintf('%s : lin.NS = %d', t, lNS));

t = 'add_constraints(om, ''mylin'', {2, 2})';
om = add_constraints(om, 'mylin', {2, 2});
t_ok(getN(om, 'lin') == lN, sprintf('%s : lin.N  = %d', t, lN));
t_ok(get(om, 'lin', 'NS') == lNS, sprintf('%s : lin.NS = %d', t, lNS));

for i = 1:2
    for j = 1:2
        t = sprintf('add_constraints(om, ''mylin'', {%d,%d}, A, l, u, vs)', i,j);
        A = sparse([1:(i+j) 1:(i+j)]', [1:(i+j) 5*ones(1,i+j)]', ...
            [ones(i+j,1);-ones(i+j,1)], i+j, 3+2+(i==2 && j==1));
        l = -ones(i+j, 1); u = [];
        vs = struct('name', {'Pg', 'x'}, 'idx', {{}, {i,j}});
        om = add_constraints(om, 'mylin', {i, j}, A, l, u, vs);
        lNS = lNS + 1; lN = lN + i+j;
        t_ok(getN(om, 'lin') == lN, sprintf('%s : lin.N  = %d', t, lN));
        t_ok(get(om, 'lin', 'NS') == lNS, sprintf('%s : lin.NS = %d', t, lNS));
    end
end

%%-----  get_idx  -----
t = 'get_idx : lin';
[vv, ll] = get_idx(om);
t_is([ll.i1.Qmis ll.iN.Qmis ll.N.Qmis], [4 6 3], 1e-14, [t ' : Qmis']);
t_is(size(ll.i1.mylin), [2, 2], 1e-14, [t ' : size(ll.i1.mylin)']);
t_is([ll.i1.mylin(2,1) ll.iN.mylin(2,1) ll.N.mylin(2,1)], [12 14 3], 1e-14, [t ' : mylin(2,1)']);

%%-----  linear_constraints  -----
t = 'linear_constraints';
[A, l, u] = linear_constraints(om);
t_ok(issparse(A), [t ' : issparse(A)']);
t_is(size(A), [lN, vN], 1e-14, [t ' : size(A)']);
t_is(length(l), lN, 1e-14, [t ' : length(l)']);
t_is(length(u), lN, 1e-14, [t ' : length(u)']);
AA = sparse([1:3 1:3 1:3]', [1:3 4:6 7 7 7]', [1 1 1 -1 -1 -1 2 3 4]', 3, vN);
t_is(full(A(ll.i1.Qmis:ll.iN.Qmis, :)), full(AA), 1e-14, [t ' : A(<Qmis>,:)']);
t_is(full(A(ll.i1.mylin(2,1):ll.iN.mylin(2,1), vv.i1.Pg:vv.iN.Pg)), eye(3,3), 1e-14, [t ' : A(<mylin(2,1)>,<Pg>)']);
t_is(full(A(ll.i1.mylin(2,1):ll.iN.mylin(2,1), vv.i1.x(2,1):vv.iN.x(2,1))), [0 -1 0;0 -1 0;0 -1 0], 1e-14, [t ' : A(<mylin(2,1)>,<x(2,1)>)']);

%%-----  add_costs  -----
t = 'add_costs';
cN = 0;
cNS = 0;
t_ok(getN(om, 'cost') == cN, sprintf('%s : cost.N  = %d', t, cN));
t_ok(get(om, 'cost', 'NS') == cNS, sprintf('%s : cost.NS = %d', t, cNS));

t = 'add_costs(om, ''ucost'', cp, {''Va'', ''Pg''})';
cp = struct('N', sparse([1:2 1:2 1:2]', [1:4 5 7]', [1 1 -1 -1 2 2]', 2,7), ...
            'Cw', [2;3]);
om = add_costs(om, 'ucost', cp, {'Va', 'Pg'});
cNS = cNS + 1; cN = cN + 2;
t_ok(getN(om, 'cost') == cN, sprintf('%s : cost.N  = %d', t, cN));
t_ok(get(om, 'cost', 'NS') == cNS, sprintf('%s : cost.NS = %d', t, cNS));

t = 'add_costs(om, ''vcost'', cp)';
cp = struct('N', sparse([1:2 1:2 1:2]', [1:4 5 7]', [1 1 -1 -1 2 2]', 2, vN), ...
            'Cw', [2;3]);
om = add_costs(om, 'vcost', cp);
cNS = cNS + 1; cN = cN + 2;
t_ok(getN(om, 'cost') == cN, sprintf('%s : cost.N  = %d', t, cN));
t_ok(get(om, 'cost', 'NS') == cNS, sprintf('%s : cost.NS = %d', t, cNS));

t = 'add_costs(om, ''wc'', {2,2})';
om = add_costs(om, 'wc', {2,2});
t_ok(getN(om, 'cost') == cN, sprintf('%s : cost.N  = %d', t, cN));
t_ok(get(om, 'cost', 'NS') == cNS, sprintf('%s : cost.NS = %d', t, cNS));

for i = 1:2
    for j = 1:2
        t = 'add_costs(om, ''wc'', {i, j}, cp, vs)';
        cp.N = sparse([1:(i+j) 1:(i+j)]', [1:(i+j) 5*ones(1,i+j)]', ...
            [ones(i+j,1);-ones(i+j,1)], i+j, 3+2+(i==2 && j==1));
        cp.Cw = (i+j:-1:1)';
        if i == 2
            cp.H = sparse((1:i+j)', (1:i+j)', (1:i+j)', i+j, i+j);
        end
        vs = struct('name', {'Pg', 'x'}, 'idx', {{}, {i,j}});
        om = add_costs(om, 'wc', {i, j}, cp, vs);
        cNS = cNS + 1; cN = cN + i+j;
        t_ok(getN(om, 'cost') == cN, sprintf('%s : cost.N  = %d', t, cN));
        t_ok(get(om, 'cost', 'NS') == cNS, sprintf('%s : cost.NS = %d', t, cNS));
    end
end

%%-----  get_cost_params/build_cost_params  -----
t = 'get_cost_params : error';
try
    cp = get_cost_params(om);
    t_ok(0, t);
catch
    t_ok(strfind(lasterr, '@opt_model/get_cost_params: must call build_cost_params first'), t);
end

t = 'build_cost_params';
om = build_cost_params(om);
cp = get_cost_params(om);
t_ok(isfield(cp, 'N'), t);

%%-----  get_idx  -----
t = 'get_idx : cost';
[vv, ll, nn, cc] = get_idx(om);
t_is([cc.i1.vcost cc.iN.vcost cc.N.vcost], [3 4 2], 1e-14, [t ' : vcost']);
t_is(size(cc.i1.wc), [2, 2], 1e-14, [t ' : size(cc.i1.wc)']);
t_is([cc.i1.wc(2,1) cc.iN.wc(2,1) cc.N.wc(2,1)], [10 12 3], 1e-14, [t ' : wc(2,1)']);

%%-----  get_cost_params  -----
t = 'get_cost_params(om, ''ucost'')';
cp = get_cost_params(om, 'ucost');
N = sparse([1:2 1:2 1:2]', [1:4 5 7]', [1 1 -1 -1 2 2]', 2, vN);
t_is(full(cp.N), full(N), 1e-14, [t, ' : N']);
t_is(cp.Cw, [2;3], 1e-14, [t, ' : Cw']);
t_is(full(cp.H), zeros(2,2), 1e-14, [t, ' : H']);
t_is(cp.dd, ones(2,1),  1e-14, [t, ' : dd']);
t_is(cp.rh, zeros(2,1), 1e-14, [t, ' : rh']);
t_is(cp.kk, zeros(2,1), 1e-14, [t, ' : kk']);
t_is(cp.mm, ones(2,1),  1e-14, [t, ' : mm']);

t = 'get_cost_params(om, ''vcost'')';
cp = get_cost_params(om, 'vcost');
N = sparse([1:2 1:2 1:2]', [1:4 5 7]', [1 1 -1 -1 2 2]', 2, vN);
t_is(full(cp.N), full(N), 1e-14, [t, ' : N']);
t_is(cp.Cw, [2;3], 1e-14, [t, ' : Cw']);
t_is(full(cp.H), zeros(2,2), 1e-14, [t, ' : H']);
t_is(cp.dd, ones(2,1),  1e-14, [t, ' : dd']);
t_is(cp.rh, zeros(2,1), 1e-14, [t, ' : rh']);
t_is(cp.kk, zeros(2,1), 1e-14, [t, ' : kk']);
t_is(cp.mm, ones(2,1),  1e-14, [t, ' : mm']);

t = 'get_cost_params(om, ''wc'') : error';
try
    cp = get_cost_params(om, 'wc')
    t_ok(0, t);
catch
    t_ok(strfind(lasterr, '@opt_model/get_cost_params: cost set ''wc'' requires an idx arg'), t);
end

t = 'get_cost_params(om, ''wc'', {1,2})';
cp = get_cost_params(om, 'wc', {1,2});
N = sparse([1:3 1:3]', [vv.i1.Pg-1+(1:3) vv.i1.x(1,2)+ones(1,3)]', [ones(3,1);-ones(3,1)], 3, vN);
t_is(full(cp.N), full(N), 1e-14, [t, ' : N']);
t_is(cp.Cw, [3;2;1], 1e-14, [t, ' : Cw']);
t_is(full(cp.H), zeros(3,3), 1e-14, [t, ' : H']);
t_is(cp.dd, ones(3,1),  1e-14, [t, ' : dd']);
t_is(cp.rh, zeros(3,1), 1e-14, [t, ' : rh']);
t_is(cp.kk, zeros(3,1), 1e-14, [t, ' : kk']);
t_is(cp.mm, ones(3,1),  1e-14, [t, ' : mm']);

t = 'get_cost_params(om, ''wc'', {2,1})';
cp = get_cost_params(om, 'wc', {2,1});
N = sparse([1:3 1:3]', [vv.i1.Pg-1+(1:3) vv.i1.x(2,1)+ones(1,3)]', [ones(3,1);-ones(3,1)], 3, vN);
t_is(full(cp.N), full(N), 1e-14, [t, ' : N']);
t_is(cp.Cw, [3;2;1], 1e-14, [t, ' : Cw']);
H = sparse(1:3, 1:3, 1:3, 3, 3);
t_is(full(cp.H), full(H), 1e-14, [t, ' : H']);
t_is(cp.dd, ones(3,1),  1e-14, [t, ' : dd']);
t_is(cp.rh, zeros(3,1), 1e-14, [t, ' : rh']);
t_is(cp.kk, zeros(3,1), 1e-14, [t, ' : kk']);
t_is(cp.mm, ones(3,1),  1e-14, [t, ' : mm']);

t = 'get_cost_params(om)';
cp = get_cost_params(om);
t_ok(issparse(cp.N), [t ' : issparse(cp.N)']);
t_is(size(cp.N), [cN, vN], 1e-14, [t ' : size(cp.N)']);
t_is(size(cp.H), [cN, cN], 1e-14, [t ' : size(cp.H)']);
t_is(length(cp.Cw), cN, 1e-14, [t ' : length(cp.Cw)']);
t_is(length(cp.dd), cN, 1e-14, [t ' : length(cp.dd)']);
t_is(length(cp.rh), cN, 1e-14, [t ' : length(cp.rh)']);
t_is(length(cp.kk), cN, 1e-14, [t ' : length(cp.kk)']);
t_is(length(cp.mm), cN, 1e-14, [t ' : length(cp.mm)']);
N = sparse([1:2 1:2]', [1:4]', [1 1 -1 -1]', 2, 4);
Cw = [2;3];
H = zeros(2,2);
t_is(full(cp.N(cc.i1.vcost:cc.iN.vcost, vv.i1.Va:vv.iN.Va)), full(N), 1e-14, [t ' : N(<vcost>,<Va>)']);
N = sparse([1:2]', [1 3]', [2 2]', 2, 3);
t_is(full(cp.N(cc.i1.vcost:cc.iN.vcost, vv.i1.Pg:vv.iN.Pg)), full(N), 1e-14, [t ' : N(<vcost>,<Pg>)']);
t_is(full(cp.Cw(cc.i1.vcost:cc.iN.vcost)), Cw, 1e-14, [t ' : Cw(<vcost>)']);
t_is(full(cp.H(cc.i1.vcost:cc.iN.vcost, cc.i1.vcost:cc.iN.vcost)), full(H), 1e-14, [t ' : H(<vcost>,<vcost>)']);
N = sparse([1:3]', [1:3]', [1 1 1]', 3, 3);
t_is(full(cp.N(cc.i1.wc(1,2):cc.iN.wc(1,2), vv.i1.Pg:vv.iN.Pg)), full(N), 1e-14, [t ' : N(<wc(1,2)>,<Pg>)']);
N = sparse([1:3]', [2 2 2]', [-1 -1 -1]', 3, 2);
Cw = [3;2;1];
H = zeros(3,3);
t_is(full(cp.N(cc.i1.wc(1,2):cc.iN.wc(1,2), vv.i1.x(1,2):vv.iN.x(1,2))), full(N), 1e-14, [t ' : N(<wc(1,2)>,<x(1,2)>)']);
t_is(full(cp.Cw(cc.i1.wc(1,2):cc.iN.wc(1,2))), Cw, 1e-14, [t ' : Cw(<wc(1,2)>)']);
t_is(full(cp.H(cc.i1.wc(1,2):cc.iN.wc(1,2), cc.i1.wc(1,2):cc.iN.wc(1,2))), full(H), 1e-14, [t ' : H(<wc(1,2)>,<wc(1,2)>)']);

N = sparse([1:3]', [1:3]', [1 1 1]', 3, 3);
t_is(full(cp.N(cc.i1.wc(2,1):cc.iN.wc(2,1), vv.i1.Pg:vv.iN.Pg)), full(N), 1e-14, [t ' : N(<wc(2,1)>,<Pg>)']);
N = sparse([1:3]', [2 2 2]', [-1 -1 -1]', 3, 3);
Cw = [3;2;1];
H = full(sparse((1:3)', (1:3)', (1:3)', 3, 3));
t_is(full(cp.N(cc.i1.wc(2,1):cc.iN.wc(2,1), vv.i1.x(2,1):vv.iN.x(2,1))), full(N), 1e-14, [t ' : N(<wc(2,1)>,<x(2,1)>)']);
t_is(full(cp.Cw(cc.i1.wc(2,1):cc.iN.wc(2,1))), Cw, 1e-14, [t ' : Cw(<wc(2,1)>)']);
t_is(full(cp.H(cc.i1.wc(2,1):cc.iN.wc(2,1), cc.i1.wc(2,1):cc.iN.wc(2,1))), full(H), 1e-14, [t ' : H(<wc(2,1)>,<wc(2,1)>)']);

t_is(cp.dd, ones(cN,1),  1e-14, [t, ' : dd']);
t_is(cp.rh, zeros(cN,1), 1e-14, [t, ' : rh']);
t_is(cp.kk, zeros(cN,1), 1e-14, [t, ' : kk']);
t_is(cp.mm, ones(cN,1),  1e-14, [t, ' : mm']);

%%-----  compute_cost  -----
% f = compute_cost(om, x, name, idx)
t = 'compute_cost(om, x)';
x = [1:7 rand(1,10) 8:(vN-10)]';
f = compute_cost(om, x);
t_is(f, 343, 1e-14, t);

t = 'compute_cost(om, ''ucost'')';
f = compute_cost(om, x, 'ucost');
t_is(f, 52, 1e-14, t);

t = 'compute_cost(om, ''wc'', {2,1})';
f = compute_cost(om, x, 'wc', {2,1});
t_is(f, 91, 1e-14, t);

% om
% om = struct(om);

t_end;
