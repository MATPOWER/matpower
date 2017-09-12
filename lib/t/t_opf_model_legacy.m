function t_opf_model_legacy(quiet)
%T_OPF_MODEL_LEGACY Tests for OPF_MODEL, using the deprecated ADD_CONSTRAINTS.

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

num_tests = 363;

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
t_ok(om.getN('var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(om.get('var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'om.add_vars(''Va'', 4)';
om.add_vars('Va', 4);
vNS = vNS + 1; vN = vN + 4;
t_ok(om.getN('var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(om.get('var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'om.add_vars(''Pg'', 3, Pg0, Pgmin, Pgmax)';
om.add_vars('Pg', 3, [2;4;6], [1;2;3], [10;20;30]);
vNS = vNS + 1; vN = vN + 3;
t_ok(om.getN('var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(om.get('var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'om.add_vars(''Vm1'', 5, V0, Vmin, Vmax, ''I'')';
V0 = [1;1;1;1;1];
Vmin = zeros(5, 1);
Vmax = 1 + 0.01*(1:5)';
vt = 'I';
om.add_vars('Vm1', 5, V0, Vmin, Vmax, vt);
vNS = vNS + 1; vN = vN + 5;
t_ok(om.getN('var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(om.get('var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'om.add_vars(''Vm2'', 5, V0, Vmin, Vmax, ''CIBIC'')';
vt = 'CIBIC';
om.add_vars('Vm2', 5, V0, Vmin, Vmax, vt);
vNS = vNS + 1; vN = vN + 5;
t_ok(om.getN('var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(om.get('var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'om.add_vars(''x'', dims)';
om.add_vars('x', {2,2});
t_ok(om.getN('var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(om.get('var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'om.add_vars(''x'', {1,1}, 2)';
om.add_vars('x', {1,1}, 2);
vNS = vNS + 1; vN = vN + 2;
t_ok(om.getN('var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(om.get('var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'om.add_vars(''x'', {1,2}, 2, x0(1,2))';
om.add_vars('x', {1,2}, 2, [-1;-2]);
vNS = vNS + 1; vN = vN + 2;
t_ok(om.getN('var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(om.get('var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'om.add_vars(''x'', {2,1}, 3)';
om.add_vars('x', {2,1}, 3);
vNS = vNS + 1; vN = vN + 3;
t_ok(om.getN('var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(om.get('var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'om.add_vars(''x'', {2,2}, 2, x0(2,2), xmin(2,2), xmax(2,2))';
om.add_vars('x', {2,2}, 2, [1;0],[0;-1],[2;1]);
vNS = vNS + 1; vN = vN + 2;
t_ok(om.getN('var') == vN, sprintf('%s : var.N  = %d', t, vN));
t_ok(om.get('var', 'NS') == vNS, sprintf('%s : var.NS = %d', t, vNS));

t = 'om.add_vars(''y'', {2,3,4})';
om.add_vars('y', {2,3,4});
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
            t = sprintf('om.add_vars(''y'', {%d,%d,%d}, y0, ymin, ymax, vt)', i,j,k);
            om.add_vars('y', {i,j,k}, n, 10*(n:-1:1)', -1*(n:-1:1)', 100+(n:-1:1)', vt);
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

%%-----  getv  -----
t = 'om.getv(''Va'')';
[v0, vl, vu] = om.getv('Va');
t_ok(~any(v0), [t ' : v0']);
t_ok(all(isinf(vl) & vl < 0), [t ' : vl']);
t_ok(all(isinf(vu) & vu > 0), [t ' : vu']);

t = 'om.getv(''Pg'')';
[v0, vl, vu] = om.getv('Pg');
t_is(v0, [2;4;6], 14, [t ' : v0']);
t_is(vl, [1;2;3], 14, [t ' : vl']);
t_is(vu, [10;20;30], 14, [t ' : vu']);

t = 'om.getv(''Vm1'')';
[v0, vl, vu, vt] = om.getv('Vm1');
t_is(double(vt), double('I'), 14, [t ' : vt']);

t = 'om.getv(''Vm2'')';
[v0, vl, vu, vt] = om.getv('Vm2');
t_is(double(vt), double('CIBIC'), 14, [t ' : vt']);

t = 'om.getv(''x'')';
[v0, vl, vu, vt] = om.getv('x');
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
            t = sprintf('om.getv(''y'', {%d,%d,%d})', i, j, k);
            [v0, vl, vu, gvt] = om.getv('y', {i,j,k});
            t_is(v0, 10*(n:-1:1)', 14, [t ' : v0']);
            t_is(vl, -1*(n:-1:1)', 14, [t ' : vl']);
            t_is(vu, 100+(n:-1:1)', 14, [t ' : vu']);
            t_is(gvt, vt, 14, [t ' : vt']);
        end
    end
end

t = 'om.getv()';
[v0, vl, vu, vt] = om.getv();
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

%%-----  add_constraints (linear)  -----
t = 'add_constraints';
lN = 0;
lNS = 0;
t_ok(om.getN('lin') == lN, sprintf('%s : lin.N  = %d', t, lN));
t_ok(om.get('lin', 'NS') == lNS, sprintf('%s : lin.NS = %d', t, lNS));

t = 'om.add_constraints(''Pmis'', A, l, u, {''Va'', ''Pg''})';
A = sparse([1:3 1:3 1:3]', [1:3 4:6 7 7 7]', [1 1 1 -1 -1 -1 2 3 4]', 3, 7);
l = -(1:3)'; u = (1:3)';
om.add_constraints('Pmis', A, l, u, {'Va', 'Pg'});
lNS = lNS + 1; lN = lN + 3;
t_ok(om.getN('lin') == lN, sprintf('%s : lin.N  = %d', t, lN));
t_ok(om.get('lin', 'NS') == lNS, sprintf('%s : lin.NS = %d', t, lNS));

t = 'om.add_constraints(''Qmis'', A, l, u)';
A = sparse([1:3 1:3 1:3]', [1:3 4:6 7 7 7]', [1 1 1 -1 -1 -1 2 3 4]', 3, vN);
om.add_constraints('Qmis', A, l, u);
lNS = lNS + 1; lN = lN + 3;
t_ok(om.getN('lin') == lN, sprintf('%s : lin.N  = %d', t, lN));
t_ok(om.get('lin', 'NS') == lNS, sprintf('%s : lin.NS = %d', t, lNS));

t = 'om.add_constraints(''mylin'', {2, 2})';
om.add_constraints('mylin', {2, 2});
% om.add_constraints('mylin', {2, 2}, 'lin');
t_ok(om.getN('lin') == lN, sprintf('%s : lin.N  = %d', t, lN));
t_ok(om.get('lin', 'NS') == lNS, sprintf('%s : lin.NS = %d', t, lNS));

for i = 1:2
    for j = 1:2
        t = sprintf('om.add_constraints(''mylin'', {%d,%d}, A, l, u, vs)', i,j);
        A = sparse([1:(i+j) 1:(i+j)]', [1:(i+j) 5*ones(1,i+j)]', ...
            [ones(i+j,1);-ones(i+j,1)], i+j, 3+2+(i==2 && j==1));
        l = -ones(i+j, 1); u = [];
        vs = struct('name', {'Pg', 'x'}, 'idx', {{}, {i,j}});
        om.add_constraints('mylin', {i, j}, A, l, u, vs);
        lNS = lNS + 1; lN = lN + i+j;
        t_ok(om.getN('lin') == lN, sprintf('%s : lin.N  = %d', t, lN));
        t_ok(om.get('lin', 'NS') == lNS, sprintf('%s : lin.NS = %d', t, lNS));
    end
end

%%-----  add_constraints (nonlinear, equality)  -----
t = 'add_constraints (nonlinear, equality)';
neN = 0;
neNS = 0;
t_ok(om.getN('nle') == neN, sprintf('%s : nle.N  = %d', t, neN));
t_ok(om.get('nle', 'NS') == neNS, sprintf('%s : nle.NS = %d', t, neNS));

t = 'om.add_constraints(''Pmise'', N, 1, fcn, hess, {''Pg'', ''Va''})';
N = 4;
fcn = @(x)my_fcn(x, N, 2);
hess = @(x, lam)my_hess(x, lam, 10);
om.add_constraints('Pmise', N, 1, fcn, hess, {'Pg', 'Va'});
neNS = neNS + 1; neN = neN + N;
t_ok(om.getN('nle') == neN, sprintf('%s : nle.N  = %d', t, neN));
t_ok(om.get('nle', 'NS') == neNS, sprintf('%s : nle.NS = %d', t, neNS));

t = 'om.add_constraints(''Qmise'', N, 1, fcn, hess)';
N = 3;
fcn = @(x)my_fcn(x, N, 2);
hess = @(x, lam)my_hess(x, lam, 10);
om.add_constraints('Qmise', N, 1, fcn, hess);
neNS = neNS + 1; neN = neN + N;
t_ok(om.getN('nle') == neN, sprintf('%s : nle.N  = %d', t, neN));
t_ok(om.get('nle', 'NS') == neNS, sprintf('%s : nle.NS = %d', t, neNS));

t = 'om.add_constraints(''mynle'', {2, 2})';
om.add_constraints('mynle', {2, 2}, 'nle');
t_ok(om.getN('nle') == neN, sprintf('%s : nle.N  = %d', t, neN));
t_ok(om.get('nle', 'NS') == neNS, sprintf('%s : nle.NS = %d', t, neNS));

for i = 1:2
    for j = 1:2
        t = sprintf('om.add_constraints(''mynle'', {%d,%d}, N, 1, fcn, hess, vs)', i,j);
        N = i+j;
        fcn = @(x)my_fcn(x, N, i);
        hess = @(x, lam)my_hess(x, lam, j);
        vs = struct('name', {'Pg', 'x'}, 'idx', {{}, {i,j}});
        om.add_constraints('mynle', {i, j}, N, 1, fcn, hess, vs);
        neNS = neNS + 1; neN = neN + N;
        t_ok(om.getN('nle') == neN, sprintf('%s : nle.N  = %d', t, neN));
        t_ok(om.get('nle', 'NS') == neNS, sprintf('%s : nle.NS = %d', t, neNS));
    end
end

%%-----  add_constraints (nonlinear, inequality)  -----
t = 'add_constraints (nonlinear, inequality)';
niN = 0;
niNS = 0;
t_ok(om.getN('nli') == niN, sprintf('%s : nli.N  = %d', t, niN));
t_ok(om.get('nli', 'NS') == niNS, sprintf('%s : nli.NS = %d', t, niNS));

t = 'om.add_constraints(''Pmisi'', N, 0, fcn, hess, {''Pg'', ''Va''})';
N = 3;
fcn = @(x)my_fcn(x, N, -2);
hess = @(x, lam)my_hess(x, lam, -10);
om.add_constraints('Pmisi', N, 0, fcn, hess, {'Pg', 'Va'});
niNS = niNS + 1; niN = niN + N;
t_ok(om.getN('nli') == niN, sprintf('%s : nli.N  = %d', t, niN));
t_ok(om.get('nli', 'NS') == niNS, sprintf('%s : nli.NS = %d', t, niNS));

t = 'om.add_constraints(''Qmisi'', N, 0, fcn, hess)';
N = 2;
fcn = @(x)my_fcn(x, N, -2);
hess = @(x, lam)my_hess(x, lam, -10);
om.add_constraints('Qmisi', N, 0, fcn, hess);
niNS = niNS + 1; niN = niN + N;
t_ok(om.getN('nli') == niN, sprintf('%s : nli.N  = %d', t, niN));
t_ok(om.get('nli', 'NS') == niNS, sprintf('%s : nli.NS = %d', t, niNS));

t = 'om.add_constraints(''mynli'', {2, 2})';
om.add_constraints('mynli', {2, 2}, 'nli');
t_ok(om.getN('nli') == niN, sprintf('%s : nli.N  = %d', t, niN));
t_ok(om.get('nli', 'NS') == niNS, sprintf('%s : nli.NS = %d', t, niNS));

for i = 1:2
    for j = 1:2
        t = sprintf('om.add_constraints(''mynli'', {%d,%d}, N, 0, fcn, hess, vs)', i,j);
        N = i+j-1;
        fcn = @(x)my_fcn(x, N, i);
        hess = @(x, lam)my_hess(x, lam, j);
        vs = struct('name', {'Pg', 'x'}, 'idx', {{}, {i,j}});
        om.add_constraints('mynli', {i, j}, N, 0, fcn, hess, vs);
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

%%-----  linear_constraints  -----
t = 'linear_constraints';
[A, l, u] = om.linear_constraints();
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
x = (1:om.getN('var'))';
[g, dg] = om.eval_nln_constraint(x, 1);
t_is(length(g), neN, 14, [t ' : length(g)']);
t_ok(issparse(dg), [t ' : issparse(dg)']);
t_is(size(dg), [neN, vN], 14, [t ' : size(dg)']);
t_is(g, [7 8 9 3 3 4 5 6 7 6 7 8 7 8 9 7 8 9 27]', 14, [t ' : g']);
e = [[  1 2 3 4 7 6 7;
        0 0 0 0 0 2 0;
        0 0 0 0 0 0 2;
        2 0 0 0 0 0 0 ] zeros(4, vN-7) ];
t_is(full(dg(1:4, :)), e, 14, [t ' : dg(1:4, :)   [Pmise]']);
e = [[3 2:vN]; [0 2 0 zeros(1, vN-3)]; [0 0 2 zeros(1, vN-3)]];
t_is(full(dg(5:7, :)), e, 14, [t ' : dg(5:7, :)   [Qmise]']);
e = [[0 0 0 0 6 6 7; 0 0 0 0 0 1 0] zeros(2, 10) [18 19; 0 0] zeros(2, vN-19)];
t_is(full(dg(8:9, :)), e, 14, [t ' : dg(8:9, :)   [mynle(1,1)]']);
e = [[0 0 0 0 6 6 7; 0 0 0 0 0 1 0; 0 0 0 0 0 0 1] zeros(3, 12) [20 21; 0 0; 0 0] zeros(3, vN-21)];
t_is(full(dg(10:12, :)), e, 14, [t ' : dg(10:12, :) [mynle(1,2)]']);
e = [[0 0 0 0 7 6 7; 0 0 0 0 0 2 0; 0 0 0 0 0 0 2] zeros(3, 14) [22 23 24; 0 0 0; 0 0 0] zeros(3, vN-24)];
t_is(full(dg(13:15, :)), e, 14, [t ' : dg(13:15, :) [mynle(2,1)]']);
e = [[0 0 0 0 7 6 7; 0 0 0 0 0 2 0; 0 0 0 0 0 0 2; 0 0 0 0 0 0 0] zeros(4, 17) [25 26; 0 0; 0 0; 2 0] zeros(4, vN-26)];
t_is(full(dg(16:19, :)), e, 14, [t ' : dg(16:19, :) [mynle(2,2)]']);
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
e = sparse([1:3 5:7 25], [1:3 5:7 25], [22.09 12.06 13.07 41.48 46.53 43.48 27.19], vN, vN);
t_is(d2G, e, 13, [t ' : d2G']);

%d2G

lam = -(1:niN)'/100;
d2H = om.eval_nln_constraint_hess(x, lam, 0);
t_ok(issparse(d2H), [t ' : issparse(d2H)']);
t_is(size(d2H), [vN, vN], 14, [t ' : size(d2H)']);
% t_is(full(d2H(27:end, :)), zeros(vN-26, vN), 14, [t ' : d2H(27:end, :)']);
% t_is(full(d2H(:, 27:end)), zeros(vN, vN-26), 14, [t ' : d2H(:, 27:end)']);
e = sparse([1:2 5:7], [1:2 5:7], [-9.04 -8.05 20.66 18.68 5.84], vN, vN);
t_is(d2H, e, 13, [t ' : d2H']);

%d2H

%%-----  add_costs  -----
t = 'add_costs';
cN = 0;
cNS = 0;
t_ok(om.getN('cost') == cN, sprintf('%s : cost.N  = %d', t, cN));
t_ok(om.get('cost', 'NS') == cNS, sprintf('%s : cost.NS = %d', t, cNS));

t = 'om.add_costs(''ucost'', cp, {''Va'', ''Pg''})';
cp = struct('N', sparse([1:2 1:2 1:2]', [1:4 5 7]', [1 1 -1 -1 2 2]', 2,7), ...
            'Cw', [2;3]);
om.add_costs('ucost', cp, {'Va', 'Pg'});
cNS = cNS + 1; cN = cN + 2;
t_ok(om.getN('cost') == cN, sprintf('%s : cost.N  = %d', t, cN));
t_ok(om.get('cost', 'NS') == cNS, sprintf('%s : cost.NS = %d', t, cNS));

t = 'om.add_costs(''vcost'', cp)';
cp = struct('N', sparse([1:2 1:2 1:2]', [1:4 5 7]', [1 1 -1 -1 2 2]', 2, vN), ...
            'Cw', [2;3]);
om.add_costs('vcost', cp);
cNS = cNS + 1; cN = cN + 2;
t_ok(om.getN('cost') == cN, sprintf('%s : cost.N  = %d', t, cN));
t_ok(om.get('cost', 'NS') == cNS, sprintf('%s : cost.NS = %d', t, cNS));

t = 'om.add_costs(''wc'', {2,2})';
om.add_costs('wc', {2,2});
t_ok(om.getN('cost') == cN, sprintf('%s : cost.N  = %d', t, cN));
t_ok(om.get('cost', 'NS') == cNS, sprintf('%s : cost.NS = %d', t, cNS));

for i = 1:2
    for j = 1:2
        t = 'om.add_costs(''wc'', {i, j}, cp, vs)';
        cp.N = sparse([1:(i+j) 1:(i+j)]', [1:(i+j) 5*ones(1,i+j)]', ...
            [ones(i+j,1);-ones(i+j,1)], i+j, 3+2+(i==2 && j==1));
        cp.Cw = (i+j:-1:1)';
        if i == 2
            cp.H = sparse((1:i+j)', (1:i+j)', (1:i+j)', i+j, i+j);
        end
        vs = struct('name', {'Pg', 'x'}, 'idx', {{}, {i,j}});
        om.add_costs('wc', {i, j}, cp, vs);
        cNS = cNS + 1; cN = cN + i+j;
        t_ok(om.getN('cost') == cN, sprintf('%s : cost.N  = %d', t, cN));
        t_ok(om.get('cost', 'NS') == cNS, sprintf('%s : cost.NS = %d', t, cNS));
    end
end

t = 'build_cost_params';
om.build_cost_params();
cp = om.get_cost_params();
t_ok(isfield(cp, 'N'), t);

%%-----  get_idx  -----
t = 'get_idx : cost';
[vv, ll, nne, nni, cc] = om.get_idx();
t_is([cc.i1.vcost cc.iN.vcost cc.N.vcost], [3 4 2], 14, [t ' : vcost']);
t_is(size(cc.i1.wc), [2, 2], 14, [t ' : size(cc.i1.wc)']);
t_is([cc.i1.wc(2,1) cc.iN.wc(2,1) cc.N.wc(2,1)], [10 12 3], 14, [t ' : wc(2,1)']);

%%-----  get_cost_params  -----
t = 'om.get_cost_params(''ucost'')';
cp = om.get_cost_params('ucost');
N = sparse([1:2 1:2 1:2]', [1:4 5 7]', [1 1 -1 -1 2 2]', 2, vN);
t_is(full(cp.N), full(N), 14, [t, ' : N']);
t_is(cp.Cw, [2;3], 14, [t, ' : Cw']);
t_is(full(cp.H), zeros(2,2), 14, [t, ' : H']);
t_is(cp.dd, ones(2,1),  14, [t, ' : dd']);
t_is(cp.rh, zeros(2,1), 14, [t, ' : rh']);
t_is(cp.kk, zeros(2,1), 14, [t, ' : kk']);
t_is(cp.mm, ones(2,1),  14, [t, ' : mm']);

t = 'om.get_cost_params(''vcost'')';
cp = om.get_cost_params('vcost');
N = sparse([1:2 1:2 1:2]', [1:4 5 7]', [1 1 -1 -1 2 2]', 2, vN);
t_is(full(cp.N), full(N), 14, [t, ' : N']);
t_is(cp.Cw, [2;3], 14, [t, ' : Cw']);
t_is(full(cp.H), zeros(2,2), 14, [t, ' : H']);
t_is(cp.dd, ones(2,1),  14, [t, ' : dd']);
t_is(cp.rh, zeros(2,1), 14, [t, ' : rh']);
t_is(cp.kk, zeros(2,1), 14, [t, ' : kk']);
t_is(cp.mm, ones(2,1),  14, [t, ' : mm']);

t = 'om.get_cost_params(''wc'') : error';
try
    cp = om.get_cost_params('wc')
    t_ok(0, t);
catch
    t_ok(strfind(lasterr, '@opt_model/get_cost_params: cost set ''wc'' requires an idx arg'), t);
end

t = 'om.get_cost_params(''wc'', {1,2})';
cp = om.get_cost_params('wc', {1,2});
N = sparse([1:3 1:3]', [vv.i1.Pg-1+(1:3) vv.i1.x(1,2)+ones(1,3)]', [ones(3,1);-ones(3,1)], 3, vN);
t_is(full(cp.N), full(N), 14, [t, ' : N']);
t_is(cp.Cw, [3;2;1], 14, [t, ' : Cw']);
t_is(full(cp.H), zeros(3,3), 14, [t, ' : H']);
t_is(cp.dd, ones(3,1),  14, [t, ' : dd']);
t_is(cp.rh, zeros(3,1), 14, [t, ' : rh']);
t_is(cp.kk, zeros(3,1), 14, [t, ' : kk']);
t_is(cp.mm, ones(3,1),  14, [t, ' : mm']);

t = 'om.get_cost_params(''wc'', {2,1})';
cp = om.get_cost_params('wc', {2,1});
N = sparse([1:3 1:3]', [vv.i1.Pg-1+(1:3) vv.i1.x(2,1)+ones(1,3)]', [ones(3,1);-ones(3,1)], 3, vN);
t_is(full(cp.N), full(N), 14, [t, ' : N']);
t_is(cp.Cw, [3;2;1], 14, [t, ' : Cw']);
H = sparse(1:3, 1:3, 1:3, 3, 3);
t_is(full(cp.H), full(H), 14, [t, ' : H']);
t_is(cp.dd, ones(3,1),  14, [t, ' : dd']);
t_is(cp.rh, zeros(3,1), 14, [t, ' : rh']);
t_is(cp.kk, zeros(3,1), 14, [t, ' : kk']);
t_is(cp.mm, ones(3,1),  14, [t, ' : mm']);

t = 'om.get_cost_params()';
cp = om.get_cost_params();
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

%%-----  compute_cost  -----
% f = om.compute_cost(x, name, idx)
t = 'om.compute_cost(x)';
x = [1:7 rand(1,10) 8:(vN-10)]';
f = om.compute_cost(x);
t_is(f, 343, 14, t);

t = 'om.compute_cost(''ucost'')';
f = om.compute_cost(x, 'ucost');
t_is(f, 52, 14, t);

t = 'om.compute_cost(''wc'', {2,1})';
f = om.compute_cost(x, 'wc', {2,1});
t_is(f, 91, 14, t);

t = 'om.compute_cost(''wc'')';
f = om.compute_cost(x, 'wc');
t_is(f, 239, 14, t);

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
