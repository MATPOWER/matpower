function t_opf_fmincon(quiet)
%T_OPF_FMINCON  Tests for FMINCON-based optimal power flow.

%   MATPOWER
%   Copyright (c) 2004-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

%% current mismatch, cartesian V
options = {
    {0, 0},
    {0, 1},
    {1, 0},
    {1, 1},
};

num_tests = 216;

t_begin(length(options)*num_tests, quiet);

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

casefile = 't_case9_opf';
if quiet
    verbose = 0;
else
    verbose = 0;
end

mpopt = mpoption('opf.violation', 1e-6);
mpopt = mpoption(mpopt, 'out.all', 0, 'verbose', verbose, 'opf.ac.solver', 'FMINCON');
mpopt = mpoption(mpopt, 'fmincon.tol_x', 1e-7, 'fmincon.tol_f', 1e-9);

%% use active-set method for MATLAB 7.6-7.9 (R2008a-R2009b)
vstr = have_fcn('matlab', 'vstr');
if strcmp(vstr, '7.6') || strcmp(vstr, '7.7') || ...
        strcmp(vstr, '7.8') || strcmp(vstr, '7.9')
    mpopt = mpoption(mpopt, 'fmincon.alg', 1);
end

for k = 1:length(options)
    if options{k}{1}, bal = 'I';  else, bal = 'S'; end  %% nodal balance
    if options{k}{2}, crd = 'c';  else, crd = 'p'; end  %% V coordinates
    t0 = sprintf('fmincon OPF (%s,%s) : ', bal, crd);

    if ~have_fcn('fmincon')
        t_skip(num_tests, 'fmincon not available');
        continue;
    end

    mpopt = mpoption(mpopt, 'opf.current_balance',  options{k}{1}, ...
                            'opf.v_cartesian',      options{k}{2} );

    %% set up indices
    ib_data     = [1:BUS_AREA BASE_KV:VMIN];
    ib_voltage  = [VM VA];
    ib_lam      = [LAM_P LAM_Q];
    ib_mu       = [MU_VMAX MU_VMIN];
    ig_data     = [GEN_BUS QMAX QMIN MBASE:APF];
    ig_disp     = [PG QG VG];
    ig_mu       = (MU_PMAX:MU_QMIN);
    ibr_data    = (1:ANGMAX);
    ibr_flow    = (PF:QT);
    ibr_mu      = [MU_SF MU_ST];
    ibr_angmu   = [MU_ANGMIN MU_ANGMAX];

    %% get solved AC OPF case from MAT-file
    load soln9_opf;     %% defines bus_soln, gen_soln, branch_soln, f_soln

    %% run OPF
    for s = 0:3
        mpopt = mpoption(mpopt, 'opf.start', s);
        t = sprintf('%s(start=%d): ', t0, s);
        [baseMVA, bus, gen, gencost, branch, f, success, et] = runopf(casefile, mpopt);
        t_ok(success, [t 'success']);
        t_is(f, f_soln, 3, [t 'f']);
        t_is(   bus(:,ib_data   ),    bus_soln(:,ib_data   ), 10, [t 'bus data']);
        t_is(   bus(:,ib_voltage),    bus_soln(:,ib_voltage),  3, [t 'bus voltage']);
        t_is(   bus(:,ib_lam    ),    bus_soln(:,ib_lam    ),  3, [t 'bus lambda']);
        t_is(   bus(:,ib_mu     ),    bus_soln(:,ib_mu     ),  2, [t 'bus mu']);
        t_is(   gen(:,ig_data   ),    gen_soln(:,ig_data   ), 10, [t 'gen data']);
        t_is(   gen(:,ig_disp   ),    gen_soln(:,ig_disp   ),  3, [t 'gen dispatch']);
        t_is(   gen(:,ig_mu     ),    gen_soln(:,ig_mu     ),  3, [t 'gen mu']);
        t_is(branch(:,ibr_data  ), branch_soln(:,ibr_data  ), 10, [t 'branch data']);
        t_is(branch(:,ibr_flow  ), branch_soln(:,ibr_flow  ),  3, [t 'branch flow']);
        t_is(branch(:,ibr_mu    ), branch_soln(:,ibr_mu    ),  2, [t 'branch mu']);
    end
    mpopt = mpoption(mpopt, 'opf.start', 0);    %% set 'opf.start' back to default

    %% run with automatic conversion of single-block pwl to linear costs
    t = [t0 '(single-block PWL) : '];
    mpc = loadcase(casefile);
    mpc.gencost(3, NCOST) = 2;
    [r, success] = runopf(mpc, mpopt);
    [f, bus, gen, branch] = deal(r.f, r.bus, r.gen, r.branch);
    t_ok(success, [t 'success']);
    t_is(f, f_soln, 3, [t 'f']);
    t_is(   bus(:,ib_data   ),    bus_soln(:,ib_data   ), 10, [t 'bus data']);
    t_is(   bus(:,ib_voltage),    bus_soln(:,ib_voltage),  3, [t 'bus voltage']);
    t_is(   bus(:,ib_lam    ),    bus_soln(:,ib_lam    ),  3, [t 'bus lambda']);
    t_is(   bus(:,ib_mu     ),    bus_soln(:,ib_mu     ),  2, [t 'bus mu']);
    t_is(   gen(:,ig_data   ),    gen_soln(:,ig_data   ), 10, [t 'gen data']);
    t_is(   gen(:,ig_disp   ),    gen_soln(:,ig_disp   ),  3, [t 'gen dispatch']);
    t_is(   gen(:,ig_mu     ),    gen_soln(:,ig_mu     ),  3, [t 'gen mu']);
    t_is(branch(:,ibr_data  ), branch_soln(:,ibr_data  ), 10, [t 'branch data']);
    t_is(branch(:,ibr_flow  ), branch_soln(:,ibr_flow  ),  3, [t 'branch flow']);
    t_is(branch(:,ibr_mu    ), branch_soln(:,ibr_mu    ),  2, [t 'branch mu']);
    if mpopt.opf.v_cartesian
        xr = [r.var.val.Vr;r.var.val.Vi;r.var.val.Pg;r.var.val.Qg;0;r.var.val.y];
    else
        xr = [r.var.val.Va;r.var.val.Vm;r.var.val.Pg;r.var.val.Qg;0;r.var.val.y];
    end
    t_is(r.x, xr, 8, [t 'raw x returned from OPF']);

    %% get solved AC OPF case from MAT-file
    load soln9_opf_Plim;       %% defines bus_soln, gen_soln, branch_soln, f_soln

    %% run OPF with active power line limits
    t = [t0 '(P line lim) : '];
    mpopt1 = mpoption(mpopt, 'opf.flow_lim', 'P');
    [baseMVA, bus, gen, gencost, branch, f, success, et] = runopf(casefile, mpopt1);
    t_ok(success, [t 'success']);
    t_is(f, f_soln, 3, [t 'f']);
    t_is(   bus(:,ib_data   ),    bus_soln(:,ib_data   ), 10, [t 'bus data']);
    t_is(   bus(:,ib_voltage),    bus_soln(:,ib_voltage),  3, [t 'bus voltage']);
    t_is(   bus(:,ib_lam    ),    bus_soln(:,ib_lam    ),  3, [t 'bus lambda']);
    t_is(   bus(:,ib_mu     ),    bus_soln(:,ib_mu     ),  2, [t 'bus mu']);
    t_is(   gen(:,ig_data   ),    gen_soln(:,ig_data   ), 10, [t 'gen data']);
    t_is(   gen(:,ig_disp   ),    gen_soln(:,ig_disp   ),  3, [t 'gen dispatch']);
    t_is(   gen(:,ig_mu     ),    gen_soln(:,ig_mu     ),  3, [t 'gen mu']);
    t_is(branch(:,ibr_data  ), branch_soln(:,ibr_data  ), 10, [t 'branch data']);
    t_is(branch(:,ibr_flow  ), branch_soln(:,ibr_flow  ),  3, [t 'branch flow']);
    t_is(branch(:,ibr_mu    ), branch_soln(:,ibr_mu    ),  2, [t 'branch mu']);

    t = [t0 '(P^2 line lim) : '];
    mpopt1 = mpoption(mpopt, 'opf.flow_lim', '2');
    [baseMVA, bus, gen, gencost, branch, f, success, et] = runopf(casefile, mpopt1);
    t_ok(success, [t 'success']);
    t_is(f, f_soln, 3, [t 'f']);
    t_is(   bus(:,ib_data   ),    bus_soln(:,ib_data   ), 10, [t 'bus data']);
    t_is(   bus(:,ib_voltage),    bus_soln(:,ib_voltage),  3, [t 'bus voltage']);
    t_is(   bus(:,ib_lam    ),    bus_soln(:,ib_lam    ),  3, [t 'bus lambda']);
    t_is(   bus(:,ib_mu     ),    bus_soln(:,ib_mu     ),  2, [t 'bus mu']);
    t_is(   gen(:,ig_data   ),    gen_soln(:,ig_data   ), 10, [t 'gen data']);
    t_is(   gen(:,ig_disp   ),    gen_soln(:,ig_disp   ),  3, [t 'gen dispatch']);
    t_is(   gen(:,ig_mu     ),    gen_soln(:,ig_mu     ),  3, [t 'gen mu']);
    t_is(branch(:,ibr_data  ), branch_soln(:,ibr_data  ), 10, [t 'branch data']);
    t_is(branch(:,ibr_flow  ), branch_soln(:,ibr_flow  ),  3, [t 'branch flow']);
    t_is(branch(:,ibr_mu    ), branch_soln(:,ibr_mu    ),  2, [t 'branch mu']);

    %%-----  test OPF with quadratic gen costs moved to generalized costs  -----
    mpc = loadcase(casefile);
    mpc.gencost = [
        2   1500    0   3   0.11    5   0;
        2   2000    0   3   0.085   1.2 0;
        2   3000    0   3   0.1225  1   0;
    ];
    [baseMVA, bus_soln, gen_soln, gencost, branch_soln, f_soln, success, et] = runopf(mpc, mpopt);
    branch_soln = branch_soln(:,1:MU_ST);

    A = sparse(0,0);
    l = [];
    u = [];
    nb = size(mpc.bus, 1);      % number of buses
    ng = size(mpc.gen, 1);      % number of gens
    thbas = 1;                thend    = thbas+nb-1;
    vbas     = thend+1;       vend     = vbas+nb-1;
    pgbas    = vend+1;        pgend    = pgbas+ng-1;
    qgbas    = pgend+1;       qgend    = qgbas+ng-1;
    nxyz = 2*nb + 2*ng;
    N = sparse((1:ng)', (pgbas:pgend)', mpc.baseMVA * ones(ng,1), ng, nxyz);
    fparm = [ 1    0   0 1;
              1 -100 100 1;
              1  -10  10 1 ];
    [junk, ix] = sort(mpc.gen(:, 1));
    H = 2 * spdiags(mpc.gencost(ix, 5), 0, ng, ng);
    Cw = mpc.gencost(ix, 6);
    mpc.gencost(:, 5:7) = 0;

    %% run OPF with quadratic gen costs moved to generalized costs
    t = [t0 'w/quadratic generalized gen cost : '];
    [r, success] = opf(mpc, A, l, u, mpopt, N, fparm, H, Cw);
    [f, bus, gen, branch] = deal(r.f, r.bus, r.gen, r.branch);
    t_ok(success, [t 'success']);
    t_is(f, f_soln, 3, [t 'f']);
    t_is(   bus(:,ib_data   ),    bus_soln(:,ib_data   ), 10, [t 'bus data']);
    t_is(   bus(:,ib_voltage),    bus_soln(:,ib_voltage),  3, [t 'bus voltage']);
    t_is(   bus(:,ib_lam    ),    bus_soln(:,ib_lam    ),  3, [t 'bus lambda']);
    t_is(   bus(:,ib_mu     ),    bus_soln(:,ib_mu     ),  2, [t 'bus mu']);
    t_is(   gen(:,ig_data   ),    gen_soln(:,ig_data   ), 10, [t 'gen data']);
    t_is(   gen(:,ig_disp   ),    gen_soln(:,ig_disp   ),  3, [t 'gen dispatch']);
    t_is(   gen(:,ig_mu     ),    gen_soln(:,ig_mu     ),  3, [t 'gen mu']);
    t_is(branch(:,ibr_data  ), branch_soln(:,ibr_data  ), 10, [t 'branch data']);
    t_is(branch(:,ibr_flow  ), branch_soln(:,ibr_flow  ),  3, [t 'branch flow']);
    t_is(branch(:,ibr_mu    ), branch_soln(:,ibr_mu    ),  2, [t 'branch mu']);
    t_is(r.cost.usr, f, 12, [t 'user cost']);

    %%-----  run OPF with legacy costs and deadzone  -----
    if mpopt.opf.v_cartesian
        t_skip(17, 'legacy cost example n/a to cartesian V case')
    else
        load soln9_opf;
        mpc = loadcase(casefile);
        mpc.N = sparse((1:nb)', (vbas:vend)', ones(nb,1), nb, nxyz);
        mpc.fparm = ones(nb,1) * [ 2 1.08 0.02 1e8 ];
        mpc.Cw = ones(nb, 1);
        t = [t0 'w/legacy cost, in deadzone : '];
        r = runopf(mpc, mpopt);
        [f, bus, gen, branch] = deal(r.f, r.bus, r.gen, r.branch);
        t_ok(r.success, [t 'success']);
        t_is(f, f_soln, 3, [t 'f']);
        t_is(   bus(:,ib_data   ),    bus_soln(:,ib_data   ), 10, [t 'bus data']);
        t_is(   bus(:,ib_voltage),    bus_soln(:,ib_voltage),  3, [t 'bus voltage']);
        t_is(   bus(:,ib_lam    ),    bus_soln(:,ib_lam    ),  3, [t 'bus lambda']);
        t_is(   bus(:,ib_mu     ),    bus_soln(:,ib_mu     ),  2, [t 'bus mu']);
        t_is(   gen(:,ig_data   ),    gen_soln(:,ig_data   ), 10, [t 'gen data']);
        t_is(   gen(:,ig_disp   ),    gen_soln(:,ig_disp   ),  3, [t 'gen dispatch']);
        t_is(   gen(:,ig_mu     ),    gen_soln(:,ig_mu     ),  3, [t 'gen mu']);
        t_is(branch(:,ibr_data  ), branch_soln(:,ibr_data  ), 10, [t 'branch data']);
        t_is(branch(:,ibr_flow  ), branch_soln(:,ibr_flow  ),  3, [t 'branch flow']);
        t_is(branch(:,ibr_mu    ), branch_soln(:,ibr_mu    ),  2, [t 'branch mu']);
        t_is(r.cost.usr, 0, 12, [t 'user cost']);

        t = [t0 'w/legacy cost, not in deadzone : '];
        mpc.fparm = ones(nb,1) * [ 2 1.08 0.01 1e8 ];
        r = runopf(mpc, mpopt);
        [f, bus, gen, branch] = deal(r.f, r.bus, r.gen, r.branch);
        t_ok(r.success, [t 'success']);
        t_is(f, 9009.0890, 3, [t 'f']);
        t_is([min(bus(:, VM)) mean(bus(:, VM)) max(bus(:, VM))], ...
            [1.066624, 1.083980, 1.091698], 5, [t 'bus voltage']);
        t_is(r.cost.usr, 1673.065465, 4, [t 'user cost']);
    end

    %%-----  run OPF with extra linear user constraints & costs  -----
    %% single new z variable constrained to be greater than or equal to
    %% deviation from 1 pu voltage at bus 1, linear cost on this z
    %% get solved AC OPF case from MAT-file
    if mpopt.opf.v_cartesian
        t_skip(14, 'lin constraint/cost example n/a to cartesian V case')
    else
        load soln9_opf_extras1;   %% defines bus_soln, gen_soln, branch_soln, f_soln
        A = sparse([1;1;2;2],[10;25;10;25],[-1;1;1;1],2,25);
        u = [Inf; Inf];
        l = [-1; 1];

        N = sparse(1, 25, 1, 1, 25);    %% new z variable only
        fparm = [1 0 0 1];              %% w = r = z
        H = sparse(1,1);                %% no quadratic term
        Cw = 100;

        t = [t0 'w/extra constraints & costs 1 : '];
        [r, success] = opf(casefile, A, l, u, mpopt, N, fparm, H, Cw);
        [f, bus, gen, branch] = deal(r.f, r.bus, r.gen, r.branch);
        t_ok(success, [t 'success']);
        t_is(f, f_soln, 3, [t 'f']);
        t_is(   bus(:,ib_data   ),    bus_soln(:,ib_data   ), 10, [t 'bus data']);
        t_is(   bus(:,ib_voltage),    bus_soln(:,ib_voltage),  3, [t 'bus voltage']);
        t_is(   bus(:,ib_lam    ),    bus_soln(:,ib_lam    ),  2, [t 'bus lambda']);
        t_is(   bus(:,ib_mu     ),    bus_soln(:,ib_mu     ),  2, [t 'bus mu']);
        t_is(   gen(:,ig_data   ),    gen_soln(:,ig_data   ), 10, [t 'gen data']);
        t_is(   gen(:,ig_disp   ),    gen_soln(:,ig_disp   ),  3, [t 'gen dispatch']);
        t_is(   gen(:,ig_mu     ),    gen_soln(:,ig_mu     ),  3, [t 'gen mu']);
        t_is(branch(:,ibr_data  ), branch_soln(:,ibr_data  ), 10, [t 'branch data']);
        t_is(branch(:,ibr_flow  ), branch_soln(:,ibr_flow  ),  3, [t 'branch flow']);
        t_is(branch(:,ibr_mu    ), branch_soln(:,ibr_mu    ),  2, [t 'branch mu']);
        t_is(r.var.val.z, 0.025419, 6, [t 'user variable']);
        t_is(r.cost.usr, 2.5419, 4, [t 'user cost']);
    end

    %%-----  test OPF with capability curves  -----
    mpc = loadcase('t_case9_opfv2');
    %% remove angle diff limits
    mpc.branch(1, ANGMAX) = 360;
    mpc.branch(9, ANGMIN) = -360;

    %% get solved AC OPF case from MAT-file
    load soln9_opf_PQcap;   %% defines bus_soln, gen_soln, branch_soln, f_soln
    
    %% run OPF with capability curves
    t = [t0 'w/capability curves : '];
    [baseMVA, bus, gen, gencost, branch, f, success, et] = runopf(mpc, mpopt);
    t_ok(success, [t 'success']);
    t_is(f, f_soln, 3, [t 'f']);
    t_is(   bus(:,ib_data   ),    bus_soln(:,ib_data   ), 10, [t 'bus data']);
    t_is(   bus(:,ib_voltage),    bus_soln(:,ib_voltage),  3, [t 'bus voltage']);
    t_is(   bus(:,ib_lam    ),    bus_soln(:,ib_lam    ),  3, [t 'bus lambda']);
    t_is(   bus(:,ib_mu     ),    bus_soln(:,ib_mu     ),  2, [t 'bus mu']);
    t_is(   gen(:,ig_data   ),    gen_soln(:,ig_data   ), 10, [t 'gen data']);
    t_is(   gen(:,ig_disp   ),    gen_soln(:,ig_disp   ),  3, [t 'gen dispatch']);
    t_is(   gen(:,ig_mu     ),    gen_soln(:,ig_mu     ),  3, [t 'gen mu']);
    t_is(branch(:,ibr_data  ), branch_soln(:,ibr_data  ), 10, [t 'branch data']);
    t_is(branch(:,ibr_flow  ), branch_soln(:,ibr_flow  ),  3, [t 'branch flow']);
    t_is(branch(:,ibr_mu    ), branch_soln(:,ibr_mu    ),  2, [t 'branch mu']);

    %%-----  test OPF with angle difference limits  -----
    mpc = loadcase('t_case9_opfv2');
    %% remove capability curves
    mpc.gen(2:3, [PC1, PC2, QC1MIN, QC1MAX, QC2MIN, QC2MAX]) = zeros(2,6);

    %% get solved AC OPF case from MAT-file
    load soln9_opf_ang;   %% defines bus_soln, gen_soln, branch_soln, f_soln
    
    %% run OPF with angle difference limits
    t = [t0 'w/angle difference limits : '];
    [baseMVA, bus, gen, gencost, branch, f, success, et] = runopf(mpc, mpopt);
    t_ok(success, [t 'success']);
    t_is(f, f_soln, 3, [t 'f']);
    t_is(   bus(:,ib_data   ),    bus_soln(:,ib_data   ), 10, [t 'bus data']);
    t_is(   bus(:,ib_voltage),    bus_soln(:,ib_voltage),  3, [t 'bus voltage']);
    t_is(   bus(:,ib_lam    ),    bus_soln(:,ib_lam    ),  3, [t 'bus lambda']);
    t_is(   bus(:,ib_mu     ),    bus_soln(:,ib_mu     ),  1, [t 'bus mu']);
    t_is(   gen(:,ig_data   ),    gen_soln(:,ig_data   ), 10, [t 'gen data']);
    t_is(   gen(:,ig_disp   ),    gen_soln(:,ig_disp   ),  3, [t 'gen dispatch']);
    t_is(   gen(:,ig_mu     ),    gen_soln(:,ig_mu     ),  3, [t 'gen mu']);
    t_is(branch(:,ibr_data  ), branch_soln(:,ibr_data  ), 10, [t 'branch data']);
    t_is(branch(:,ibr_flow  ), branch_soln(:,ibr_flow  ),  3, [t 'branch flow']);
    t_is(branch(:,ibr_mu    ), branch_soln(:,ibr_mu    ),  2, [t 'branch mu']);
    t_is(branch(:,ibr_angmu ), branch_soln(:,ibr_angmu ),  2, [t 'branch angle mu']);

    %%-----  test OPF with ignored angle difference limits  -----
    %% get solved AC OPF case from MAT-file
    load soln9_opf;   %% defines bus_soln, gen_soln, branch_soln, f_soln

    %% run OPF with ignored angle difference limits
    t = [t0 'w/ignored angle difference limits : '];
    mpopt1 = mpoption(mpopt, 'opf.ignore_angle_lim', 1);
    [baseMVA, bus, gen, gencost, branch, f, success, et] = runopf(mpc, mpopt1);
    %% ang limits are not in this solution data, so let's remove them
    branch(1, ANGMAX) = 360;
    branch(9, ANGMIN) = -360;
    t_ok(success, [t 'success']);
    t_is(f, f_soln, 3, [t 'f']);
    t_is(   bus(:,ib_data   ),    bus_soln(:,ib_data   ), 10, [t 'bus data']);
    t_is(   bus(:,ib_voltage),    bus_soln(:,ib_voltage),  3, [t 'bus voltage']);
    t_is(   bus(:,ib_lam    ),    bus_soln(:,ib_lam    ),  3, [t 'bus lambda']);
    t_is(   bus(:,ib_mu     ),    bus_soln(:,ib_mu     ),  2, [t 'bus mu']);
    t_is(   gen(:,ig_data   ),    gen_soln(:,ig_data   ), 10, [t 'gen data']);
    t_is(   gen(:,ig_disp   ),    gen_soln(:,ig_disp   ),  3, [t 'gen dispatch']);
    t_is(   gen(:,ig_mu     ),    gen_soln(:,ig_mu     ),  3, [t 'gen mu']);
    t_is(branch(:,ibr_data  ), branch_soln(:,ibr_data  ), 10, [t 'branch data']);
    t_is(branch(:,ibr_flow  ), branch_soln(:,ibr_flow  ),  3, [t 'branch flow']);
    t_is(branch(:,ibr_mu    ), branch_soln(:,ibr_mu    ),  2, [t 'branch mu']);

    %% angle bounded above by 0, unbounded below
    %% for issue/18
    t = [t0 'w/angle difference limit = 0 : '];
    mpc = loadcase(casefile);
    b = 5;
    mpc.branch(b, ANGMAX) = 0;
    r = runopf(mpc, mpopt);
    t_ok(success, [t 'success']);
    diff = r.bus(r.branch(b, F_BUS), VA) - r.bus(r.branch(b, T_BUS), VA);
    t_is(diff, 0, 5, [t 'angle diff']);

    %%-----  OPF with ref bus not = bus 1, ref angle not = 0  -----
    t = [t0 'ref bus ~= 1, ref ang ~= 0 : '];
    mpc = loadcase(casefile);
    mpc.bus([1;3], BUS_TYPE) = [PV; REF];   %% swap bus types
    bus_soln([1;3], BUS_TYPE) = bus_soln([3;1], BUS_TYPE);   %% swap bus types
    mpc.bus(3, VA) = 3.3014277;
    r = runopf(mpc, mpopt);
    [success, f, bus, gen, branch] = deal(r.success, r.f, r.bus, r.gen, r.branch);
    t_ok(success, [t 'success']);
    t_is(f, f_soln, 3, [t 'f']);
    t_is(   bus(:,ib_data   ),    bus_soln(:,ib_data   ), 10, [t 'bus data']);
    t_is(   bus(:,ib_voltage),    bus_soln(:,ib_voltage),  3, [t 'bus voltage']);
    t_is(   bus(:,ib_lam    ),    bus_soln(:,ib_lam    ),  3, [t 'bus lambda']);
    t_is(   bus(:,ib_mu     ),    bus_soln(:,ib_mu     ),  2, [t 'bus mu']);
    t_is(   gen(:,ig_data   ),    gen_soln(:,ig_data   ), 10, [t 'gen data']);
    t_is(   gen(:,ig_disp   ),    gen_soln(:,ig_disp   ),  3, [t 'gen dispatch']);
    t_is(   gen(:,ig_mu     ),    gen_soln(:,ig_mu     ),  3, [t 'gen mu']);
    t_is(branch(:,ibr_data  ), branch_soln(:,ibr_data  ), 10, [t 'branch data']);
    t_is(branch(:,ibr_flow  ), branch_soln(:,ibr_flow  ),  3, [t 'branch flow']);
    t_is(branch(:,ibr_mu    ), branch_soln(:,ibr_mu    ),  2, [t 'branch mu']);

    %%-----  test OPF with opf.use_vg  -----
    %% get solved AC OPF case from MAT-file
    load soln9_opf_vg;  %% defines bus_soln, gen_soln, branch_soln, f_soln
                        %%    and bus_soln1, gen_soln1, branch_soln1, f_soln1
    
    %% run with opf.use_vg = 1
    t = [t0 'w/opf.use_vg = 1 : '];
    mpc = loadcase(casefile);
    mpc.gen = mpc.gen([1 2 1 3], :);
    mpc.gencost = mpc.gencost([1 2 1 3], :);
    mpc.gen([1 3], [PMAX PMIN]) = mpc.gen([1 3], [PMAX PMIN]) / 2;
    mpc.gen(3, [QMIN, QMAX]) = 0;   %% no reactive capability for gen 3
    mpc.gencost([1 3], COST:end) = mpc.gencost([1 3], COST:end) / 2;
    mpc.gen(1, VG) = 1.05;
    mpc.gen(3, VG) = 1.06;
    mpopt1 = mpoption(mpopt, 'opf.use_vg', 1);
    r = runopf(mpc, mpopt1);
    t_ok(r.success, [t 'success']);
    t_is(r.f, f_soln, 3, [t 'f']);
    t_is(   r.bus(:,ib_data   ),    bus_soln(:,ib_data   ), 10, [t 'bus data']);
    t_is(   r.bus(:,ib_voltage),    bus_soln(:,ib_voltage),  3, [t 'bus voltage']);
    t_is(   r.bus(:,ib_lam    ),    bus_soln(:,ib_lam    ),  3, [t 'bus lambda']);
    t_is(   r.bus(:,ib_mu     ),    bus_soln(:,ib_mu     ),  2, [t 'bus mu']);
    t_is(   r.gen(:,ig_data   ),    gen_soln(:,ig_data   ), 10, [t 'gen data']);
    t_is(   r.gen(:,ig_disp   ),    gen_soln(:,ig_disp   ),  3, [t 'gen dispatch']);
    t_is(   r.gen(:,ig_mu     ),    gen_soln(:,ig_mu     ),  3, [t 'gen mu']);
    t_is(r.branch(:,ibr_data  ), branch_soln(:,ibr_data  ), 10, [t 'branch data']);
    t_is(r.branch(:,ibr_flow  ), branch_soln(:,ibr_flow  ),  3, [t 'branch flow']);
    t_is(r.branch(:,ibr_mu    ), branch_soln(:,ibr_mu    ),  2, [t 'branch mu']);

    %% run with opf.use_vg = 0.9
    t = [t0 'w/opf.use_vg = 0.9 : '];
    mpopt1 = mpoption(mpopt, 'opf.use_vg', 0.9);
    r = runopf(mpc, mpopt1);
    t_ok(r.success, [t 'success']);
    t_is(r.f, f_soln1, 3, [t 'f']);
    t_is(   r.bus(:,ib_data   ),    bus_soln1(:,ib_data   ), 10, [t 'bus data']);
    t_is(   r.bus(:,ib_voltage),    bus_soln1(:,ib_voltage),  3, [t 'bus voltage']);
    t_is(   r.bus(:,ib_lam    ),    bus_soln1(:,ib_lam    ),  3, [t 'bus lambda']);
    t_is(   r.bus(:,ib_mu     ),    bus_soln1(:,ib_mu     ),  2, [t 'bus mu']);
    t_is(   r.gen(:,ig_data   ),    gen_soln1(:,ig_data   ), 10, [t 'gen data']);
    t_is(   r.gen(:,ig_disp   ),    gen_soln1(:,ig_disp   ),  3, [t 'gen dispatch']);
    t_is(   r.gen(:,ig_mu     ),    gen_soln1(:,ig_mu     ),  3, [t 'gen mu']);
    t_is(r.branch(:,ibr_data  ), branch_soln1(:,ibr_data  ), 10, [t 'branch data']);
    t_is(r.branch(:,ibr_flow  ), branch_soln1(:,ibr_flow  ),  3, [t 'branch flow']);
    t_is(r.branch(:,ibr_mu    ), branch_soln1(:,ibr_mu    ),  2, [t 'branch mu']);

    t = [t0 'hi-deg polynomial costs : '];
    mpc = loadcase(casefile);
    mpc.gencost = [
        2   1500    0   6   1e-6/5  0   0   0   0   0;
        2   2000    0   3   1/2     0   0   0   0   0;
        2   3000    0   5   1e-4/4  0   0   0   0   0;
    ];
    r = runopf(mpc, mpopt);
    [f, bus, gen, branch] = deal(r.f, r.bus, r.gen, r.branch);
    t_ok(r.success, [t 'success']);
    t_is(f, 11899.4652, 4, [t 'f']);
    t_is(gen(:, PG), [100.703628; 128.679485; 88.719864], 5, [t 'f']);
    t_is([min(bus(:, VM)) mean(bus(:, VM)) max(bus(:, VM))], ...
        [1.059191 1.079404 1.1], 5, [t 'bus voltage']);

    %% OPF with user-defined nonlinear constraints
    t = [t0 'w/nonlin eq constraint : '];
    mpc = loadcase('case30');
    mpc.user_constraints.nle = {
        {'Pg_usr', 1, 'opf_nle_fcn1', 'opf_nle_hess1', {'Pg'}, {}}
    };
    r = runopf(mpc, mpopt);
    t_ok(r.success, [t 'success']);
    t_is(r.gen(1, PG) * r.gen(2, PG) / 100, r.gen(6, PG), 8, t);
    t_is(r.gen(6, PG), 20.751163, 5, t);

    %% OPF with all buses isolated
    t = [t0 'all buses isolated : '];
    mpc.bus(:, BUS_TYPE) = NONE;
    try
        r = runopf(mpc, mpopt);
        t_is(r.success, 0, 12, [t 'success = 0']);
    catch
        t_ok(0, [t 'unexpected fatal error']);
    end

    %% OPF with no branch limits
    t = [t0 'w/no branch limits : '];
    mpc = loadcase(casefile);
    mpc.branch(:, RATE_A) = 0;
    r = runopf(mpc, mpopt);
    t_ok(r.success, [t 'success']);
    t_is(r.f, 5496.128635, 4, [t 'f']);
    t_is(r.gen(:, PG), [90; 10; 220.463932], 5, [t 'Pg']);
    t_is([min(r.bus(:, VM)) mean(r.bus(:, VM)) max(r.bus(:, VM))], ...
        [1.070692 1.090449 1.1], 5, [t 'bus voltage']);
end

t_end;
