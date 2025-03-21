function t_most_uc(quiet, create_plots, create_pdfs, savedir)
% t_most_uc - Tests of deteministic unit commitment optimizations
% ::
%
%   T_MOST_UC(QUIET, CREATE_PLOTS, CREATE_PDFS, SAVEDIR)
%   Can generate summary plots and save them as PDFs in a directory of
%   your choice.
%   E.g. t_most_uc(0, 1, 1, '~/Downloads/uc_plots')

%   MOST
%   Copyright (c) 2015-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

if nargin < 4
    savedir = '.';              %% save in current working directory by default
    if nargin < 3
        create_pdfs = 0;        %% do NOT save plots to PDF files
        if nargin < 2
            create_plots = 0;   %% do NOT create summary plots of results
            if nargin < 1
                quiet = 0;      %% verbose by default
            end
        end
    end
end
if create_plots
    if create_pdfs
        fname = 'uc-ex';
    else
        fname = '';
    end
    pp = 0;     %% plot counter
end

solvers = {'CPLEX', 'GLPK', 'GUROBI', 'MOSEK', 'OT'};
fcn = {'cplex', 'glpk', 'gurobi', 'mosek', 'intlinprog'};
% solvers = {'OT'};
% fcn = {'intlinprog'};
% solvers = {'GUROBI'};
% fcn = {'gurobi'};
% solvers = {'MOSEK'};
% fcn = {'mosek'};
ntests = 69;
t_begin(ntests*length(solvers), quiet);

if quiet
    verbose = 0;
else
    verbose = 0;
end
% verbose = 2;

if have_feature('octave')
    if have_feature('octave', 'vnum') >= 4
        file_in_path_warn_id = 'Octave:data-file-in-path';
    else
        file_in_path_warn_id = 'Octave:load-file-in-path';
    end
    s1 = warning('query', file_in_path_warn_id);
    warning('off', file_in_path_warn_id);
end

casefile = 'ex_case3b';
solnfile =  't_most_uc_soln';
soln = load(solnfile);
mpopt = mpoption;
mpopt = mpoption(mpopt, 'out.gen', 1);
mpopt = mpoption(mpopt, 'verbose', verbose);
% mpopt = mpoption(mpopt, 'opf.violation', 1e-6, 'mips.gradtol', 1e-8, ...
%         'mips.comptol', 1e-8, 'mips.costtol', 1e-8);
mpopt = mpoption(mpopt, 'model', 'DC');
mpopt = mpoption(mpopt, 'most.price_stage_warn_tol', 1e1);

%% solver options
if have_feature('cplex')
    %mpopt = mpoption(mpopt, 'cplex.lpmethod', 0);       %% automatic
    %mpopt = mpoption(mpopt, 'cplex.lpmethod', 1);       %% primal simplex
    mpopt = mpoption(mpopt, 'cplex.lpmethod', 2);       %% dual simplex
    %mpopt = mpoption(mpopt, 'cplex.lpmethod', 3);       %% network simplex
    %mpopt = mpoption(mpopt, 'cplex.lpmethod', 4);       %% barrier
    mpopt = mpoption(mpopt, 'cplex.opts.mip.tolerances.mipgap', 0);
    mpopt = mpoption(mpopt, 'cplex.opts.mip.tolerances.absmipgap', 0);
    mpopt = mpoption(mpopt, 'cplex.opts.simplex.tolerances.optimality', 1e-9);
    mpopt = mpoption(mpopt, 'cplex.opts.simplex.tolerances.feasibility', 1e-9);
    mpopt = mpoption(mpopt, 'cplex.opts.emphasis.numerical', 1);
    mpopt = mpoption(mpopt, 'cplex.opts.threads', 2);
end
if have_feature('glpk')
    mpopt = mpoption(mpopt, 'glpk.opts.mipgap', 0);
    mpopt = mpoption(mpopt, 'glpk.opts.tolint', 1e-10);
    mpopt = mpoption(mpopt, 'glpk.opts.tolobj', 1e-10);
end
if have_feature('gurobi')
    %mpopt = mpoption(mpopt, 'gurobi.method', -1);       %% automatic
    %mpopt = mpoption(mpopt, 'gurobi.method', 0);        %% primal simplex
    mpopt = mpoption(mpopt, 'gurobi.method', 1);        %% dual simplex
    %mpopt = mpoption(mpopt, 'gurobi.method', 2);        %% barrier
    mpopt = mpoption(mpopt, 'gurobi.threads', 4);
    mpopt = mpoption(mpopt, 'gurobi.opts.MIPGap', 0);
    mpopt = mpoption(mpopt, 'gurobi.opts.MIPGapAbs', 0);
end
if have_feature('mosek')
    sc = mosek_symbcon;
    %mpopt = mpoption(mpopt, 'mosek.lp_alg', sc.MSK_OPTIMIZER_FREE);            %% default
    %mpopt = mpoption(mpopt, 'mosek.lp_alg', sc.MSK_OPTIMIZER_INTPNT);          %% interior point
    %mpopt = mpoption(mpopt, 'mosek.lp_alg', sc.MSK_OPTIMIZER_PRIMAL_SIMPLEX);  %% primal simplex
    mpopt = mpoption(mpopt, 'mosek.lp_alg', sc.MSK_OPTIMIZER_DUAL_SIMPLEX);     %% dual simplex
    %mpopt = mpoption(mpopt, 'mosek.lp_alg', sc.MSK_OPTIMIZER_FREE_SIMPLEX);    %% automatic simplex
    %mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_MIO_TOL_X', 0);
    mpopt = mpoption(mpopt, 'mosek.opts.MSK_IPAR_MIO_NODE_OPTIMIZER', sc.MSK_OPTIMIZER_DUAL_SIMPLEX);
    mpopt = mpoption(mpopt, 'mosek.opts.MSK_IPAR_MIO_ROOT_OPTIMIZER', sc.MSK_OPTIMIZER_DUAL_SIMPLEX);
    mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_MIO_TOL_ABS_RELAX_INT', 1e-9);
    %mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_MIO_TOL_REL_RELAX_INT', 0);
    mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_MIO_TOL_REL_GAP', 0);
    mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_MIO_TOL_ABS_GAP', 0);
end
if have_feature('intlinprog')
    mpopt = mpoption(mpopt, 'intlinprog.TolCon', 1e-9);
    mpopt = mpoption(mpopt, 'intlinprog.TolGapAbs', 0);
    mpopt = mpoption(mpopt, 'intlinprog.TolGapRel', 0);
    if have_feature('intlinprog', 'vnum') < 24
        %mpopt = mpoption(mpopt, 'linprog.Algorithm', 'interior-point');
        %mpopt = mpoption(mpopt, 'linprog.Algorithm', 'active-set');
        %mpopt = mpoption(mpopt, 'linprog.Algorithm', 'simplex');
        mpopt = mpoption(mpopt, 'linprog.Algorithm', 'dual-simplex');
        %mpopt = mpoption(mpopt, 'intlinprog.RootLPAlgorithm', 'primal-simplex');
        mpopt = mpoption(mpopt, 'intlinprog.RootLPAlgorithm', 'dual-simplex');
        mpopt = mpoption(mpopt, 'intlinprog.TolInteger', 1e-6);
        %% next line is to work around a bug in intlinprog
        % (Technical Support Case #01841662)
        % (except actually in this case it triggers it rather than working
        %  around it, so we comment it out)
        %mpopt = mpoption(mpopt, 'intlinprog.LPPreprocess', 'none');
    elseif have_feature('intlinprog', 'vnum') == 24.001
        mpopt = mpoption(mpopt, 'intlinprog.LPPreprocess', 'none');
        s2 = warning('query', 'optim:intlinprog:IgnoreOptions');
        warning('off', 'optim:intlinprog:IgnoreOptions');
    elseif have_feature('intlinprog', 'vnum') == 24.002
        mpopt = mpoption(mpopt, 'intlinprog.Presolve', 'off');
    end
end
if ~verbose
    mpopt = mpoption(mpopt, 'out.all', 0);
end
% mpopt = mpoption(mpopt, 'out.all', -1);

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
[CT_LABEL, CT_PROB, CT_TABLE, CT_TBUS, CT_TGEN, CT_TBRCH, CT_TAREABUS, ...
    CT_TAREAGEN, CT_TAREABRCH, CT_ROW, CT_COL, CT_CHGTYPE, CT_REP, ...
    CT_REL, CT_ADD, CT_NEWVAL, CT_TLOAD, CT_TAREALOAD, CT_LOAD_ALL_PQ, ...
    CT_LOAD_FIX_PQ, CT_LOAD_DIS_PQ, CT_LOAD_ALL_P, CT_LOAD_FIX_P, ...
    CT_LOAD_DIS_P, CT_TGENCOST, CT_TAREAGENCOST, CT_MODCOST_F, ...
    CT_MODCOST_X] = idx_ct;

%% load base case file
mpc = loadcase(casefile);

nb = size(mpc.bus, 1);
nl = size(mpc.branch, 1);
ng = size(mpc.gen, 1);

xgd = loadxgendata('ex_xgd_uc', mpc);
[iwind, mpc, xgd] = addwind('ex_wind_uc', mpc, xgd);
profiles = getprofiles('ex_wind_profile_d', iwind);
profiles = getprofiles('ex_load_profile', profiles);
nt = size(profiles(1).values, 1);

mpc_full = mpc;
xgd_full = xgd;

mpc.gencost(:, [STARTUP SHUTDOWN]) = 0; % remove startup/shutdown costs
xgd.MinUp(2) = 1;                       % remove min up-time constraint
xgd.PositiveLoadFollowReserveQuantity(3) = 250; % remove ramp reserve
xgd.PositiveLoadFollowReservePrice(3) = 1e-6;   % constraint and costs
xgd.NegativeLoadFollowReservePrice(3) = 1e-6;
mpc0 = mpc;
xgd0 = xgd;

for s = 1:length(solvers)
    if ~have_feature(fcn{s})    %% check if we have the solver
        t_skip(ntests, sprintf('%s not installed', solvers{s}));
    else
        mpopt = mpoption(mpopt, 'opf.dc.solver', solvers{s});
        mpopt = mpoption(mpopt, 'most.solver', mpopt.opf.dc.solver);
        mpopt = mpoption(mpopt, 'most.storage.cyclic', 1);

        t = sprintf('%s : base (econ disp, no network) : ', solvers{s});
        mpc = mpc0;
        xgd = xgd0;
        mpopt = mpoption(mpopt, 'most.dc_model', 0);
        mdi = loadmd(mpc, nt, xgd, [], [], profiles);
        mdo = most(mdi, mpopt);
        ms = most_summary(mdo);
        t_ok(mdo.QP.exitflag > 0, [t 'success']);
        ex = soln.ed;
        t_is(ms.f, ex.f, 8, [t 'f']);
        t_is(ms.Pg, ex.Pg, 8, [t 'Pg']);
        t_is(ms.Rup, ex.Rup, 8, [t 'Rup']);
        t_is(ms.Rdn, ex.Rdn, 8, [t 'Rdn']);
        t_is(ms.Pf, ex.Pf, 8, [t 'Pf']);
        t_is(ms.u, ex.u, 8, [t 'u']);
        t_is(ms.lamP, ex.lamP, 5, [t 'lamP']);
        t_is(ms.muF, ex.muF, 8, [t 'muF']);
        % ed = most_summary(mdo);
        if create_plots
            pp = pp + 1;
            plot_case('Base : No Network', mdo, ms, 500, 100, savedir, pp, fname);
        end
        % keyboard;

        t = sprintf('%s : + DC OPF constraints : ', solvers{s});
        mpc = mpc0;
        % mpc.gen(iwind, PMAX) = 50;
        mpopt = mpoption(mpopt, 'most.dc_model', 1);
        mdi = loadmd(mpc, nt, xgd, [], [], profiles);
        mdo = most(mdi, mpopt);
        ms = most_summary(mdo);
        t_ok(mdo.QP.exitflag > 0, [t 'success']);
        ex = soln.dcopf;
        t_is(ms.f, ex.f, 8, [t 'f']);
        t_is(ms.Pg, ex.Pg, 8, [t 'Pg']);
        t_is(ms.Rup, ex.Rup, 8, [t 'Rup']);
        t_is(ms.Rdn, ex.Rdn, 8, [t 'Rdn']);
        t_is(ms.Pf, ex.Pf, 8, [t 'Pf']);
        t_is(ms.u, ex.u, 8, [t 'u']);
        t_is(ms.lamP, ex.lamP, 5, [t 'lamP']);
        t_is(ms.muF, ex.muF, 8, [t 'muF']);
        % dcopf = most_summary(mdo);
        if create_plots
            pp = pp + 1;
            plot_case('+ DC Network', mdo, ms, 500, 100, savedir, pp, fname);
        end
        % keyboard;
        mdi.Delta_T = 2;
        mdo = most(mdi, mpopt);
        ms = most_summary(mdo);
        t_is(ms.f, 2 * ex.f, 8, [t '(Delta_T = 2) : f']);

        t = sprintf('%s : + startup/shutdown costs : ', solvers{s});
        if mpopt.out.all
            fprintf('Add STARTUP and SHUTDOWN costs\n');
        end
        mpc = mpc_full;
        % mpc.gencost(3, STARTUP)  = 3524.9944997;    %% CPLEX, GLPK
        % mpc.gencost(3, SHUTDOWN) = 3524.9944997;
        % mpc.gencost(3, STARTUP)  = 3524.99499778;    %% Gurobi
        % mpc.gencost(3, SHUTDOWN) = 3524.9949978;
        % mpc.gencost(3, STARTUP)  = 3524.9949988;    %% MOSEK
        % mpc.gencost(3, SHUTDOWN)  = 3524.9949986;    %% MOSEK
        mdi = loadmd(mpc, nt, xgd, [], [], profiles);
        mdo = most(mdi, mpopt);
        ms = most_summary(mdo);
        t_ok(mdo.QP.exitflag > 0, [t 'success']);
        ex = soln.wstart;
        t_is(ms.f, ex.f, 8, [t 'f']);
        t_is(ms.Pg, ex.Pg, 8, [t 'Pg']);
        t_is(ms.Rup, ex.Rup, 8, [t 'Rup']);
        t_is(ms.Rdn, ex.Rdn, 8, [t 'Rdn']);
        t_is(ms.Pf, ex.Pf, 8, [t 'Pf']);
        t_is(ms.u, ex.u, 8, [t 'u']);
        t_is(ms.lamP, ex.lamP, 8, [t 'lamP']);
        t_is(ms.muF, ex.muF, 8, [t 'muF']);
        % wstart = most_summary(mdo);
        if create_plots
            pp = pp + 1;
            plot_case('+ Startup/Shutdown Costs', mdo, ms, 500, 100, savedir, pp, fname);
        end
        % keyboard;

        t = sprintf('%s : + min up/down time constraints : ', solvers{s});
        if mpopt.out.all
            fprintf('Add MinUp time\n');
        end
        xgd.MinUp(2) = 3;
        mdi = loadmd(mpc, nt, xgd, [], [], profiles);
        mdo = most(mdi, mpopt);
        ms = most_summary(mdo);
        t_ok(mdo.QP.exitflag > 0, [t 'success']);
        ex = soln.wminup;
        t_is(ms.f, ex.f, 6, [t 'f']);
        t_is(ms.Pg, ex.Pg, 7, [t 'Pg']);
        t_is(ms.Rup, ex.Rup, 7, [t 'Rup']);
        t_is(ms.Rdn, ex.Rdn, 8, [t 'Rdn']);
        t_is(ms.Pf, ex.Pf, 8, [t 'Pf']);
        t_is(ms.u, ex.u, 8, [t 'u']);
        t_is(ms.lamP, ex.lamP, 8, [t 'lamP']);
        t_is(ms.muF, ex.muF, 8, [t 'muF']);
        % wminup = most_summary(mdo);
        if create_plots
            pp = pp + 1;
            plot_case('+ Min Up/Down Time Constraints', mdo, ms, 500, 100, savedir, pp, fname);
        end
        % keyboard;

        t = sprintf('%s : + ramp constraint/ramp res cost : ', solvers{s});
        if mpopt.out.all
            fprintf('Restrict ramping and add ramp reserve costs\n');
        end
        xgd = xgd_full;
        mdi = loadmd(mpc, nt, xgd, [], [], profiles);
        mdo = most(mdi, mpopt);
        ms = most_summary(mdo);
        t_ok(mdo.QP.exitflag > 0, [t 'success']);
        ex = soln.wramp;
        t_is(ms.f, ex.f, 8, [t 'f']);
        t_is(ms.Pg, ex.Pg, 8, [t 'Pg']);
        t_is(ms.Rup, ex.Rup, 8, [t 'Rup']);
        t_is(ms.Rdn, ex.Rdn, 8, [t 'Rdn']);
        t_is(ms.Pf, ex.Pf, 8, [t 'Pf']);
        t_is(ms.u, ex.u, 8, [t 'u']);
        t_is(ms.lamP, ex.lamP, 8, [t 'lamP']);
        t_is(ms.muF, ex.muF, 8, [t 'muF']);
        ex = soln.rramp;
        t_is(mdo.results.GenPrices, ex.GenPrices, 8, [t 'GenPrices']);
        t_is(mdo.results.GenTLMP, ex.GenTLMP, 5, [t 'TLMP']);
        % rramp = mdo.results;
        % wramp = most_summary(mdo);
        if create_plots
            pp = pp + 1;
            plot_case('+ Ramping Constraints/Ramp Reserve Costs', mdo, ms, 500, 100, savedir, pp, fname);
        end
        % keyboard;

        t = sprintf('%s : + storage : ', solvers{s});
        if mpopt.out.all
            fprintf('Add storage\n');
        end
        [iess, mpc, xgd, sd] = addstorage('ex_storage', mpc, xgd);
        mdi = loadmd(mpc, nt, xgd, sd, [], profiles);
        mdo = most(mdi, mpopt);
        ms = most_summary(mdo);
        t_ok(mdo.QP.exitflag > 0, [t 'success']);
        ex = soln.wstorage;
        t_is(ms.f, ex.f, 8, [t 'f']);
        t_is(ms.Pg, ex.Pg, 8, [t 'Pg']);
        t_is(ms.Rup, ex.Rup, 8, [t 'Rup']);
        t_is(ms.Rdn, ex.Rdn, 8, [t 'Rdn']);
        t_is(ms.Pf, ex.Pf, 8, [t 'Pf']);
        t_is(ms.u, ex.u, 8, [t 'u']);
        % t_is(ms.lamP, ex.lamP, 5, [t 'lamP']);
        % t_is(ms.muF, ex.muF, 5, [t 'muF']);
        % wstorage = most_summary(mdo);
        if create_plots
            pp = pp + 1;
            plot_case('+ Storage', mdo, ms, 500, 100, savedir, pp, fname);
            create_plots = 0;   %% don't do them again
        end
        % keyboard;

        t = sprintf('%s : + storage2 : ', solvers{s});
        if mpopt.out.all
            fprintf('Add storage\n');
        end
        mpopt = mpoption(mpopt, 'most.storage.cyclic', 0);
        mdi.Storage.rho = 1;
        mdi.Storage.InitialStorageUpperBound = 0;
        mdo = most(mdi, mpopt);
        ms = most_summary(mdo);
        t_ok(mdo.QP.exitflag > 0, [t 'success']);
        ex = soln.wstorage3;
        t_is(ms.f, ex.f, 8, [t 'f']);
        t_is(ms.Pg, ex.Pg, 8, [t 'Pg']);
        t_is(ms.Rup, ex.Rup, 8, [t 'Rup']);
        t_is(ms.Rdn, ex.Rdn, 8, [t 'Rdn']);
        t_is(ms.Pf, ex.Pf, 8, [t 'Pf']);
        t_is(ms.u, ex.u, 8, [t 'u']);
        % t_is(ms.lamP, ex.lamP, 5, [t 'lamP']);
        % t_is(ms.muF, ex.muF, 5, [t 'muF']);
        % wstorage2 = most_summary(mdo);
        if create_plots
            pp = pp + 1;
            plot_case('+ Storage2', mdo, ms, 500, 100, savedir, pp, fname);
            create_plots = 0;   %% don't do them again
        end
        % keyboard;

        t = sprintf('%s : + storage3 : ', solvers{s});
        if mpopt.out.all
            fprintf('Add storage\n');
        end
        mdi.Storage.rho = 0;
        mdo = most(mdi, mpopt);
        ms = most_summary(mdo);
        t_ok(mdo.QP.exitflag > 0, [t 'success']);
        ex = soln.wstorage3;
        t_is(ms.f, ex.f, 8, [t 'f']);
        t_is(ms.Pg, ex.Pg, 8, [t 'Pg']);
        t_is(ms.Rup, ex.Rup, 3.5, [t 'Rup']);
        t_is(ms.Rdn, ex.Rdn, 3.5, [t 'Rdn']);
        t_is(ms.Pf, ex.Pf, 8, [t 'Pf']);
        t_is(ms.u, ex.u, 8, [t 'u']);
        % t_is(ms.lamP, ex.lamP, 5, [t 'lamP']);
        % t_is(ms.muF, ex.muF, 5, [t 'muF']);
        % wstorage3 = most_summary(mdo);
        if create_plots
            pp = pp + 1;
            plot_case('+ Storage3', mdo, ms, 500, 100, savedir, pp, fname);
            create_plots = 0;   %% don't do them again
        end
        % keyboard;
    end
end

if have_feature('octave')
    warning(s1.state, file_in_path_warn_id);
elseif have_feature('intlinprog') && have_feature('intlinprog', 'vnum') == 24.001
    warning(s2.state, 'optim:intlinprog:IgnoreOptions');
end

t_end;

% save t_most_uc_soln ed dcopf wstart wminup wramp wstorage wstorage2 wstorage3 rramp

function h = plot_case(label, md, ms, maxq, maxp, mypath, pp, fname)

if nargin < 8
    fname = '';
end

%% colors:  blue     red               yellow           purple            green
cc = {[0 0.45 0.74], [0.85 0.33 0.1], [0.93 0.69 0.13], [0.49 0.18 0.56], [0.47 0.67 0.19]};

ig = (1:3)';
id = 4;
iw = 5;
is = 6;

subplot(3, 1, 1);
md.mpc = rmfield(md.mpc, 'genfuel');
plot_uc(md, [], 'title', label);
ylabel('Unit Commitment', 'FontSize', 16);
ah = gca;
ah.YAxisLocation = 'left';

subplot(3, 1, 2);
x = (1:ms.nt)';
y1 = ms.Pg(ig, :)';
if ms.ng == 6
    y1 = [y1 max(-ms.Pg(is, :), 0)' max(ms.Pg(is, :), 0)'];
end
y2 = -sum(ms.Pg([id; iw], :), 1)';
[ah1, h1, h2] = plotyy(x, y1, x, y2);
axis(ah1(1), [0.5 12.5 0 maxq]);
axis(ah1(2), [0.5 12.5 0 maxq]);
% ah1(1).XLim = [0.5 12.5];
% ah1(2).XLim = [0.5 12.5];
% ah1(1).YLim = [0 300];
% ah1(2).YLim = [0 450];
ah1(1).YTickMode = 'auto';
ah1(2).YTickMode = 'auto';
ah1(1).XTick = 1:12;
nn = 3;
for j = 1:3
    h1(j).LineWidth = 2;
    h1(j).Color = cc{j};
end
if ms.ng == 6
    h1(4).LineWidth = 2;
    h1(4).Color = cc{5};
    h1(4).LineStyle = ':';
    h1(5).LineWidth = 2;
    h1(5).Color = cc{5};
end
h2.LineWidth = 2;
h2.Color = cc{4};
h2.LineStyle = ':';
ah1(2).YColor = cc{4};
%title('Generation & Net Load', 'FontSize', 16);
ylabel(ah1(1), 'Generation, MW', 'FontSize', 16);
ylabel(ah1(2), 'Net Load, MW', 'FontSize', 16);
xlabel('Period', 'FontSize', 16);
set(ah1(1), 'FontSize', 14);
set(ah1(2), 'FontSize', 14);
if ms.ng == 6
    legend('Gen 1', 'Gen 2', 'Gen 3', 'Storage Charge', 'Storage Discharge', 'Location', [0.7 0.6 0 0]);
else
    legend('Gen 1', 'Gen 2', 'Gen 3', 'Location', [0.7 0.58 0 0]);
end

subplot(3, 1, 3);
y1 = ms.lamP';
plot(x, y1, 'LineWidth', 2);
% title('Nodal Price', 'FontSize', 16);
ylabel('Nodal Price, $/MWh', 'FontSize', 16);
xlabel('Period', 'FontSize', 16);
axis([0.5 12.5 0 maxp]);
ah = gca;
set(ah, 'FontSize', 14);
ah.XTick = 1:12;
legend('Bus 1', 'Bus 2', 'Bus 3', 'Location', [0.7 0.28 0 0]);

if nargin > 7 && ~isempty(fname)
    h = gcf;
    set(h, 'PaperSize', [11 8.5]);
    set(h, 'PaperPosition', [0.25 0.25 10.5 8]);
    print('-dpdf', fullfile(mypath, sprintf('%s-%d', fname, pp)));
end
