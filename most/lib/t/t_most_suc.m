function t_most_suc(quiet)
%T_MOST_SUC  Tests of stochastic unit commitment optimizations.

%   MOST
%   Copyright (c) 2015-2016 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   $Id$
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

solvers = {'CPLEX', 'GLPK', 'GUROBI', 'MOSEK', 'OT'};
fcn = {'cplex', 'glpk', 'gurobi', 'mosek', 'intlinprog'};
% solvers = {'OT'};
% fcn = {'intlinprog'};
% solvers = {'GUROBI'};
% fcn = {'gurobi'};
% solvers = {'DEFAULT'};
% fcn = {'most'};
ntests = 37;
t_begin(ntests*length(solvers), quiet);

if quiet
    verbose = 0;
else
    verbose = 0;
end
% verbose = 2;

if have_fcn('octave')
    if have_fcn('octave', 'vnum') >= 4
        file_in_path_warn_id = 'Octave:data-file-in-path';
    else
        file_in_path_warn_id = 'Octave:load-file-in-path';
    end
    s1 = warning('query', file_in_path_warn_id);
    warning('off', file_in_path_warn_id);
end

casefile = 'ex_case3b';
solnfile =  't_most_suc_soln';
soln = load(solnfile);
mpopt = mpoption;
mpopt = mpoption(mpopt, 'out.gen', 1);
mpopt = mpoption(mpopt, 'verbose', verbose);
% mpopt = mpoption(mpopt, 'opf.violation', 1e-6, 'mips.gradtol', 1e-8, ...
%         'mips.comptol', 1e-8, 'mips.costtol', 1e-8);
mpopt = mpoption(mpopt, 'model', 'DC');
mpopt = mpoption(mpopt, 'most.price_stage_warn_tol', 10);

%% solver options
if have_fcn('cplex')
    %mpopt = mpoption(mpopt, 'cplex.lpmethod', 0);       %% automatic
    %mpopt = mpoption(mpopt, 'cplex.lpmethod', 1);       %% primal simplex
    mpopt = mpoption(mpopt, 'cplex.lpmethod', 2);       %% dual simplex
    %mpopt = mpoption(mpopt, 'cplex.lpmethod', 3);       %% network simplex
    %mpopt = mpoption(mpopt, 'cplex.lpmethod', 4);       %% barrier
    mpopt = mpoption(mpopt, 'cplex.opts.mip.tolerances.mipgap', 0);
    mpopt = mpoption(mpopt, 'cplex.opts.mip.tolerances.absmipgap', 0);
    mpopt = mpoption(mpopt, 'cplex.opts.threads', 2);
end
if have_fcn('glpk')
    mpopt = mpoption(mpopt, 'glpk.opts.mipgap', 0);
    mpopt = mpoption(mpopt, 'glpk.opts.tolint', 1e-10);
    mpopt = mpoption(mpopt, 'glpk.opts.tolobj', 1e-10);
end
if have_fcn('gurobi')
    %mpopt = mpoption(mpopt, 'gurobi.method', -1);       %% automatic
    %mpopt = mpoption(mpopt, 'gurobi.method', 0);        %% primal simplex
    mpopt = mpoption(mpopt, 'gurobi.method', 1);        %% dual simplex
    %mpopt = mpoption(mpopt, 'gurobi.method', 2);        %% barrier
    mpopt = mpoption(mpopt, 'gurobi.threads', 2);
    mpopt = mpoption(mpopt, 'gurobi.opts.MIPGap', 0);
    mpopt = mpoption(mpopt, 'gurobi.opts.MIPGapAbs', 0);
end
if have_fcn('mosek')
    %mpopt = mpoption(mpopt, 'mosek.lp_alg', 0);         %% automatic
    %mpopt = mpoption(mpopt, 'mosek.lp_alg', 1);         %% interior point
    %mpopt = mpoption(mpopt, 'mosek.lp_alg', 3);         %% primal simplex
    mpopt = mpoption(mpopt, 'mosek.lp_alg', 4);         %% dual simplex
    sc = mosek_symbcon;
    %mpopt = mpoption(mpopt, 'mosek.lp_alg', 5);         %% primal dual simplex
    %mpopt = mpoption(mpopt, 'mosek.lp_alg', 6);         %% automatic simplex
    %mpopt = mpoption(mpopt, 'mosek.lp_alg', 7);         %% network primal simplex
    %mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_MIO_TOL_X', 0);
    mpopt = mpoption(mpopt, 'mosek.opts.MSK_IPAR_MIO_NODE_OPTIMIZER', sc.MSK_OPTIMIZER_DUAL_SIMPLEX);
    mpopt = mpoption(mpopt, 'mosek.opts.MSK_IPAR_MIO_ROOT_OPTIMIZER', sc.MSK_OPTIMIZER_DUAL_SIMPLEX);
    mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_MIO_TOL_ABS_RELAX_INT', 1e-9);
    mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_MIO_TOL_REL_RELAX_INT', 0);
    mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_MIO_TOL_REL_GAP', 0);
    mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_MIO_TOL_ABS_GAP', 0);
end
if have_fcn('intlinprog')
    %mpopt = mpoption(mpopt, 'linprog.Algorithm', 'interior-point');
    %mpopt = mpoption(mpopt, 'linprog.Algorithm', 'active-set');
    %mpopt = mpoption(mpopt, 'linprog.Algorithm', 'simplex');
    mpopt = mpoption(mpopt, 'linprog.Algorithm', 'dual-simplex');
    mpopt = mpoption(mpopt, 'intlinprog.RootLPAlgorithm', 'primal-simplex');
    % mpopt = mpoption(mpopt, 'intlinprog.RootLPAlgorithm', 'dual-simplex');
    mpopt = mpoption(mpopt, 'intlinprog.TolCon', 1e-9);
    mpopt = mpoption(mpopt, 'intlinprog.TolGapAbs', 0);
    mpopt = mpoption(mpopt, 'intlinprog.TolGapRel', 0);
    mpopt = mpoption(mpopt, 'intlinprog.TolInteger', 1e-6);
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


profiles = getprofiles('ex_load_profile');
nt = size(profiles.values, 1);

nb = size(mpc.bus, 1);
nl = size(mpc.branch, 1);
ng = size(mpc.gen, 1);

xgd = loadxgendata('ex_xgd', mpc);
[iwind, mpc, xgd] = addwind('ex_wind', mpc, xgd);
profiles_s = getprofiles('ex_wind_profile', iwind);
profiles_s = getprofiles('ex_load_profile', profiles_s);

mpc0 = mpc;
xgd0 = xgd;

for s = 1:length(solvers)
    if ~have_fcn(fcn{s})     %% check if we have the solver
        t_skip(ntests, sprintf('%s not installed', solvers{s}));
    else
        mpopt = mpoption(mpopt, 'opf.dc.solver', solvers{s});
        mpopt = mpoption(mpopt, 'most.solver', mpopt.opf.dc.solver);
        mpopt = mpoption(mpopt, 'most.storage.cyclic', 1);

        mpc = mpc0;
        xgd = xgd0;

        t = sprintf('%s : deterministic : ', solvers{s});
        mdi = loadmd(mpc, nt, xgd, [], [], profiles);
        mdo = most(mdi, mpopt);
        ms = most_summary(mdo);
        t_ok(mdo.QP.exitflag > 0, [t 'success']);
        ex = soln.determ;
        t_is(ms.f, ex.f, 8, [t 'f']);
        t_is(ms.Pg, ex.Pg, 8, [t 'Pg']);
        t_is(ms.Rup, ex.Rup, 8, [t 'Rup']);
        t_is(ms.Rdn, ex.Rdn, 8, [t 'Rdn']);
        t_is(ms.Pf, ex.Pf, 8, [t 'Pf']);
        t_is(ms.u, ex.u, 8, [t 'u']);
        t_is(ms.lamP, ex.lamP, 8, [t 'lamP']);
        t_is(ms.muF, ex.muF, 8, [t 'muF']);
        % determ = most_summary(mdo);
        % keyboard;

        t = sprintf('%s : individual trajectories : ', solvers{s});
        transmat_s = cell(1, nt);
        I = speye(3);
        [transmat_s{:}] = deal(I);
        transmat_s{1} = [0.2; 0.6; 0.2];
        mdi = loadmd(mpc, transmat_s, xgd, [], [], profiles_s);
        mdi = filter_ramp_transitions(mdi, 0.1);
        mdo = most(mdi, mpopt);
        ms = most_summary(mdo);
        t_ok(mdo.QP.exitflag > 0, [t 'success']);
        ex = soln.transprob1;
        t_is(ms.f, ex.f, 5, [t 'f']);
        t_is(ms.Pg, ex.Pg, 6, [t 'Pg']);
        t_is(ms.Rup, ex.Rup, 8, [t 'Rup']);
        t_is(ms.Rdn, ex.Rdn, 8, [t 'Rdn']);
        t_is(ms.Pf, ex.Pf, 8, [t 'Pf']);
        t_is(ms.u, ex.u, 8, [t 'u']);
        % t_is(ms.lamP, ex.lamP, 5, [t 'lamP']);
        % t_is(ms.muF, ex.muF, 5, [t 'muF']);
        % transprob1 = most_summary(mdo);
        % keyboard;

        t = sprintf('%s : full transition probabilities : ', solvers{s});
%        transmat_sf = transmat_s;
        transmat_sf = ex_transmat(nt);
%         transmat_sf = cell(1, nt);
%         T = [ 0.158655253931457; 0.682689492137086; 0.158655253931457 ];
%         [transmat_sf{:}] = deal(T * ones(1,3));
%         transmat_sf{1} = T;
        mdi = loadmd(mpc, transmat_sf, xgd, [], [], profiles_s);
%        mdi = filter_ramp_transitions(mdi, 0.9);
        mdo = most(mdi, mpopt);
        ms = most_summary(mdo);
        t_ok(mdo.QP.exitflag > 0, [t 'success']);
        ex = soln.transprobfull;
        t_is(ms.f, ex.f, 3, [t 'f']);
        t_is(ms.Pg, ex.Pg, 3, [t 'Pg']);
        t_is(ms.Rup, ex.Rup, 3, [t 'Rup']);
        t_is(ms.Rdn, ex.Rdn, 3, [t 'Rdn']);
        t_is(ms.Pf, ex.Pf, 8, [t 'Pf']);
        t_is(ms.u, ex.u, 8, [t 'u']);
        % t_is(ms.lamP, ex.lamP, 5, [t 'lamP']);
        % t_is(ms.muF, ex.muF, 5, [t 'muF']);
        % transprobfull = most_summary(mdo);
        % keyboard;

        t = sprintf('%s : full transition probabilities + cont : ', solvers{s});
        mdi = loadmd(mpc, transmat_sf, xgd, [], 'ex_contab', profiles_s);
%        mdi = filter_ramp_transitions(mdi, 0.9);
        mdo = most(mdi, mpopt);
        ms = most_summary(mdo);
        t_ok(mdo.QP.exitflag > 0, [t 'success']);
        ex = soln.transprobcont;
        t_is(ms.f, ex.f, 4, [t 'f']);
        t_is(ms.Pg, ex.Pg, 6, [t 'Pg']);
        t_is(ms.Rup, ex.Rup, 6, [t 'Rup']);
        t_is(ms.Rdn, ex.Rdn, 6, [t 'Rdn']);
        t_is(ms.Pf, ex.Pf, 6, [t 'Pf']);
        t_is(ms.u, ex.u, 8, [t 'u']);
        % t_is(ms.lamP, ex.lamP, 5, [t 'lamP']);
        % t_is(ms.muF, ex.muF, 5, [t 'muF']);
        % transprobcont = most_summary(mdo);
        % keyboard;

        t = sprintf('%s : + storage : ', solvers{s});
        if mpopt.out.all
            fprintf('Add storage\n');
        end
        [iess, mpc, xgd, sd] = addstorage('ex_storage', mpc, xgd);
        mdi = loadmd(mpc, transmat_sf, xgd, sd, 'ex_contab', profiles_s);
        mdo = most(mdi, mpopt);
        ms = most_summary(mdo);
        t_ok(mdo.QP.exitflag > 0, [t 'success']);
        ex = soln.wstorage;
        t_is(ms.f, ex.f, 3, [t 'f']);
        t_is(ms.Pg, ex.Pg, 3, [t 'Pg']);
        t_is(ms.Rup, ex.Rup, 8, [t 'Rup']);
        t_is(ms.Rdn, ex.Rdn, 8, [t 'Rdn']);
        t_is(ms.Pf, ex.Pf, 3, [t 'Pf']);
        t_is(ms.u, ex.u, 8, [t 'u']);
        % t_is(ms.lamP, ex.lamP, 5, [t 'lamP']);
        % t_is(ms.muF, ex.muF, 5, [t 'muF']);
        % wstorage = most_summary(mdo);
        % keyboard;
    end
end

if have_fcn('octave')
    warning(s1.state, file_in_path_warn_id);
end

t_end;

% save t_most_suc_soln determ transprob1 transprobfull transprobcont wstorage
