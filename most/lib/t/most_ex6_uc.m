function most_ex6_uc(quiet)
%MOST_EX6_UC  Examples of deterministic unit commitment problems.

%   MOST
%   Copyright (c) 2015-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

%% set up options
define_constants;
verbose = 1;
mpopt = mpoption('verbose', verbose);
mpopt = mpoption(mpopt, 'out.gen', 1);
mpopt = mpoption(mpopt, 'model', 'DC');
mpopt = mpoption(mpopt, 'most.solver', 'DEFAULT');
if ~verbose
    mpopt = mpoption(mpopt, 'out.all', 0);
end
mpopt = mpoption(mpopt, 'most.price_stage_warn_tol', 1e-5);

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
if have_fcn('intlinprog')
    %mpopt = mpoption(mpopt, 'linprog.Algorithm', 'interior-point');
    %mpopt = mpoption(mpopt, 'linprog.Algorithm', 'active-set');
    %mpopt = mpoption(mpopt, 'linprog.Algorithm', 'simplex');
    mpopt = mpoption(mpopt, 'linprog.Algorithm', 'dual-simplex');
    %mpopt = mpoption(mpopt, 'intlinprog.RootLPAlgorithm', 'primal-simplex');
    mpopt = mpoption(mpopt, 'intlinprog.RootLPAlgorithm', 'dual-simplex');
    mpopt = mpoption(mpopt, 'intlinprog.TolCon', 1e-9);
    mpopt = mpoption(mpopt, 'intlinprog.TolGapAbs', 0);
    mpopt = mpoption(mpopt, 'intlinprog.TolGapRel', 0);
    mpopt = mpoption(mpopt, 'intlinprog.TolInteger', 1e-6);
    %% next line is to work around a bug in intlinprog
    % (Technical Support Case #01841662)
    % (except actually in this case it triggers it rather than working
    %  around it, so we comment it out)
    %mpopt = mpoption(mpopt, 'intlinprog.LPPreprocess', 'none');
end

casefile = 'ex_case3b';
mpc = loadcase(casefile);
xgd = loadxgendata('ex_xgd_uc', mpc);
[iwind, mpc, xgd] = addwind('ex_wind_uc', mpc, xgd);
profiles = getprofiles('ex_wind_profile_d', iwind);
profiles = getprofiles('ex_load_profile', profiles);
nt = size(profiles(1).values, 1);       % number of periods
mpc_full = mpc;
xgd_full = xgd;
mpc.gencost(:, [STARTUP SHUTDOWN]) = 0; % remove startup/shutdown costs
xgd.MinUp(2) = 1;                       % remove min up-time constraint
xgd.PositiveLoadFollowReserveQuantity(3) = 250; % remove ramp reserve
xgd.PositiveLoadFollowReservePrice(3) = 1e-6;   % constraint and costs
xgd.NegativeLoadFollowReservePrice(3) = 1e-6;


%%-----  Base : No Network  -----
mpopt = mpoption(mpopt, 'most.dc_model', 0);    % use model with no network
mdi = loadmd(mpc, nt, xgd, [], [], profiles);
mdo = most(mdi, mpopt);
if verbose
    ms = most_summary(mdo);
end

%%-----  + DC Network  -----
mpopt = mpoption(mpopt, 'most.dc_model', 1);    % use DC network model (default)
mdo = most(mdi, mpopt);
if verbose
    ms = most_summary(mdo);
end

%%-----  + startup/shutdown costs  -----
mpc.gencost(2, [STARTUP SHUTDOWN]) = [ 200 200];
mpc.gencost(3, [STARTUP SHUTDOWN]) = [3000 600];
% equivalent to doing: mpc = mpc_full;
mdi = loadmd(mpc, nt, xgd, [], [], profiles);
mdo = most(mdi, mpopt);
if verbose
    ms = most_summary(mdo);
end

%%-----  + min up/down time constraints  -----
xgd.MinUp(2) = 3;
mdi = loadmd(mpc, nt, xgd, [], [], profiles);
mdo = most(mdi, mpopt);
if verbose
    ms = most_summary(mdo);
end

%%-----  + ramp constraint/ramp res cost  -----
xgd.PositiveLoadFollowReserveQuantity(3) = 100; % restore ramp reserve
xgd.PositiveLoadFollowReservePrice(3) = 10;     % constraint and costs
xgd.NegativeLoadFollowReservePrice(3) = 10;
% equivalent to doing: xgd = xgd_full;
mdi = loadmd(mpc, nt, xgd, [], [], profiles);
mdo = most(mdi, mpopt);
if verbose
    ms = most_summary(mdo);
end

%%-----  + storage  -----
mpopt = mpoption(mpopt, 'most.storage.cyclic', 1);
[iess, mpc, xgd, sd] = addstorage('ex_storage', mpc, xgd);
mdi = loadmd(mpc, nt, xgd, sd, [], profiles);
mdo = most(mdi, mpopt);
if verbose
    ms = most_summary(mdo);
end
