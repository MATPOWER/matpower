function t_most_mpopf(quiet)
% t_most_mpopf - Tests of deterministic multiperiod DC OPF problems
%
%   Cases taken from most_ex5_mpopf.

%   MOST
%   Copyright (c) 2015-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

if nargin < 1
    quiet = 0;      %% verbose by default
end

ntests = 42;
t_begin(ntests, quiet);

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
    near_sing_matrix_warn_id = 'Octave:nearly-singular-matrix';
else
    near_sing_matrix_warn_id = 'MATLAB:nearlySingularMatrix';
end
s = warning('query', near_sing_matrix_warn_id);
warning('off', near_sing_matrix_warn_id);

casefile = 'ex_case3b';
solnfile =  't_most_mpopf_soln';
soln = load(solnfile);
mpopt = mpoption('verbose', verbose);
mpopt = mpoption(mpopt, 'out.gen', 1);
mpopt = mpoption(mpopt, 'model', 'DC');
mpopt = mpoption(mpopt, 'most.solver', 'MIPS');
mpopt = mpoption(mpopt, 'most.dc_model', 1);
% mpopt = mpoption(mpopt, 'opf.violation', 1e-6, 'mips.gradtol', 1e-8, ...
%         'mips.comptol', 1e-8, 'mips.costtol', 1e-8);
if have_feature('octave') && have_feature('octave', 'vnum') > 7
    mpopt = mpoption(mpopt, 'mips.linsolver', 'LU');
end
if ~verbose
    mpopt = mpoption(mpopt, 'out.all', 0);
end
% mpopt = mpoption(mpopt, 'out.all', -1);

%% load base case file
mpc = loadcase(casefile);

nb = size(mpc.bus, 1);
nl = size(mpc.branch, 1);
ng = size(mpc.gen, 1);

xgd = loadxgendata('ex_xgd_ramp', mpc);
[iwind, mpc, xgd] = addwind('ex_wind', mpc, xgd);
profiles = getprofiles('ex_wind_profile_d', iwind);
profiles = getprofiles('ex_load_profile', profiles);
nt = size(profiles(1).values, 1);       % number of periods

mpc_full = mpc;
xgd_full = xgd;

%%-----  Multiperiod DC OPF (w/ramp)  -----
t = 'Multiperiod DC OPF w/ramp : ';
mdi = loadmd(mpc, nt, xgd, [], [], profiles);
mdo = most(mdi, mpopt);
% EPg1 = mdo.results.ExpectedDispatch;    % active generation
% Elam1 = mdo.results.GenPrices;          % nodal energy price
t_is(mdo.results.success, 1, 12, [t 'success']);
ms = most_summary(mdo);
ex = soln.wramp;
t_is(ms.f, ex.f, 1.8, [t 'f']);
t_is(ms.Pg, ex.Pg, 3.5, [t 'Pg']);
t_is(ms.Rup, ex.Rup, 2.5, [t 'Rup']);
t_is(ms.Rdn, ex.Rdn, 2.5, [t 'Rdn']);
t_is(ms.Pf, ex.Pf, 4, [t 'Pf']);
t_is(ms.u, ex.u, 8, [t 'u']);
t_is(ms.lamP, ex.lamP, 5, [t 'lamP']);
t_is(ms.muF, ex.muF, 5, [t 'muF']);
ex = soln.rramp;
t_is(mdo.results.GenPrices, ex.GenPrices, 3, [t 'GenPrices']);
t_is(mdo.results.GenTLMP, ex.GenTLMP, 2, [t 'TLMP']);
% rramp = mdo.results;
% wramp = most_summary(mdo);
% keyboard;

%%-----  Multiperiod DC OPF (w/ramp+wear/tear)  -----
t = 'Multiperiod DC OPF w/ramp+wear/tear : ';
xgd.RampWearCostCoeff(1:3) = 1;
mdi = loadmd(mpc, nt, xgd, [], [], profiles);
mdo = most(mdi, mpopt);
% EPg1 = mdo.results.ExpectedDispatch;    % active generation
% Elam1 = mdo.results.GenPrices;          % nodal energy price
t_is(mdo.results.success, 1, 12, [t 'success']);
ms = most_summary(mdo);
ex = soln.wwear;
t_is(ms.f, ex.f, 0.5, [t 'f']);
t_is(ms.Pg, ex.Pg, 2.8, [t 'Pg']);
t_is(ms.Rup, ex.Rup, 2.5, [t 'Rup']);
t_is(ms.Rdn, ex.Rdn, 2.5, [t 'Rdn']);
t_is(ms.Pf, ex.Pf, 3, [t 'Pf']);
t_is(ms.u, ex.u, 8, [t 'u']);
t_is(ms.lamP, ex.lamP, 2.5, [t 'lamP']);
t_is(ms.muF, ex.muF, 2, [t 'muF']);
ex = soln.rwear;
t_is(mdo.results.GenPrices, ex.GenPrices, 2.8, [t 'GenPrices']);
t_is(mdo.results.GenTLMP, ex.GenTLMP, 1, [t 'TLMP']);
% rwear = mdo.results;
% wwear = most_summary(mdo);
% keyboard;

t = 'build model without solving : ';
mpopt = mpoption(mpopt, 'most.build_model', 1, ...
                        'most.solve_model', 0, ...
                        'most.resolve_new_cost', 0);
mdo = most(mdi, mpopt);
t_ok(isfield(mdo, 'om'), [t '''om'' field']);
t_ok(isfield(mdo, 'QP'), [t '''QP'' field']);
t_ok(~isfield(mdo.QP, 'x'), [t 'no ''QP.x'' field']);
t_ok(~isfield(mdo.QP, 'f'), [t 'no ''QP.f'' field']);
t_ok(~isfield(mdo.QP, 'exitflag'), [t 'no ''QP.exitflag'' field']);
t_ok(~isfield(mdo.QP, 'output'), [t 'no ''QP.output'' field']);
t_ok(isfield(mdo.results, 'SetupTime'), [t '''results.SetupTime'' field']);
t_ok(~isfield(mdo.results, 'success'), [t 'no ''results.success'' field']);

t = 'solve previously built model : ';
mpopt = mpoption(mpopt, 'most.build_model', 0, ...
                        'most.solve_model', 1, ...
                        'most.resolve_new_cost', 1);
mdi1 = mdo;
mdo = most(mdi1, mpopt);
t_is(mdo.results.success, 1, 12, [t 'success']);
ms = most_summary(mdo);
ex = soln.wwear;
t_is(ms.f, ex.f, 0.5, [t 'f']);
t_is(ms.Pg, ex.Pg, 2.8, [t 'Pg']);
t_is(ms.Rup, ex.Rup, 2.5, [t 'Rup']);
t_is(ms.Rdn, ex.Rdn, 2.5, [t 'Rdn']);
t_is(ms.Pf, ex.Pf, 3, [t 'Pf']);
t_is(ms.u, ex.u, 8, [t 'u']);
t_is(ms.lamP, ex.lamP, 2.5, [t 'lamP']);
t_is(ms.muF, ex.muF, 2, [t 'muF']);
ex = soln.rwear;
t_is(mdo.results.GenPrices, ex.GenPrices, 2.8, [t 'GenPrices']);
t_is(mdo.results.GenTLMP, ex.GenTLMP, 1, [t 'TLMP']);
% wwear = most_summary(mdo);
% keyboard;

t = 'output model is copy of input model';
mdo.om.add_var('test', 10);
t_is(mdo.om.var.N, mdi1.om.var.N+10, 12, t);

if have_feature('octave')
    warning(s1.state, file_in_path_warn_id);
end
warning(s.state, near_sing_matrix_warn_id);


t_end;

% save t_most_mpopf_soln wramp wwear rramp rwear
