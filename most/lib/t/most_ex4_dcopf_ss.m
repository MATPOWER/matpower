function most_ex4_dcopf_ss(quiet)
%MOST_EX4_DCOPF_SS  Examples of secure and stochastic DC OPF problems.

%   MOST
%   Copyright (c) 2015-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

%% set up options
define_constants;
verbose = 0;
mpopt = mpoption('verbose', verbose);
mpopt = mpoption(mpopt, 'out.gen', 1);
mpopt = mpoption(mpopt, 'model', 'DC');
mpopt = mpoption(mpopt, 'opf.dc.solver', 'MIPS');
mpopt = mpoption(mpopt, 'most.solver', 'DEFAULT');
if ~verbose
    mpopt = mpoption(mpopt, 'out.all', 0);
end

casefile = 'ex_case3a';
mpc = loadcase(casefile);
xgd = loadxgendata('ex_xgd_res', mpc);

%%-----  Secure DC OPF (w/cont,res,ramp)  -----
%% most
mdi = loadmd(mpc, [], xgd, [], 'ex_contab');
mdo = most(mdi, mpopt);
EPg1 = mdo.results.ExpectedDispatch;     % active generation
Elam1 = mdo.results.GenPrices;           % nodal energy price
if verbose
    most_summary(mdo);
end


%%-----  Stochastic DC OPF (w/wind,res)  -----
[iwind, mpc, xgd] = addwind('ex_wind', mpc, xgd);
transmat = {[0.158655253931457; 0.682689492137086; 0.158655253931457]};
nt = 1;     % number of periods
nj = 3;     % number of scenarios
profiles = getprofiles(uniformwindprofile(nt, nj), iwind);

%% most
mdi = loadmd(mpc, transmat, xgd, [], [], profiles);
mdo = most(mdi, mpopt);
EPg2 = mdo.results.ExpectedDispatch;     % active generation
Elam2 = mdo.results.GenPrices;           % nodal energy price
if verbose
    most_summary(mdo);
end

%%-----  Secure Stochastic DC OPF (w/wind,cont,res,ramp)  -----
%% most
mdi = loadmd(mpc, transmat, xgd, [], 'ex_contab', profiles);
mdo = most(mdi, mpopt);
EPg3 = mdo.results.ExpectedDispatch;     % active generation
Elam3 = mdo.results.GenPrices;           % nodal energy price
if verbose
    most_summary(mdo);
end

%%-----  Secure Stochastic UC / DC OPF (w/wind,cont,res,ramp)  -----
%% most
casefile = 'ex_case3b';
mpc = loadcase(casefile);
xgd = loadxgendata('ex_xgd_uc', mpc);
[iwind, mpc, xgd] = addwind('ex_wind_uc', mpc, xgd);
mpc = scale_load(350, mpc, [], struct('scale', 'QUANTITY'));
mpc.gencost(:, STARTUP) = 0;
mpc.gencost(:, SHUTDOWN) = 0;
mpc.reserves.zones = [mpc.reserves.zones 0];

mdi = loadmd(mpc, transmat, xgd, [], 'ex_contab', profiles);
mdo = most(mdi, mpopt);
u = mdo.UC.CommitSched;                 % commitment status
EPg4 = mdo.results.ExpectedDispatch;     % active generation
Elam4 = mdo.results.GenPrices;           % nodal energy price
if verbose
    most_summary(mdo);
end

EPg = [[EPg1; NaN] EPg2 EPg3 EPg4]
Elam = [[Elam1; NaN] Elam2 Elam3 Elam4]
