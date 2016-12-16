function most_ex5_mpopf(quiet)
%MOST_EX5_MPOPF  Examples of deterministic multiperiod DC OPF problems.

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
mpopt = mpoption(mpopt, 'opf.dc.solver', 'MIPS');
mpopt = mpoption(mpopt, 'most.solver', 'MIPS');
mpopt = mpoption(mpopt, 'most.dc_model', 1);
if ~verbose
    mpopt = mpoption(mpopt, 'out.all', 0);
end

casefile = 'ex_case3b';
mpc = loadcase(casefile);
xgd = loadxgendata('ex_xgd_ramp', mpc);
[iwind, mpc, xgd] = addwind('ex_wind', mpc, xgd);
profiles = getprofiles('ex_wind_profile_d', iwind);
profiles = getprofiles('ex_load_profile', profiles);
nt = size(profiles(1).values, 1);       % number of periods

%%-----  Multiperiod DC OPF (w/ramp)  -----
mdi = loadmd(mpc, nt, xgd, [], [], profiles);
mdo = most(mdi, mpopt);
EPg1 = mdo.results.ExpectedDispatch;    % active generation
Elam1 = mdo.results.GenPrices;          % nodal energy price
if verbose
    most_summary(mdo);
end

%%-----  Multiperiod DC OPF (w/ramp+wear/tear)  -----
xgd.RampWearCostCoeff(1:3) = 1;
mdi = loadmd(mpc, nt, xgd, [], [], profiles);
mdo = most(mdi, mpopt);
EPg1 = mdo.results.ExpectedDispatch;    % active generation
Elam1 = mdo.results.GenPrices;          % nodal energy price
if verbose
    most_summary(mdo);
end
