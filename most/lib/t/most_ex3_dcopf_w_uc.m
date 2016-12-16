function most_ex3_dcopf_w_uc(quiet)
%MOST_EX3_DCOPF_W_UC  Examples of deterministic DC optimal power flow with UC.

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

%%-----  DC OPF (w/UC)  -----
casefile = 'ex_case3b';
mpc = loadcase(casefile);
xgd_table.colnames = { 'CommitKey' };
xgd_table.data = [
    1;
    1;
    1;
    2;
];
xgd = loadxgendata(xgd_table, mpc);
[iwind, mpc, xgd] = addwind('ex_wind_uc', mpc, xgd);
mpc = scale_load(499, mpc, [], struct('scale', 'QUANTITY'));
mpc.gencost(:, STARTUP) = 0;    % ignore STARTUP and SHUTDOWN
mpc.gencost(:, SHUTDOWN) = 0;   % costs for this example
mpc.reserves.zones = [mpc.reserves.zones 0];
mpc0 = mpc;
xgd0 = xgd;

%% runduopf
mpc = mpc0;
r1 = runuopf(mpc, mpopt);
u1 = r1.gen(:, GEN_STATUS); % commitment status
Pg1 = r1.gen(:, PG);        % active generation
lam1 = r1.bus(:, LAM_P);    % nodal energy price

%% most
mdi = loadmd(mpc, [], xgd);
mdo = most(mdi, mpopt);
r2 = mdo.flow.mpc;
u2 = r2.gen(:, GEN_STATUS); % commitment status
u2 = mdo.UC.CommitSched;    % commitment status
Pg2 = r2.gen(:, PG);        % active generation
lam2 = r2.bus(:, LAM_P);    % nodal energy price
if verbose
    printpf(r2, [], mpopt);
    most_summary(mdo);
end

%% comparison
u = [u1 u2]
Pg = [Pg1 Pg2]
lam = [lam1 lam2]
