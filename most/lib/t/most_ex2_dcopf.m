function most_ex2_dcopf(quiet)
%MOST_EX2_DCOPF  Examples of deterministic DC optimal power flow.

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
casefile = 'ex_case3a';
mpopt = mpoption('verbose', verbose);
mpopt = mpoption(mpopt, 'out.gen', 1);
mpopt = mpoption(mpopt, 'model', 'DC');
mpopt = mpoption(mpopt, 'opf.dc.solver', 'MIPS');
mpopt = mpoption(mpopt, 'most.solver', mpopt.opf.dc.solver);
if ~verbose
    mpopt = mpoption(mpopt, 'out.all', 0);
end

%%-----  DC OPF  -----
%% runopf
mpc = loadcase(casefile);
r1 = rundcopf(mpc, mpopt);
Pg1 = r1.gen(:, PG);        % active generation
lam1 = r1.bus(:, LAM_P);    % nodal energy price

%% most
mpc = loadcase(casefile);
mpopt = mpoption(mpopt, 'most.dc_model', 1);    % use DC network model (default)
mdi = loadmd(mpc);
mdo = most(mdi, mpopt);
ms = most_summary(mdo);
r2 = mdo.flow.mpc;
Pg2 = r2.gen(:, PG);        % active generation
lam2 = r2.bus(:, LAM_P);    % nodal energy price

%% comparison
Pg = [Pg1 Pg2]
lam = [lam1 lam2]

%%-----  DC OPF (w/reserves)  -----
%% runopf_w_res
mpc = loadcase(casefile);
r1 = runopf_w_res(mpc, mpopt);
Pg1 = r1.gen(:, PG);        % active generation
lam1 = r1.bus(:, LAM_P);    % nodal energy price
R1 = r1.reserves.R;         % reserve quantity
prc1 = r1.reserves.prc;     % reserve price

%% most
mpc = loadcase(casefile);
mpopt = mpoption(mpopt, 'most.dc_model', 1);    % use DC network model (default)
mdi = loadmd(mpc);
mdi.FixedReserves = mpc.reserves;   % include fixed zonal reserves
mdo = most(mdi, mpopt);
ms = most_summary(mdo);
r2 = mdo.flow.mpc;
Pg2 = r2.gen(:, PG);        % active generation
lam2 = r2.bus(:, LAM_P);    % nodal energy price
R2 = r2.reserves.R;         % reserve quantity
prc2 = r2.reserves.prc;     % reserve price

%% comparison
Pg = [Pg1 Pg2]
lam = [lam1 lam2]
R = [R1 R2]
prc = [prc1 prc2]
