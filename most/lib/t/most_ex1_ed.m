function most_ex1_ed(quiet)
%MOST_EX1_ED  Examples of deterministic economic dispatch.

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

%%-----  economic dispatch (no network)  -----
%% runopf
mpc = loadcase(casefile);
mpc.branch(:, RATE_A) = 0;  % disable line flow limits (mimic no network case)
r1 = rundcopf(mpc, mpopt);
Pg1 = r1.gen(:, PG);        % active generation
lam1 = r1.bus(:, LAM_P);    % nodal energy price

%% most
mpc = loadcase(casefile);
mpopt = mpoption(mpopt, 'most.dc_model', 0);    % use model with no network
mdi = loadmd(mpc);
mdo = most(mdi, mpopt);
ms = most_summary(mdo);
r2 = mdo.flow.mpc;
Pg2 = r2.gen(:, PG);        % active generation
lam2 = r2.bus(:, LAM_P);    % nodal energy price

%% comparison
Pg = [Pg1 Pg2]
lam = [lam1 lam2]

%%-----  economic dispatch (w/reserves)  -----
%% runopf_w_res
mpc = loadcase(casefile);
mpc.branch(:, RATE_A) = 0;  % disable line flow limits (mimic no network case)
r1 = runopf_w_res(mpc, mpopt);
Pg1 = r1.gen(:, PG);        % active generation
lam1 = r1.bus(:, LAM_P);    % nodal energy price
R1 = r1.reserves.R;         % reserve quantity
prc1 = r1.reserves.prc;     % reserve price

%% most
mpc = loadcase(casefile);
mpopt = mpoption(mpopt, 'most.dc_model', 0);    % use model with no network
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
