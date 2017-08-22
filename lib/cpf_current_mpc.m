function mpc = cpf_current_mpc(mpc, mpct, Ybus, Yf, Yt, ref, pv, pq, V, lam, mpopt)
%CPF_CURRENT_MPC  Construct MPC for current continuation step.
%   MPC = CPF_CURRENT_MPC(MPC_BASE, MPC_TARGET, YBUS, YF, YT, REF, PV, PQ, V, LAM, MPOPT)
%
%   Constructs the MATPOWER case struct for the current continuation step
%   based on the MPC_BASE and MPC_TARGET cases and the value of LAM.

%   MATPOWER
%   Copyright (c) 2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% update bus and gen matrices to reflect the loading and generation
mpc.bus(:, PD) = mpc.bus(:, PD) + lam * (mpct.bus(:, PD) - mpc.bus(:, PD));
mpc.bus(:, QD) = mpc.bus(:, QD) + lam * (mpct.bus(:, QD) - mpc.bus(:, QD));
mpc.gen(:, PG) = mpc.gen(:, PG) + lam * (mpct.gen(:, PG) - mpc.gen(:, PG));

%% update data matrices with solution
[mpc.bus, mpc.gen, mpc.branch] = pfsoln(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, Ybus, Yf, Yt, V, ref, pv, pq, mpopt);
