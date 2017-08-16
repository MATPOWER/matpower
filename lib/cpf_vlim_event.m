function ef = cpf_vlim_event(cb_data, cx)
%CPF_VLIM_EVENT  Event function to detect a bus voltage limits violations
%   EF = CPF_VLIM_EVENT(CB_DATA, CX)
%
%   CPF event function to detect bus voltage limits violations,
%   i.e. Vm <= Vmin or Vm >= Vmax.
%
%   Inputs:
%       CB_DATA : struct of data for callback functions
%       CX : struct containing info about current point (continuation soln)
%
%   Outputs:
%       EF : event function value

%   MATPOWER
%   Copyright (c) 2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Ahmad Abubakar Sadiq, Federal University of Technology Minna, Nigeria
%   and Shrirang Abhyankar, Argonne National Laboratory

%%   This file is not yet part of MATPOWER.
%   It is not yet covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% event function value is 2 nb x 1 vector equal to:
%%      [ Vmin - Vm ]
%%      [ Vm - Vmax ]

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% get updated MPC
d = cb_data;
mpc = cpf_current_mpc(d.mpc_base, d.mpc_target, ...
    d.Ybus, d.Yf, d.Yt, d.ref, d.pv, d.pq, cx.V, cx.lam, d.mpopt);

%% Voltage magnitude violations
v_Vmin = mpc.bus(:, VMIN) - mpc.bus(:, VM); 
v_Vmax = mpc.bus(:, VM) - mpc.bus(:, VMAX); 

%% assemble event function value
ef = [v_Vmin;v_Vmax];
