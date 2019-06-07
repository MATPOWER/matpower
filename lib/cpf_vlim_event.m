function ef = cpf_vlim_event(cb_data, cx)
%CPF_VLIM_EVENT  Event function to detect bus voltage limit violations
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
%   Copyright (c) 2016-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Ahmad Abubakar Sadiq, Federal University of Technology Minna, Nigeria
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% event function value is 2 nb x 1 vector equal to:
%%      [ Vmin - Vm ]
%%      [ Vm - Vmax ]

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

%% get updated MPC
d = cb_data;
mpc = cpf_current_mpc(d.mpc_base, d.mpc_target, ...
    d.Ybus, d.Yf, d.Yt, d.ref, d.pv, d.pq, cx.V, cx.lam, d.mpopt);

%% voltage magnitude violations
v_Vmin = mpc.bus(:, VMIN) - mpc.bus(:, VM);
v_Vmax = mpc.bus(:, VM) - mpc.bus(:, VMAX);

%% assemble event function value
ef = [v_Vmin;v_Vmax];
