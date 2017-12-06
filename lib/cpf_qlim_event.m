function ef = cpf_qlim_event(cb_data, cx)
%CPF_QLIM_EVENT  Event function to detect gen reactive power limit violations
%   EF = CPF_QLIM_EVENT(CB_DATA, CX)
%
%   CPF event function to detect generator reactive power limit violations,
%   i.e. Qg <= Qmin or Qg >= Qmax.
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
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% event function value is 2 ng x 1 vector equal to:
%%      [ Qg - Qmax ]
%%      [ Qmin - Qg ]

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

%% compute Qg violations for on-line gens, not at PQ buses
nb = size(mpc.bus, 1);
ng = size(mpc.gen, 1);
on = find(mpc.gen(:, GEN_STATUS) > 0 & ...  %% which generators are on?
          mpc.bus(mpc.gen(:, GEN_BUS), BUS_TYPE) ~= PQ);  %% ... and are not PQ buses
gbus = mpc.gen(on, GEN_BUS);                %% what buses are they at?
ngon = size(on, 1);

%% build connection matrix, element i, j is 1 if gen on(i) at bus j is ON
Cg = sparse((1:ngon)', gbus, ones(ngon, 1), ngon, nb);
C = Cg * Cg';

%% violations are based on total violation at bus, not individual violations
%% (see https://github.com/MATPOWER/matpower/issues/26)
v_Qmax = NaN(ng, 1);
v_Qmin = v_Qmax;
v_Qmax(on) = C * (mpc.gen(on, QG) - mpc.gen(on, QMAX));
v_Qmin(on) = C * (mpc.gen(on, QMIN) - mpc.gen(on, QG));

%% assemble event function value
ef = [v_Qmax; v_Qmin];
