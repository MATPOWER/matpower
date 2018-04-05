function [J, Ybus, Yf, Yt] = makeJac(baseMVA, bus, branch, gen, fullJac)
%MAKEJAC  Forms the power flow Jacobian.
%   J = MAKEJAC(MPC)
%   J = MAKEJAC(MPC, FULLJAC)
%   J = MAKEJAC(BASEMVA, BUS, BRANCH, GEN)
%   J = MAKEJAC(BASEMVA, BUS, BRANCH, GEN, FULLJAC)
%   [J, YBUS, YF, YT] = MAKEJAC(MPC)
%
%   Returns the power flow Jacobian and, optionally, the system admittance
%   matrices. Inputs can be a MATPOWER case struct or individual BASEMVA,
%   BUS, BRANCH and GEN values. Bus numbers must be consecutive beginning
%   at 1 (i.e. internal ordering). If the FULLJAC argument is present and
%   true, it returns the full Jacobian (sensitivities of all bus injections
%   w.r.t all voltage angles/magnitudes) as opposed to the reduced version
%   used in the Newton power flow updates. The units for all quantities are
%   in per unit with radians for voltage angles.
%
%   Note: This function builds the Jacobian from scratch, rebuilding the
%         YBUS matrix in the process. You probably don't want to use this
%         in performance critical code.
%
%   See also MAKEYBUS, EXT2INT

%   MATPOWER
%   Copyright (c) 1996-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 4
    mpc     = baseMVA;
    if nargin > 1
        fullJac = bus;
    else
        fullJac = 0;
    end
    baseMVA = mpc.baseMVA;
    bus     = mpc.bus;
    branch  = mpc.branch;
    gen     = mpc.gen;
elseif nargin < 5
    fullJac = 0;
end

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% build Ybus
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

%% extract voltage
V = bus(:, VM) .* exp(sqrt(-1) * pi/180 * bus(:, VA));

%% make sure we use generator setpoint voltage for PV and slack buses
on = find(gen(:, GEN_STATUS) > 0);      %% which generators are on?
gbus = gen(on, GEN_BUS);                %% what buses are they at?
k = find(bus(gbus, BUS_TYPE) == PV | bus(gbus, BUS_TYPE) == REF);
V(gbus(k)) = gen(on(k), VG) ./ abs(V(gbus(k))).* V(gbus(k));

%% build Jacobian
[dSbus_dVa, dSbus_dVm] = dSbus_dV(Ybus, V);
if fullJac
    j11 = real(dSbus_dVa);
    j12 = real(dSbus_dVm);
    j21 = imag(dSbus_dVa);
    j22 = imag(dSbus_dVm);
else
    %% get bus index lists of each type of bus
    [ref, pv, pq] = bustypes(bus, gen);

    j11 = real(dSbus_dVa([pv; pq], [pv; pq]));
    j12 = real(dSbus_dVm([pv; pq], pq));
    j21 = imag(dSbus_dVa(pq, [pv; pq]));
    j22 = imag(dSbus_dVm(pq, pq));
end

J = [   j11 j12;
        j21 j22;    ];
