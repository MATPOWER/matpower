function [J, Ybus, Yf, Yt] = makeJac(baseMVA, bus, branch, gen)
%MAKEJAC  Forms the power flow Jacobian.
%   J = MAKEJAC(MPC)
%   J = MAKEJAC(BASEMVA, BUS, BRANCH, GEN)
%   [J, YBUS, YF, YT] = MAKEJAC(MPC)
%
%   Returns the power flow Jacobian and, optionally, the system admittance
%   matrices. Inputs can be a MATPOWER case struct or individual BASEMVA,
%   BUS, BRANCH and GEN values. Bus numbers must be consecutive beginning
%   at 1 (internal ordering).
%
%   Note: This function builds the Jacobian from scratch, rebuilding the
%         YBUS matrix in the process. You probably don't want to use this
%         in performance critical code.
%
%   See also MAKEYBUS, EXT2INT

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2010 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

if nargin < 4
    mpc     = baseMVA;
    baseMVA = mpc.baseMVA;
    bus     = mpc.bus;
    branch  = mpc.branch;
    gen     = mpc.gen;
end

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% get bus index lists of each type of bus
[ref, pv, pq] = bustypes(bus, gen);

%% build Ybus
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

%% extract voltage
V = bus(:, VM) .* exp(sqrt(-1) * pi/180 * bus(:, VA));
on = find(gen(:, GEN_STATUS) > 0);      %% which generators are on?
gbus = gen(on, GEN_BUS);                %% what buses are they at?
V(gbus) = gen(on, VG) ./ abs(V(gbus)).* V(gbus);

%% build Jacobian
[dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, V);
j11 = real(dSbus_dVa([pv; pq], [pv; pq]));
j12 = real(dSbus_dVm([pv; pq], pq));
j21 = imag(dSbus_dVa(pq, [pv; pq]));
j22 = imag(dSbus_dVm(pq, pq));

J = [   j11 j12;
        j21 j22;    ];
