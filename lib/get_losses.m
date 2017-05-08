function [loss, fchg, tchg, dloss_dV, dchg_dVm] = get_losses(baseMVA, bus, branch)
%GET_LOSSES   Returns series losses (and reactive injections) per branch.
%
%   LOSS = GET_LOSSES(RESULTS)
%   LOSS = GET_LOSSES(BASEMVA, BUS, BRANCH)
%
%   [LOSS, CHG] = GET_LOSSES(RESULTS)
%   [LOSS, FCHG, TCHG] = GET_LOSSES(RESULTS)
%   [LOSS, FCHG, TCHG, DLOSS_DV] = GET_LOSSES(RESULTS)
%   [LOSS, FCHG, TCHG, DLOSS_DV, DCHG_DVM] = GET_LOSSES(RESULTS)
%
%   Computes branch series losses, and optionally reactive injections from
%   line charging, as functions of bus voltages and branch parameters, using the
%   following formulae:
%
%       loss = abs( Vf / tau - Vt ) ^ 2 / (Rs - j Xs)
%       fchg = abs( Vf / tau ) ^ 2 * Bc / 2
%       tchg = abs( Vt ) ^ 2 * Bc / 2
%
%   Optionally, computes the partial derivatives of the line losses with
%   respect to voltage angles and magnitudes.
%
%   Input:
%       RESULTS - a MATPOWER case struct with bus voltages corresponding to
%                 a valid power flow solution.
%                 (Can optionally be specified as individual fields BASEMVA,
%                  BUS, and BRANCH.)
%
%   Output(s):
%       LOSS - complex NL x 1 vector of losses (in MW), where NL is the number
%              of branches in the system, representing only the losses in the
%              series impedance element of the PI model for each branch.
%       CHG -  NL x 1 vector of total reactive injection for each line
%              (in MVAr), representing the line charging injections of both
%              of the shunt elements of PI model for each branch.
%       FCHG - Same as CHG, but for the element at the "from" end of the
%              branch only.
%       TCHG - Same as CHG, but for the element at the "to" end of the branch.
%       DLOSS_DV - Struct with partial derivatives of LOSS with respect to bus
%              voltages, with fields:
%           .a  - Partial with respect to bus voltage angles.
%           .m  - Partial with respect to bus voltage magnitudes.
%       DCHG_DVM - Struct with partial derivatives of FCHG and TCHG with
%              respect to bus voltage magnitudes, with fields:
%           .f  - Partial of FCHG with respect to bus voltage magnitudes.
%           .t  - Partial of TCHG with respect to bus voltage magnitudes.
%
%   Example:
%       results = runpf(mycase);
%       total_system_real_losses = sum(real(get_losses(results)));
%
%       [loss, fchg, tchg, dloss_dV] = get_losses(results);

%   MATPOWER
%   Copyright (c) 1996-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% default arguments
if isstruct(baseMVA)
    mpc = baseMVA;
    [baseMVA, bus, branch] = deal(mpc.baseMVA, mpc.bus, mpc.branch);
end

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% create map of external bus numbers to bus indices
i2e = bus(:, BUS_I);
e2i = sparse(max(i2e), 1);
e2i(i2e) = (1:size(bus, 1))';
out = find(branch(:, BR_STATUS) == 0);          %% out-of-service branches

%% sizes of things
nb = size(bus, 1);      %% number of buses
nl = size(branch, 1);   %% number of branches

%% construct complex bus voltage vector
V = bus(:, VM) .* exp(sqrt(-1) * pi/180 * bus(:, VA));

%% parameters
Cf = sparse(1:nl, e2i(branch(:, F_BUS)), branch(:, BR_STATUS), nl, nb);
Ct = sparse(1:nl, e2i(branch(:, T_BUS)), branch(:, BR_STATUS), nl, nb);
tap = ones(nl, 1);                              %% default tap ratio = 1 for lines
xfmr = find(branch(:, TAP));                    %% indices of transformers
tap(xfmr) = branch(xfmr, TAP);                  %% include transformer tap ratios
tap = tap .* exp(1j*pi/180 * branch(:, SHIFT)); %% add phase shifters
A = spdiags(1 ./ tap, 0, nl, nl) * Cf - Ct;
Ysc = 1 ./ (branch(:, BR_R) - 1j * branch(:, BR_X));
Vdrop = A * V;      %% vector of voltage drop across series impedance element
loss = baseMVA * Ysc .* Vdrop .* conj(Vdrop);
% loss = baseMVA * abs(V(e2i(branch(:, F_BUS))) ./ tap - V(e2i(branch(:, T_BUS)))) .^ 2 ./ ...
%             (branch(:, BR_R) - 1j * branch(:, BR_X));
% loss(out) = 0;

if nargout > 1
    Vf = Cf * V;
    Vt = Ct * V;
    fchg = real(baseMVA / 2 * branch(:, BR_B) .* Vf .* conj(Vf) ./ (tap .* conj(tap)));
    tchg = real(baseMVA / 2 * branch(:, BR_B) .* Vt .* conj(Vt));
%     fchg = abs(V(e2i(branch(:, F_BUS))) ./ tap) .^ 2 .* branch(:, BR_B) * baseMVA / 2;
%     tchg = abs(V(e2i(branch(:, T_BUS)))       ) .^ 2 .* branch(:, BR_B) * baseMVA / 2;
    fchg(out) = 0;
    tchg(out) = 0;

    if nargout == 2
        fchg = fchg + tchg;
    end
end

if nargout > 3
    B = spdiags(A * V, 0, nl, nl) * conj(A) * spdiags(conj(V), 0, nb, nb);
    dYsc = spdiags(Ysc, 0, nl, nl);
    dloss_dV = struct(...
        'a', -1j * baseMVA * dYsc * (B - conj(B)), ...
        'm',       baseMVA * dYsc * (B + conj(B)) * spdiags(1 ./ abs(V), 0, nb, nb) ...
    );
    if nargout > 4
        Bc = spdiags(branch(:, BR_B), 0, nl, nl);
        tt = spdiags(1 ./ (tap .* conj(tap)), 0, nl, nl);
        dchg_dVm = struct(...
            'f', baseMVA * Bc * tt * spdiags(Cf * bus(:, VM), 0, nl, nl) * Cf, ...
            't', baseMVA * Bc      * spdiags(Ct * bus(:, VM), 0, nl, nl) * Ct ...
        );
    end
end
