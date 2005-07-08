function [Bbus, Bf, Pbusinj, Pfinj] = makeBdc(baseMVA, bus, branch)
%MAKEBBUS   Builds the B matrices and phase shift injections for DC power flow.
%   [Bbus, Bf, Pbusinj, Pfinj] = makeBdc(baseMVA, bus, branch) returns the
%   B matrices and phase shift injection vectors needed for a DC power flow.
%   The bus real power injections are related to bus voltage angles by
%       P = Bbus * Va + Pbusinj
%   The real power flows at the from end the lines are related to the bus
%   voltage angles by
%       Pf = Bf * Va + Pfinj
%   Does appropriate conversions to p.u.

%   MATPOWER
%   $Id$
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Autonoma de Manizales
%   and Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% constants
j = sqrt(-1);
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of lines

%% define named indices into bus, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% check that bus numbers are equal to indices to bus (one set of bus numbers)
if any(bus(:, BUS_I) ~= [1:nb]')
    error('buses must appear in order by bus number')
end

%% for each branch, compute the elements of the branch B matrix and the phase
%% shift "quiescent" injections, where
%%
%%      | Pf |   | Bff  Bft |   | Vaf |   | Pfinj |
%%      |    | = |          | * |     | + |       |
%%      | Pt |   | Btf  Btt |   | Vat |   | Ptinj |
%%
stat = branch(:, BR_STATUS);                    %% ones at in-service branches
b = stat ./ branch(:, BR_X);                    %% series susceptance
tap = ones(nl, 1);                              %% default tap ratio = 1
i = find(branch(:, TAP));                       %% indices of non-zero tap ratios
tap(i) = branch(i, TAP);                        %% assign non-zero tap ratios
b = b ./ tap;

%% build Bbus
f = branch(:, F_BUS);                           %% list of "from" buses
t = branch(:, T_BUS);                           %% list of "to" buses
Cf = sparse(f, 1:nl, ones(nl, 1), nb, nl);      %% connection matrix for line & from buses
Ct = sparse(t, 1:nl, ones(nl, 1), nb, nl);      %% connection matrix for line & to buses
Bbus =  Cf * spdiags(b, 0, nl, nl) * Cf' + ...  %% Bff term of branch admittance
        Cf * spdiags(-b, 0, nl, nl) * Ct' + ... %% Bft term of branch admittance
        Ct * spdiags(-b, 0, nl, nl) * Cf' + ... %% Btf term of branch admittance
        Ct * spdiags(b, 0, nl, nl) * Ct';       %% Btt term of branch admittance

%% build phase shift injection vectors
Pfinj = b .* (branch(:, SHIFT) * pi/180);       %% injected at the from bus ...
    % Ptinj = -Pfinj;                           %% ... and extracted at the to bus
Pbusinj = (Cf - Ct) * Pfinj;                    %% Pbusinj = Cf * Pfinj + Ct * Ptinj;

%% Build Bf such that Bf * Va is the vector of real branch powers injected
%% at each branch's "from" bus
if nargout > 1
    i = [[1:nl]'; [1:nl]'];     %% double set of row indices    
    Bf = sparse(i, [f; t], [b; -b]);
end

return;
