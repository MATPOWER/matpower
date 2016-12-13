function [Bp, Bpp] = makeB(baseMVA, bus, branch, alg)
%MAKEB   Builds the FDPF matrices, B prime and B double prime.
%   [BP, BPP] = MAKEB(MPC, ALG)
%   [BP, BPP] = MAKEB(BASEMVA, BUS, BRANCH, ALG)
%
%   Returns the two matrices B prime and B double prime used in the fast
%   decoupled power flow. Does appropriate conversions to p.u. ALG is either
%   'FDXB' or 'FDBX', the corresponding value of MPOPT.pf.alg option
%   specifying the power flow algorithm.
%   Bus numbers must be consecutive beginning at 1 (i.e. internal ordering).
%
%   Note: For backward compatibility, ALG can also take on a value of
%   2 or 3, corresponding to values of the old PF_ALG option. This usage
%   is deprecated and will be removed in a future version.
%
%   Example:
%       [Bp, Bpp] = makeB(baseMVA, bus, branch, 'FDXB');
%
%   See also FDPF.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% extract from MPC if necessary
if nargin < 3
    mpc     = baseMVA;
    if nargin == 2
        alg = bus;
    end
    baseMVA = mpc.baseMVA;
    bus     = mpc.bus;
    branch  = mpc.branch;
end

%% constants
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of lines

%% backward compatiblility (deprecated)
if ~ischar(alg)
    if alg == 2
        alg = 'FDXB';
    elseif alg == 3
        alg = 'FDBX';
    end
end

%% check for valid ALG value
alg = upper(alg);
if ~strcmp(alg, 'FDXB') && ~strcmp(alg, 'FDBX')
    error('makeB: ''%s'' is not a valid value for ALG', alg);
end

%% define named indices into bus, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%%-----  form Bp (B prime)  -----
temp_branch = branch;                       %% modify a copy of branch
temp_bus = bus;                             %% modify a copy of bus
temp_bus(:, BS) = zeros(nb, 1);             %% zero out shunts at buses
temp_branch(:, BR_B) = zeros(nl, 1);        %% zero out line charging shunts
temp_branch(:, TAP) = ones(nl, 1);          %% cancel out taps
if strcmp(alg, 'FDXB')                      %% if XB method
    temp_branch(:, BR_R) = zeros(nl, 1);        %% zero out line resistance
end
Bp = -imag( makeYbus(baseMVA, temp_bus, temp_branch) );

%%-----  form Bpp (B double prime)  -----
if nargout == 2
    temp_branch = branch;                       %% modify a copy of branch
    temp_branch(:, SHIFT) = zeros(nl, 1);       %% zero out phase shifters
    if strcmp(alg, 'FDBX')                      %% if BX method
        temp_branch(:, BR_R) = zeros(nl, 1);        %% zero out line resistance
    end
    Bpp = -imag( makeYbus(baseMVA, bus, temp_branch) );
end
