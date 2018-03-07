function delta = calc_branch_angle(mpc)
%%% calculates the angle difference (in degrees) accross all active 
%%% branches in the matpower case. Angles are calculated a the difference 
%%% between the FROM bus and the TO bus.
%%%
%%% INPUT
%%%         mpc     matpower case with fields bus and branch.
%%% OUTPUT
%%%         delta   vector of size nl x 1. Letting the angle of the FROM 
%%%                 bus be Af and the angle of the TO bus At, the values in
%%%                 delta are Af-At. Note that for out of service branches
%%%                 the value will be 0 by default.
%%%
%%% Eran Schweitzer 2018. 
%%% Eran.Schweitzer@asu.edu
%% define named indices into data matrices
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

%% 
status = mpc.branch(:,BR_STATUS);
nl   = size(mpc.branch,1);
nb   = size(mpc.bus,1);
nmap = sparse(mpc.bus(:,BUS_I),1,1:nb,nb,1); % Map from external numbering to 1:nb numbering
bf   = full(nmap(mpc.branch(:,F_BUS)));      % FROM bus ids
bt   = full(nmap(mpc.branch(:,T_BUS)));      % TO bus ids
delta = sparse([1:nl,1:nl]',[bf;bt], [ones(nl,1).*status;-ones(nl,1).*status], nl,nb)*mpc.bus(:,VA); %angle differences