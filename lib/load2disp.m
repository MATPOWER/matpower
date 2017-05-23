function mpc1 = load2disp(mpc0, fname, idx, voll)
%LOAD2DISP Converts fixed loads to dispatchable.
%   MPC = LOAD2DISP(MPC0);
%   MPC = LOAD2DISP(MPC0, FNAME);
%   MPC = LOAD2DISP(MPC0, FNAME, IDX);
%   MPC = LOAD2DISP(MPC0, FNAME, IDX, VOLL);
%
%   Takes a MATPOWER case file or struct and converts fixed loads to
%   dispatchable loads and returns the resulting case struct. Inputs
%   are as follows:
%
%   MPC0 - File name or struct with initial MATPOWER case.
%
%   FNAME (optional) - Name to use to save resulting MATPOWER case. If empty,
%       the case will not be saved to a file.
%
%   IDX (optional) - Vector of bus indices of loads to be converted. If empty
%       or not supplied, it will convert all loads with positive real
%       power demand.
%
%   VOLL (optional) - Scalar or vector specifying the value of lost
%       load to use as the value for the dispatchable loads. If it is
%       a scalar it is used for all loads, if a vector, the dimension
%       must match that of IDX. Default is $5000 per MWh.

%   MATPOWER
%   Copyright (c) 2010-2016, Power Systems Engineering Research Center (PSERC)
%   by Alberto Lamadrid, PSERC Cornell
%   modified by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define constants
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

mpc = loadcase(mpc0);

%% which loads will be converted?
if nargin < 3 || isempty(idx)
    idx = find(mpc.bus(:, PD) > 0); %% by default, all with PD > 0
end

%% set some defaults
voll0   = 5000;             %% default value of lost load
mBase   = 100;              %% generator MVA base
nld     = length(idx);
v1      = ones(nld, 1);     %% vector of ones
v0      = zeros(nld, 1);    %% vector of zeros

%% gen table
gen = [
    mpc.bus(idx, BUS_I), ...        %% GEN_BUS
    -mpc.bus(idx, PD), ...          %% PG
    -mpc.bus(idx, QD), ...          %% QG
    max(0, -mpc.bus(idx, QD)), ...  %% QMAX
    min(0, -mpc.bus(idx, QD)), ...  %% QMIN
    mpc.bus(idx, VM), ...           %% VG
    mBase * v1, ...                 %% MBASE
    v1, ...                         %% GEN_STATUS
    max(0, -mpc.bus(idx, PD)), ...  %% PMAX
    min(0, -mpc.bus(idx, PD)), ...  %% PMIN
    zeros(nld, 6), ...              %% capability curve
    Inf(nld, 4), ...                %% ramp rates
    zeros(nld, 1), ...              %% participation factor
];
[ng, nc] = size(mpc.gen);
mpc.gen = [mpc.gen; zeros(nld, nc)];
mpc.gen(ng+1:ng+nld, 1:APF) = gen;   %% add dispatchable loads

%% bus table
mpc.bus(idx, [PD, QD]) = 0; %% zero out fixed loads

%% gencost table
nc = size(mpc.gencost, 2);
if nargin < 4
    voll = voll0 * v1;
elseif length(voll) == 1
    voll = voll * v1;
end
gencost = [             %% use a linear, polynomial cost format
    POLYNOMIAL*v1, ...      %% MODEL
    zeros(nld, 2), ...      %% STARTUP, SHUTDOWN
    2 * v1, ...             %% NCOST
    voll, ...               %% COST, linear term
    zeros(nld, nc-5) ...    %% constant term and zero-padding
];
mpc.gencost = [mpc.gencost; gencost];

%% (optional) generator fuel types
if isfield(mpc, 'genfuel') && iscell(mpc.genfuel)
    genfuel = cell(nld, 1);
    for k = 1:nld
        genfuel{k} = 'dl';
    end
    mpc.genfuel = [mpc.genfuel; genfuel];
end

%% (optional) generator unit types
if isfield(mpc, 'gentype') && iscell(mpc.gentype)
    gentype = cell(nld, 1);
    for k = 1:nld
        gentype{k} = 'DL';
    end
    mpc.gentype = [mpc.gentype; gentype];
end

%% save case, if filename is given
if nargin > 1 && ~isempty(fname)
    savecase(fname, mpc, '2');
end
if nargout > 0
    mpc1 = mpc;
end
