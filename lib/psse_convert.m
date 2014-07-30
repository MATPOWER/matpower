function [mpc, warns] = psse_convert(warns, data, verbose)
%PSSE_CONVERT  Converts data read from PSS/E RAW file to MATPOWER case.
%   [MPC, WARNINGS] = PSSE_CONVERT(WARNINGS, DATA)
%   [MPC, WARNINGS] = PSSE_CONVERT(WARNINGS, DATA, VERBOSE)
%
%   Converts data read from a version RAW data file into a
%   MATPOWER case struct.
%
%   Input:
%       WARNINGS :  cell array of strings containing accumulated
%                   warning messages
%       DATA : struct read by PSSE_READ (see PSSE_READ for details).
%       VERBOSE :   1 to display progress info, 0 (default) otherwise
%
%   Output:
%       MPC : a MATPOWER case struct created from the PSS/E data
%       WARNINGS :  cell array of strings containing updated accumulated
%                   warning messages
%
%   See also PSSE_READ.

%   MATPOWER
%   $Id$
%   by Yujia Zhu, PSERC ASU
%   and Ray Zimmerman, PSERC Cornell
%   Based on mpraw2mp.m, written by: Yujia Zhu, Jan 2014, yzhu54@asu.edu.
%   Copyright (c) 2014 by Power System Engineering Research Center (PSERC)
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

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% options
sort_buses = 1;

%% defaults
if nargin < 3
    verbose = 0;
end
haveVlims = 0;
Vmin = 0.9;
Vmax = 1.1;

%% PSS/E id record
%% IC, SBASE, REV, XFRRAT, NXFRAT, BASFRQ
baseMVA = data.id.SBASE;
rev = data.id.REV;

%% PSS/E bus data
%% v1: (http://www.ee.washington.edu/research/pstca/formats/pti.txt)
%%  I, IDE, PL, QL, GL, BL, IA, VM, VA, 'NAME', BASKL, ZONE
%% v29-30:
%%  I, 'NAME', BASKV, IDE, GL, BL, AREA, ZONE, VM, VA, OWNER
%% v31:
%%  I, 'NAME', BASKV, IDE, AREA, ZONE, OWNER, VM, VA
%% v33:
%%  I, 'NAME', BASKV, IDE, AREA, ZONE, OWNER, VM, VA, NVHI, NVLO, EVHI, EVLO
%% Note: Some v33 files seem to follow v31 format
%% NVHI/NVLO: Normal voltage high/low limits
%% EVHI/EVLO: Emergency voltage high/low limits
numbus = data.bus.num;
[nb, ncols] = size(numbus); %% number of buses, number of cols
bus = zeros(nb, VMIN);      %% initialize bus matrix
if rev == 1
    bus_name_col = 10;
else
    bus_name_col = 2;
end
if sort_buses
    [numbus, i] = sortrows(numbus, 1);
    bus_name = data.bus.txt(i, bus_name_col);
else
    bus_name = data.bus.txt(:, bus_name_col);
end
if rev == 1
    bus(:, [BUS_I BUS_TYPE PD QD GS BS BUS_AREA VM VA BASE_KV ZONE]) = ...
        numbus(:, [1:9 11:12]);
elseif rev < 31     %% includes GL, BL
    bus(:, [BUS_I BASE_KV BUS_TYPE GS BS BUS_AREA ZONE VM VA]) = ...
        numbus(:, [1 3 4 5 6 7 8 9 10]);
else
    bus(:, [BUS_I BASE_KV BUS_TYPE BUS_AREA ZONE VM VA]) = ...
        numbus(:, [1 3 4 5 6 8 9]);
    if ncols >= 11 && all(all(~isnan(numbus(:, [10 11]))))
        haveVlims = 1;
        bus(:, [VMAX VMIN]) = numbus(:, [10 11]);
    end
end
if ~haveVlims  %% add default voltage magnitude limits if not provided
    warns{end+1} = sprintf('Using default voltage magnitude limits: VMIN = %g p.u., VMAX = %g p.u.', Vmin, Vmax);
    if verbose
        fprintf('WARNING: No bus voltage magnitude limits provided.\n         Using defaults: VMIN = %g p.u., VMAX = %g p.u.\n', Vmin, Vmax);
    end
    bus(:, VMIN) = Vmin;
    bus(:, VMAX) = Vmax;
end

%% create map of external bus numbers to bus indices
i2e = bus(:, BUS_I);
e2i = sparse(i2e, ones(nb, 1), 1:nb, max(i2e), 1);

%% PSS/E load data
%% v29-31:
%%  I, ID, STATUS, AREA, ZONE, PL, QL, IP, IQ, YP, YQ, OWNER
%% v33:
%%  I, ID, STATUS, AREA, ZONE, PL, QL, IP, IQ, YP, YQ, OWNER, SCALE, INTRPT
%% Note: Some v33 files seem to end with SCALE
if rev > 1
    nld = size(data.load.num, 1);
    loadbus = e2i(data.load.num(:,1));
    %% PSS/E loads are divided into:
    %%  1. constant MVA, (I=1)
    %%  2. constant current (I=I)
    %%  3. constant reactance/resistance (I = I^2)
    %% NOTE: reactive power component of constant admittance load is negative
    %%       quantity for inductive load and positive for capacitive load
    Pd = data.load.num(:,6) + data.load.num(:,8) .* bus(loadbus, VM) ...
            + data.load.num(:,10) .* bus(loadbus, VM).^2;
    Qd = data.load.num(:,7) + data.load.num(:,9) .* bus(loadbus, VM) ...
            - data.load.num(:,11) .* bus(loadbus, VM).^2;
    Cld = sparse(1:nld, loadbus, data.load.num(:,3), nld, nb);    %% only in-service-loads
    bus(:, [PD QD]) = Cld' * [Pd Qd];
end

%% PSS/E fixed shunt data
%% v31-33:
%%  I, ID, STATUS, GL, BL
if isfield(data, 'shunt')   %% rev > 30
    nsh = size(data.shunt.num, 1);
    shuntbus = e2i(data.shunt.num(:,1));
    Csh = sparse(1:nsh, shuntbus, data.shunt.num(:,3), nsh, nb);  %% only in-service shunts
    bus(:, [GS BS]) = Csh' * data.shunt.num(:, 4:5);
end

%% PSS/E switched shunt data
%% v1: (http://www.ee.washington.edu/research/pstca/formats/pti.txt)
%%  I, MODSW, VSWHI, VSWLO, SWREM, BINIT, N1, B1, N2, B2, ... N8, B8
%% v29:
%%  I, MODSW, VSWHI, VSWLO, SWREM, RMIDNT, BINIT, N1, B1, N2, B2, ... N8, B8
%% v30:
%%  I, MODSW, VSWHI, VSWLO, SWREM, RMPCT, 'RMIDNT', BINIT, N1, B1, N2, B2, ... N8, B8
%% v31:
%%  I, MODSW, VSWHI, VSWLO, SWREM, RMPCT, 'RMIDNT', BINIT, N1, B1, N2, B2, ... N8, B8
%% v32:
%%  ? (we will assume v32 is the same as v33, until we find out otherwise)
%% v33:
%%  I, MODSW, ADJM, STAT, VSWHI, VSWLO, SWREM, RMPCT, 'RMIDNT', BINIT, N1, B1, N2, B2, ... N8, B8
nswsh = size(data.swshunt.num, 1);
swshuntbus = e2i(data.swshunt.num(:,1));
Cswsh = sparse(1:nswsh, swshuntbus, 1, nswsh, nb);
if rev == 1
    bus(:, BS) = bus(:, BS) + Cswsh' * data.swshunt.num(:, 6);
elseif rev <= 29
    bus(:, BS) = bus(:, BS) + Cswsh' * data.swshunt.num(:, 7);
elseif rev < 32
    bus(:, BS) = bus(:, BS) + Cswsh' * data.swshunt.num(:, 8);
else
    bus(:, BS) = bus(:, BS) + Cswsh' * data.swshunt.num(:, 10);
end

%% PSS/E branch data
%% v1: (http://www.ee.washington.edu/research/pstca/formats/pti.txt)
%%  I,J,CKT,R,X,B,RATEA,RATEB,RATEC,RATIO,ANGLE,GI,BI,GJ,BJ,ST
%% v29-30:
%%  I,J,CKT,R,X,B,RATEA,RATEB,RATEC,GI,BI,GJ,BJ,ST,LEN,O1,F1,...,O4,F4
%% v31-33:
%%  I,J,CKT,R,X,B,RATEA,RATEB,RATEC,GI,BI,GJ,BJ,ST,MET,LEN,O1,F1,...,O4,F4
nbr = size(data.branch.num, 1);
branch = zeros(nbr, ANGMAX);
branch(:, ANGMIN) = -360;
branch(:, ANGMAX) = 360;
branch(:, [F_BUS BR_R BR_X BR_B RATE_A RATE_B RATE_C]) = ...
    data.branch.num(:, [1 4 5 6 7 8 9]);
branch(:, T_BUS) = abs(data.branch.num(:, 2));  %% can be negative to indicate metered end
if rev == 1
    branch(:, BR_STATUS) = data.branch.num(:, 16);
    branch(~isnan(data.branch.num(:, 10)), TAP) = ...
        data.branch.num(~isnan(data.branch.num(:, 10)), 10);
    branch(~isnan(data.branch.num(:, 11)), SHIFT) = ...
        data.branch.num(~isnan(data.branch.num(:, 11)), 11);
else
    branch(:, BR_STATUS)    = data.branch.num(:, 14);
end
%% integrate branch shunts (explicit shunts, not line-charging)
ibr = (1:nbr)';
fbus = e2i(branch(:, F_BUS));
tbus = e2i(branch(:, T_BUS));
nzf = find(fbus);               %% ignore branches with bad bus numbers
nzt = find(tbus);
if length(nzf) < nbr
    warns{end+1} = sprintf('%d branches have bad ''from'' bus numbers', nbr-length(nzf));
    if verbose
        fprintf('WARNING: %d branches have bad ''from'' bus numbers\n', nbr-length(nzf));
    end
end
if length(nzt) < nbr
    warns{end+1} = sprintf('%d branches have bad ''to'' bus numbers', nbr-length(nzt));
    if verbose
        fprintf('WARNING: %d branches have bad ''to'' bus numbers\n', nbr-length(nzt));
    end
end
Cf = sparse(ibr(nzf), fbus(nzf), branch(nzf, BR_STATUS), nbr, nb);  %% only in-service branches
Ct = sparse(ibr(nzt), tbus(nzt), branch(nzt, BR_STATUS), nbr, nb);  %% only in-service branches
if rev == 1
    bus(:, [GS BS]) = bus(:, [GS BS]) + ...
        Cf' * data.branch.num(:, 12:13)*baseMVA + ...
        Ct' * data.branch.num(:, 14:15)*baseMVA;
else
    bus(:, [GS BS]) = bus(:, [GS BS]) + ...
        Cf' * data.branch.num(:, 10:11)*baseMVA + ...
        Ct' * data.branch.num(:, 12:13)*baseMVA;
end

%% PSS/E generator data
%% v1: (http://www.ee.washington.edu/research/pstca/formats/pti.txt)
%%  I,ID,PG,QG,QT,QB,VS,IREG,MBASE,ZR,ZX,RT,XT,GTAP,STAT,RMPCT,PT,PB
%% v29-30:
%%  I,ID,PG,QG,QT,QB,VS,IREG,MBASE,ZR,ZX,RT,XT,GTAP,STAT,RMPCT,PT,PB,O1,F1,...,O4,F4
%% v31-33:
%%  I,ID,PG,QG,QT,QB,VS,IREG,MBASE,ZR,ZX,RT,XT,GTAP,STAT,RMPCT,PT,PB,O1,F1,...,O4,F4,WMOD,WPF
ng = size(data.gen.num, 1);
genbus = e2i(data.gen.num(:,1));
gen = zeros(ng, APF);
gen = zeros(ng, MU_QMIN);
gen(:, [GEN_BUS PG QG QMAX QMIN VG MBASE GEN_STATUS PMAX PMIN]) = ...
    data.gen.num(:, [1 3 4 5 6 7 9 15 17 18]);

%% PSS/E transformer data
if rev > 1
    [transformer, bus, warns, bus_name] = psse_convert_xfmr(warns, data.trans2.num, data.trans3.num, verbose, baseMVA, bus, bus_name);
    branch = [branch; transformer];
end

%% PSS/E two terminal HVDC line data
dcline = psse_convert_hvdc(data.twodc.num, bus);

%% assemble MPC
mpc = struct( ...
    'baseMVA',  baseMVA, ...
    'bus', bus, ...
    'bus_name', {bus_name}, ...
    'branch', branch, ...
    'gen', gen ...
);
if ~isempty(dcline)
    mpc.dcline = dcline;
    mpc = toggle_dcline(mpc, 'on');
end
