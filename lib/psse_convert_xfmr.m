function [xfmr, bus, warns, bus_name] = psse_convert_xfmr(warns, trans2, trans3, verbose, baseMVA, bus, bus_name)
%PSSE_CONVERT_XFMR Convert transformer data from PSS/E RAW to MATPOWER
%   [XFMR, BUS, WARNINGS] = PSSE_CONVERT_XFMR(WARNINGS, TRANS2, TRANS3, ...
%                                   VERBOSE, BASEMVA, BUS)
%   [XFMR, BUS, WARNINGS, BUS_NAME] = PSSE_CONVERT_XFMR(WARNINGS, TRANS2, ...
%                                   TRANS3, VERBOSE, BASEMVA, BUS, BUS_NAME)
%
%   Convert all transformer data read from a PSS/E RAW data file
%   into MATPOWER format. Returns a branch matrix corresponding to
%   the transformers and an updated bus matrix, with additional buses
%   added for the star points of three winding transformers.
%
%   Inputs:
%       WARNINGS :  cell array of strings containing accumulated
%                   warning messages
%       TRANS2  : matrix of raw two winding transformer data returned
%                 by PSSE_READ in data.trans2.num
%       TRANS3  : matrix of raw three winding transformer data returned
%                 by PSSE_READ in data.trans3.num
%       VERBOSE :   1 to display progress info, 0 (default) otherwise
%       BASEMVA : system MVA base
%       BUS     : MATPOWER bus matrix
%       BUS_NAME: (optional) cell array of bus names
%
%   Outputs:
%       XFMR    : MATPOWER branch matrix of transformer data
%       BUS     : updated MATPOWER bus matrix, with additional buses
%                 added for star points of three winding transformers
%       WARNINGS :  cell array of strings containing updated accumulated
%                   warning messages
%       BUS_NAME: (optional) updated cell array of bus names
%
%   See also PSSE_CONVERT.

%   MATPOWER
%   Copyright (c) 2014-2016, Power Systems Engineering Research Center (PSERC)
%   by Yujia Zhu, PSERC ASU
%   and Ray Zimmerman, PSERC Cornell
%   Based on mptransin.m and mptransficbus.m, written by:
%       Yujia Zhu, Jan 2014, yzhu54@asu.edu.
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% compatibility
use_winding_baseV = 1;  %% If true, will use winding base V for any conversions
                        %% instead of bus base V if they are different.
                        %% Turn off for compatibility with, apparently buggy,
                        %% PSS/E implementation (checked by YZ both for
                        %% PSS/E internal model and conversion from RAW to
                        %% IEEE format)

%% sizes
nb = size(bus, 1);
nt2 = size(trans2, 1);
nt3 = size(trans3, 1);

%%-----  create fictitious buses for star point of 3 winding transformers  -----
starbus = zeros(nt3, VMIN);     %% initialize additional bus matrix
if nt3 > 0
    mb =  max(bus(:, BUS_I));       %% find maximum bus number
    b = 10^ceil(log10(mb+1));       %% start new numbers at next order of magnitude
    wind1bus = trans3(:,1);         %% winding1 bus number
    e2i = sparse(bus(:, BUS_I), ones(nb, 1), 1:nb, max(bus(:, BUS_I)), 1);
    starbus(:, BUS_I) = (b+1:b+nt3)';   %% bus numbers follow originals
    starbus(:, BUS_TYPE) = PQ;
    starbus(trans3(:,12)==0, BUS_TYPE) = NONE;  %% isolated if transformer is off-line
    starbus(:, VA) = trans3(:,31);  %% VA = star point voltage angle
    starbus(:, VM) = trans3(:,30);  %% VM = star point voltage magnitude (PU)
    starbus(:, [BUS_AREA, ZONE]) = bus(e2i(wind1bus), [7,11]); %% wind1 bus area, zone
    starbus(:, BASE_KV) = 1;        %% baseKV = 1 kV (RayZ: why?)
    starbus(:, VMAX) = 1.1;
    starbus(:, VMIN) = 0.9;
end

%% two winding transformer data
[tf,fbus] = ismember(trans2(:,1), bus(:, BUS_I));      %% I
[tf,tbus] = ismember(trans2(:,2), bus(:, BUS_I));      %% J
%% check for bad bus numbers
k = find(fbus == 0 | tbus == 0);
if ~isempty(k)
    warns{end+1} = sprintf('Ignoring %d two-winding transformers with bad bus numbers', length(k));
    if verbose
        fprintf('WARNING: Ignoring %d two-winding transformers with bad bus numbers', length(k));
    end
    fbus(k) = [];
    tbus(k) = [];
    trans2(k, :) = [];
    nt2 = nt2 - length(k);
end
if use_winding_baseV
    Zbs = bus(:, BASE_KV).^2 / baseMVA;     %% system impedance base
end
if nt2 > 0
    cw2 = find(trans2(:,5) == 2);   %% CW = 2
    cw3 = find(trans2(:,5) == 3);   %% CW = 3
    cw23 = [cw2;cw3];               %% CW = 2 or 3
    cz2 = find(trans2(:,6) == 2);   %% CZ = 2
    cz3 = find(trans2(:,6) == 3);   %% CZ = 3
    cz23 = [cz2;cz3];               %% CZ = 2 or 3

    R = trans2(:,21);
    X = trans2(:,22);
    if use_winding_baseV
        Zb = ones(nt2, 1);
        Zb(cz23) = trans2(cz23,25).^2 ./ trans2(cz23,23);
    end
    R(cz3) = 1e-6 * R(cz3, 1) ./ trans2(cz3,23);
    X(cz3) = sqrt(X(cz3).^2 - R(cz3).^2);   %% R, X for cz3, pu on winding bases
    if use_winding_baseV
        R(cz23) = R(cz23) .* Zb(cz23) ./ Zbs(fbus(cz23));
        X(cz23) = X(cz23) .* Zb(cz23) ./ Zbs(fbus(cz23));
    else    %% use bus base V (even if winding base V is different)
        R(cz23) = baseMVA * R(cz23, 1) ./ trans2(cz23,23);
        X(cz23) = baseMVA * X(cz23, 1) ./ trans2(cz23,23);
    end
    tap = trans2(:,24) ./ trans2(:,40);                 %% WINDV1/WINDV2
    tap(cw23) = tap(cw23, 1) .* bus(tbus(cw23), BASE_KV)./bus(fbus(cw23), BASE_KV);
    tap(cw3)  = tap(cw3, 1)  .* trans2(cw3,25)./trans2(cw3,41);
    shift = trans2(:, 26);
end

%% three winding transformer data
[tf,ind1] = ismember(trans3(:,1), bus(:, BUS_I));
[tf,ind2] = ismember(trans3(:,2), bus(:, BUS_I));
[tf,ind3] = ismember(trans3(:,3), bus(:, BUS_I));
%% check for bad bus numbers
k = find(ind1 == 0 | ind2 == 0 | ind3 == 0);
if ~isempty(k)
    warns{end+1} = sprintf('Ignoring %d three-winding transformers with bad bus numbers', length(k));
    if verbose
        fprintf('WARNING: Ignoring %d three-winding transformers with bad bus numbers', length(k));
    end
    ind1(k) = [];
    ind2(k) = [];
    ind3(k) = [];
    trans3(k, :) = [];
    starbus(k, :) = [];
    nt3 = nt3 - length(k);
end
%% Each three winding transformer will be converted into 3 branches:
%% The branches will be in the order of
%% # winding1 -> # winding2
%% # winding2 -> # winding3
%% # winding3 -> # winding1
if nt3 > 0
    cw2 = find(trans3(:,5) == 2);   %% CW = 2
    cw3 = find(trans3(:,5) == 3);   %% CW = 3
    cw23 = [cw2;cw3];               %% CW = 2 or 3
    cz2 = find(trans3(:,6) == 2);   %% CZ = 2
    cz3 = find(trans3(:,6) == 3);   %% CZ = 3
    cz23 = [cz2;cz3];               %% CZ = 2 or 3

    tap1 = trans3(:, 32);   %% off nominal tap ratio of branch 1
    tap2 = trans3(:, 48);   %% off nominal tap ratio of branch 2
    tap3 = trans3(:, 64);   %% off nominal tap ratio of branch 3
    tap1(cw23) = tap1(cw23, 1) ./ bus(ind1(cw23), BASE_KV);
    tap2(cw23) = tap2(cw23, 1) ./ bus(ind2(cw23), BASE_KV);
    tap3(cw23) = tap3(cw23, 1) ./ bus(ind3(cw23), BASE_KV);
    tap1(cw3)  = tap1(cw3, 1)  .* trans3(cw3, 33);
    tap2(cw3)  = tap2(cw3, 1)  .* trans3(cw3, 49);
    tap3(cw3)  = tap3(cw3, 1)  .* trans3(cw3, 65);
    shift1 = trans3(:, 34);
    shift2 = trans3(:, 50);
    shift3 = trans3(:, 66);

    %% replace winding base voltage with bus base voltage
    % commented out: Yujia thinks this is wrong
    % trans3(cz3, 33) = bus(ind1(cz3), BASE_KV);
    % trans3(cz3, 49) = bus(ind1(cz3), BASE_KV);
    % trans3(cz3, 65) = bus(ind1(cz3), BASE_KV);

    R12 = trans3(:, 21);
    X12 = trans3(:, 22);
    R23 = trans3(:, 24);
    X23 = trans3(:, 25);
    R31 = trans3(:, 27);
    X31 = trans3(:, 28);
    Zb1 = trans3(:, 33).^2 ./ trans3(:, 23);
    Zb2 = trans3(:, 49).^2 ./ trans3(:, 26);
    Zb3 = trans3(:, 65).^2 ./ trans3(:, 29);

    R12(cz3) = 1e-6 * R12(cz3, 1) ./ trans3(cz3,23);
    X12(cz3) = sqrt(X12(cz3).^2 - R12(cz3).^2);
    R23(cz3) = 1e-6 * R23(cz3, 1) ./ trans3(cz3,26);
    X23(cz3) = sqrt(X23(cz3).^2 - R23(cz3).^2);
    R31(cz3) = 1e-6 * R31(cz3, 1) ./ trans3(cz3,29);
    X31(cz3) = sqrt(X31(cz3).^2 - R31(cz3).^2);

    if use_winding_baseV
        R12(cz23) = R12(cz23) .* Zb1(cz23) ./ Zbs(ind1(cz23));
        X12(cz23) = X12(cz23) .* Zb1(cz23) ./ Zbs(ind1(cz23));
        R23(cz23) = R23(cz23) .* Zb2(cz23) ./ Zbs(ind2(cz23));
        X23(cz23) = X23(cz23) .* Zb2(cz23) ./ Zbs(ind2(cz23));
        R31(cz23) = R31(cz23) .* Zb3(cz23) ./ Zbs(ind3(cz23));
        X31(cz23) = X31(cz23) .* Zb3(cz23) ./ Zbs(ind3(cz23));
    else    %% use bus base V (even if winding base V is different)
        R12(cz23) = baseMVA * R12(cz23, 1) ./ trans3(cz23,23);
        X12(cz23) = baseMVA * X12(cz23, 1) ./ trans3(cz23,23);
        R23(cz23) = baseMVA * R23(cz23, 1) ./ trans3(cz23,26);
        X23(cz23) = baseMVA * X23(cz23, 1) ./ trans3(cz23,26);
        R31(cz23) = baseMVA * R31(cz23, 1) ./ trans3(cz23,29);
        X31(cz23) = baseMVA * X31(cz23, 1) ./ trans3(cz23,29);
    end

    R1 = (R12+R31-R23) ./ 2;
    R2 = (R12+R23-R31) ./ 2;
    R3 = (R31+R23-R12) ./ 2;
    X1 = (X12+X31-X23) ./ 2;
    X2 = (X12+X23-X31) ./ 2;
    X3 = (X31+X23-X12) ./ 2;
end

%%-----  assemble transformer data into MATPOWER branch format  -----
%%% two winding transformers %%%
xfmr2 = zeros(nt2, ANGMAX);
if nt2 > 0
    xfmr2(:, [F_BUS T_BUS]) = trans2(:,[1,2]);
    xfmr2(:, [BR_R BR_X TAP SHIFT]) = [R X tap shift];
    xfmr2(:, [RATE_A RATE_B RATE_C]) = trans2(:,[27,28,29]);
    xfmr2(:, BR_STATUS) = trans2(:,12);
    xfmr2(:, ANGMIN) = -360;
    xfmr2(:, ANGMAX) = 360;
end

%%% three winding transformers %%%
xfmr3 = zeros(3*nt3, ANGMAX);
if nt3 > 0
    idx1 = (1:3:3*nt3)';        %% indices of winding 1
    idx2 = idx1+1;              %% indices of winding 2
    idx3 = idx1+2;              %% indices of winding 3
    %% bus numbers
    xfmr3(idx1, [F_BUS T_BUS]) = [trans3(:,1), starbus(:,1)];
    xfmr3(idx2, [F_BUS T_BUS]) = [trans3(:,2), starbus(:,1)];
    xfmr3(idx3, [F_BUS T_BUS]) = [trans3(:,3), starbus(:,1)];
    %% impedances, tap ratios & phase shifts
    xfmr3(idx1, [BR_R BR_X TAP SHIFT]) = [R1 X1 tap1 shift1];
    xfmr3(idx2, [BR_R BR_X TAP SHIFT]) = [R2 X2 tap2 shift2];
    xfmr3(idx3, [BR_R BR_X TAP SHIFT]) = [R3 X3 tap3 shift3];
    %% ratings
    xfmr3(idx1, [RATE_A RATE_B RATE_C]) = trans3(:,[35,36,37]);
    xfmr3(idx2, [RATE_A RATE_B RATE_C]) = trans3(:,[51,52,53]);
    xfmr3(idx3, [RATE_A RATE_B RATE_C]) = trans3(:,[67,68,69]);
    xfmr3(:, ANGMIN) = -360;        %% angle limits
    xfmr3(:, ANGMAX) =  360;
    %% winding status
    xfmr3(:, BR_STATUS) = 1;        %% initialize to all in-service
    status = trans3(:, 12);
    k1 = find(status == 0 | status == 4);   %% winding 1 out-of-service
    k2 = find(status == 0 | status == 2);   %% winding 2 out-of-service
    k3 = find(status == 0 | status == 3);   %% winding 3 out-of-service
    if ~isempty(k1)
        xfmr3(idx1(k1), BR_STATUS) = 0;
    end
    if ~isempty(k2)
        xfmr3(idx2(k2), BR_STATUS) = 0;
    end
    if ~isempty(k3)
        xfmr3(idx3(k3), BR_STATUS) = 0;
    end
end

%% combine 2-winding and 3-winding transformer data
xfmr = [xfmr2; xfmr3];

% %% delete out-of-service windings
% k = find(xfmr(:, BR_STATUS) == 0);
% xfmr(k, :) = [];
% k = find(trans3(:,12) <= 0);
% starbus(k, :) = [];

%% finish adding the star point bus
bus = [bus; starbus];
if nt3 > 0
    warns{end+1} = sprintf('Added buses %d-%d as star-points for 3-winding transformers.', ...
        starbus(1, BUS_I), starbus(end, BUS_I));
    if verbose
        fprintf('Added buses %d-%d as star-points for 3-winding transformers.\n', ...
            starbus(1, BUS_I), starbus(end, BUS_I));
    end
end
if nargin > 6 && nargout > 3
    starbus_name = cell(nt3, 1);
    for k = 1:nt3
        starbus_name{k} = sprintf('STR_PT_XF_%-2d', k);
    end
    bus_name = [bus_name; starbus_name];
end
