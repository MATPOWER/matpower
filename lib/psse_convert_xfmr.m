function [xfmr, bus, bus_name] = psse_convert_xfmr(trans2, trans3, baseMVA, bus, bus_name)
%PSSE_CONVERT_XFMR Convert transformer data from PSS/E RAW to MATPOWER
%   [XFMR, BUS] = PSSE_CONVERT_XFMR(TRANS2, TRANS3, BASEMVA, BUS)
%   [XFMR, BUS, BUS_NAME] = PSSE_CONVERT_XFMR(TRANS2, TRANS3, BASEMVA, ...
%                                                       BUS, BUS_NAME)
%
%   Convert all transformer data read from a PSS/E RAW data file
%   into MATPOWER format. Returns a branch matrix corresponding to
%   the transformers and an updated bus matrix, with additional buses
%   added for the star points of three winding transformers.
%
%   Inputs:
%       TRANS2  : matrix of raw two winding transformer data returned
%                 by PSSE_READ in data.trans2.num
%       TRANS3  : matrix of raw three winding transformer data returned
%                 by PSSE_READ in data.trans3.num
%       BASEMVA : system MVA base
%       BUS     : MATPOWER bus matrix
%       BUS_NAME: (optional) cell array of bus names
%
%   Outputs:
%       XFMR    : MATPOWER branch matrix of transformer data
%       BUS     : updated MATPOWER bus matrix, with additional buses
%                 added for star points of three winding transformers
%       BUS_NAME: (optional) updated cell array of bus names
%
%   See also PSSE_CONVERT.

%   MATPOWER
%   $Id$
%   by Yujia Zhu, PSERC ASU
%   and Ray Zimmerman, PSERC Cornell
%   Based on mptransin.m and mptransficbus.m, written by:
%       Yujia Zhu, Jan 2014, yzhu54@asu.edu.
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

%% create fictitious buses for star point of 3 winding transformers
nb = size(bus, 1);
nt3 = size(trans3, 1);
starbus = zeros(nt3, VMIN);     %% initialize output matrix
mb =  max(bus(:, BUS_I));       %% find maximum bus number
%% for area and zone data, use the data for winding 1
wind1bus = trans3(:,1);         %% winding1 bus number
ib = interp1q(bus(:,1), (1:nb)', wind1bus); %% corresponding bus index
starbus(:, BUS_I) = (mb+1:mb+nt3)'; %% bus numbers follow originals
starbus(:, BUS_TYPE) = PQ;
starbus(trans3(:,12)==0, BUS_TYPE) = NONE;  %% isolated if transformer is off-line
starbus(:, VA) = trans3(:,31);  %% VA = star point phase angle
starbus(:, VM) = trans3(:,30);  %% VM = star point voltage magnitude (PU)
starbus(:, [BUS_AREA, ZONE]) = bus(ib, [7,11]); %% wind1 bus area, zone
starbus(:, BASE_KV) = 1;        %% baseKV = 1 kV
starbus(:, VMAX) = 1.1;
starbus(:, VMIN) = 0.9;
bus = [bus; starbus];
if nargin > 4 && nargout > 2
    starbus_name = cell(nt3, 1);
    for k = 1:nt3
        starbus_name{k} = sprintf('STAR_POINT_XFMR_%d', k);
    end
    bus_name = [bus_name; starbus_name];
end

%% eliminate out-of-service transformers (but, why?)
status2 = trans2(:,12); % find the status of the 2 winding transformers
status3 = trans3(:,12); % find the status of the 3 winding transformers
idx2 = find(status2 > 0);
idx3 = find(status3 > 0); %find the transformers in service
trans2 = trans2(idx2,:);
trans3 = trans3(idx3,:);
starbus = starbus(idx3,:);

nt2 = size(trans2, 1);
nt3 = size(trans3, 1);

%% PSS/E two winding transformer data
%% I,J,K,CKT,CW,CZ,CM,MAG1,MAG2,NMETR,’NAME’,STAT,O1,F1,...,O4,F4,VECGRP (1-21)
%% R1-2,X1-2,SBASE1-2 (22-24)
%% WINDV1,NOMV1,ANG1,RATA1,RATB1,RATC1,COD1,CONT1,RMA1,RMI1,VMA1,VMI1,NTP1,TAB1,CR1,CX1,CNXA1
%% WINDV2,NOMV2
[tf,fbus] = ismember(trans2(:,1), bus(:, BUS_I));      %% I
[tf,tbus] = ismember(trans2(:,2), bus(:, BUS_I));      %% J
Zbs = bus(:, BASE_KV).^2 / baseMVA;     %% system impedance base
if ~isempty(trans2)
    cw2 = find(trans2(:,5) == 2);   %% CW = 2
    cw3 = find(trans2(:,5) == 3);   %% CW = 3
    cw23 = [cw2;cw3];               %% CW = 2 or 3
    cz2 = find(trans2(:,6) == 2);   %% CZ = 2
    cz3 = find(trans2(:,6) == 3);   %% CZ = 3
    cz23 = [cz2;cz3];               %% CZ = 2 or 3

    R = trans2(:,21);               %% 
    X = trans2(:,22);
    Zb = ones(nt2, 1);
    Zb(cz23) = trans2(cz23,25).^2 ./ trans2(cz23,23);
    R(cz2)  = R(cz2)  .* Zb(cz2)  ./ Zbs(fbus(cz2));
    X(cz23) = X(cz23) .* Zb(cz23) ./ Zbs(fbus(cz23));
    R(cz3)  = trans2(cz3,25).^2 ./ trans2(cz3,21) ./ Zbs(fbus(cz3));
    tap = trans2(:,24) ./ trans2(:,40);
    tap(cw23) = tap(cw23) .* bus(tbus(cw23), BASE_KV)./bus(fbus(cw23), BASE_KV);
    tap(cw3) = tap(cw3) .* trans2(cw3,25)./trans2(cw3,41);
    shift = trans2(:, 26);
else
    disp('No two winding transformer in the network');
end

%% PSS/E three winding transformer data
%% I,J,K,CKT,CW,CZ,CM,MAG1,MAG2,NMETR,’NAME’,STAT,O1,F1,...,O4,F4,VECGRP
%% R1-2,X1-2,SBASE1-2,R2-3,X2-3,SBASE2-3,R3-1,X3-1,SBASE3-1,VMSTAR,ANSTAR
%% WINDV1,NOMV1,ANG1,RATA1,RATB1,RATC1,COD1,CONT1,RMA1,RMI1,VMA1,VMI1,NTP1,TAB1,CR1,CX1,CNXA1
%% WINDV2,NOMV2,ANG2,RATA2,RATB2,RATC2,COD2,CONT2,RMA2,RMI2,VMA2,VMI2,NTP2,TAB2,CR2,CX2,CNXA2
%% WINDV3,NOMV3,ANG3,RATA3,RATB3,RATC3,COD3,CONT3,RMA3,RMI3,VMA3,VMI3,NTP3,TAB3,CR3,CX3,CNXA3
[tf,ind1] = ismember(trans3(:,1), bus(:, BUS_I));
[tf,ind2] = ismember(trans3(:,2), bus(:, BUS_I));
[tf,ind3] = ismember(trans3(:,3), bus(:, BUS_I));
% Each three winding transformer will be converted into 3 branches:
% The branches will be in the order of
% # winding1 -> # winding2
% # winding2 -> # winding3
% # winding3 -> # winding1
if ~isempty(trans3)
    cw2 = find(trans3(:,5) == 2);   %% CW = 2
    cw3 = find(trans3(:,5) == 3);   %% CW = 3
    cw23 = [cw2;cw3];               %% CW = 2 or 3
    cz2 = find(trans3(:,6) == 2);   %% CZ = 2
    cz3 = find(trans3(:,6) == 3);   %% CZ = 3
    cz23 = [cz2;cz3];               %% CZ = 2 or 3

    tap1 = trans3(:, 32);   %% off nominal tap ratio of branch 1
    tap2 = trans3(:, 48);   %% off nominal tap ratio of branch 2
    tap3 = trans3(:, 64);   %% off nominal tap ratio of branch 3
    tap1(cw23) = tap1(cw23) ./ bus(ind1(cw23), BASE_KV);
    tap2(cw23) = tap2(cw23) ./ bus(ind2(cw23), BASE_KV);
    tap3(cw23) = tap3(cw23) ./ bus(ind3(cw23), BASE_KV);
    tap1(cw3)  = tap1(cw3)  .* trans3(cw3, 33);
    tap2(cw3)  = tap2(cw3)  .* trans3(cw3, 49);
    tap3(cw3)  = tap3(cw3)  .* trans3(cw3, 65);
    shift1 = trans3(:, 34);
    shift2 = trans3(:, 50);
    shift3 = trans3(:, 66);

    R12 = trans3(:, 21);
    X12 = trans3(:, 22);
    R23 = trans3(:, 24);
    X23 = trans3(:, 25);
    R31 = trans3(:, 27);
    X31 = trans3(:, 28);
    Zb1 = trans3(:, 33).^2 ./ trans3(:, 23);
    Zb2 = trans3(:, 49).^2 ./ trans3(:, 26);
    Zb3 = trans3(:, 65).^2 ./ trans3(:, 29);

    R12(cz2) = R12(cz2) .* Zb1(cz2) ./ Zbs(ind1(cz2));
    X12(cz2) = X12(cz2) .* Zb1(cz2) ./ Zbs(ind1(cz2));
    R23(cz2) = R23(cz2) .* Zb2(cz2) ./ Zbs(ind2(cz2));
    X23(cz2) = X23(cz2) .* Zb2(cz2) ./ Zbs(ind2(cz2));
    R31(cz2) = R31(cz2) .* Zb3(cz2) ./ Zbs(ind3(cz2));
    X31(cz2) = X31(cz2) .* Zb3(cz2) ./ Zbs(ind3(cz2));

    R12(cz3) = trans3(cz3,33) .^ 2 ./ trans3(cz3,21) ./ Zbs(ind1(cz3));
    X12(cz3) = trans3(cz3,22) .* Zb1(cz3) ./ Zbs(ind3(cz3));
    R23(cz3) = trans3(cz3,49) .^ 2 ./ trans3(cz3,24) ./ Zbs(ind2(cz3));
    X23(cz3) = trans3(cz3,25) .* Zb2(cz3) ./ Zbs(ind3(cz3));
    R31(cz3) = trans3(cz3,65) .^ 2 ./ trans3(cz3,27) ./ Zbs(ind3(cz3));
    X31(cz3) = trans3(cz3,28) .* Zb3(cz3) ./ Zbs(ind3(cz3));

    R1 = (R12+R31-R23) ./ 2;
    R2 = (R12+R23-R31) ./ 2;
    R3 = (R31+R23-R12) ./ 2;
    X1 = (X12+X31-X23) ./ 2;
    X2 = (X12+X23-X31) ./ 2;
    X3 = (X31+X23-X12) ./ 2;
else
    disp('No three winding transformer in the network');
end

%% Put the the transformer data together based on MATPOWER data format

%%% two winding transformers %%%
xfmr2 = zeros(nt2, ANGMAX); % Initialize the info matrix of two winding transformer
xfmr2(:, [F_BUS T_BUS]) = trans2(:,[1,2]);
xfmr2(:, [BR_R BR_X]) = [R,X];
xfmr2(:, [RATE_A RATE_B RATE_C]) = trans2(:,[27,28,29]);
xfmr2(:, [TAP SHIFT]) = [tap,shift];
xfmr2(:, BR_STATUS) = trans2(:,12);
xfmr2(:, ANGMIN) = -360;
xfmr2(:, ANGMAX) = 360;

%%% three winding transformers %%%
xfmr3 = zeros(3*nt3, ANGMAX);
del_idx3 = [];
status3 = trans3(:, 12);
for i = 1:nt3
    % check the status vector in RAW data and find the windings that out of
    % services. Branches connecting the out of service winding will be
    % eliminated
    % to determine the branches that will be eliminated
    status = status3(i);
    if status == 2
        del_idx3 = [del_idx3;(i-1)*3+2]; % record the index of winding 2
    elseif status == 3
        del_idx3 = [del_idx3;(i-1)*3+3]; % record the index of winding 3
    elseif status == 4
        del_idx3 = [del_idx3;(i-1)*3+1]; % record the index of winding 1
    end
    % begin to assign values
    % numbering the winding buses
    xfmr3(((i-1)*3+1),[F_BUS T_BUS]) = [trans3(i,1),starbus(i,1)]; % 
    xfmr3(((i-1)*3+2),[F_BUS T_BUS]) = [trans3(i,2),starbus(i,1)];
    xfmr3(((i-1)*3+3),[F_BUS T_BUS]) = [trans3(i,3),starbus(i,1)];
    % impedances
    xfmr3(((i-1)*3+1),[BR_R BR_X]) = [R1(i),X1(i)];
    xfmr3(((i-1)*3+2),[BR_R BR_X]) = [R2(i),X2(i)];
    xfmr3(((i-1)*3+3),[BR_R BR_X]) = [R3(i),X3(i)];
    % ratings
    xfmr3(((i-1)*3+1),[RATE_A RATE_B RATE_C]) = trans3(i,[35,36,37]);
    xfmr3(((i-1)*3+2),[RATE_A RATE_B RATE_C]) = trans3(i,[51,52,53]);
    xfmr3(((i-1)*3+3),[RATE_A RATE_B RATE_C]) = trans3(i,[67,68,69]);
    % tap ratio & phase shifts
    xfmr3(((i-1)*3+1),[TAP SHIFT]) = [tap1(i),shift1(i)];
    xfmr3(((i-1)*3+2),[TAP SHIFT]) = [tap2(i),shift2(i)];
    xfmr3(((i-1)*3+3),[TAP SHIFT]) = [tap3(i),shift3(i)];
end
xfmr3(del_idx3,:) = []; % delete the off line windings
xfmr3(:, ANGMIN) = -360*ones(size(xfmr3,1),1); % Angle limits
xfmr3(:, ANGMAX) = 360*ones(size(xfmr3,1),1); % Adding the angle limit on the transformer data
xfmr3(:, BR_STATUS) = ones(size(xfmr3,1),1); % only in-service windings will be retained
%% integrate all transformer data
xfmr = [xfmr2; xfmr3];
