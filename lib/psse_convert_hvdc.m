function dcline = psse_convert_hvdc(dc, bus)
%PSSE_CONVERT_HVDC Convert HVDC data from PSS/E RAW to MATPOWER
%   DCLINE = PSSE_CONVERT_HVDC(DC, BUS)
%
%   Convert all two terminal HVDC line data read from a PSS/E
%   RAW data file into MATPOWER format. Returns a dcline matrix for
%   inclusion in a MATPOWER case struct.
%
%   Inputs:
%       DC  : matrix of raw two terminal HVDC line data returned by
%             PSSE_READ in data.twodc.num
%       BUS : MATPOWER bus matrix
%
%   Output:
%       DCLINE : a MATPOWER dcline matrix suitable for inclusion in
%                a MATPOWER case struct.
%
%   See also PSSE_CONVERT.

%   MATPOWER
%   Copyright (c) 2014-2016, Power Systems Engineering Research Center (PSERC)
%   by Yujia Zhu, PSERC ASU
%   and Ray Zimmerman, PSERC Cornell
%   Based on mpdcin.m and mpqhvdccal.m, written by:
%       Yujia Zhu, Jan 2014, yzhu54@asu.edu.
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
c = idx_dcline;

nb = size(bus, 1);
ndc = size(dc, 1);
e2i = sparse(bus(:, BUS_I), ones(nb, 1), 1:nb, max(bus(:, BUS_I)), 1);
if ~ndc
    dcline = [];
    return;
end

%% extract data
MDC = dc(:,2); % Control mode
SETVL = dc(:,4); % depend on control mode: current or power demand
VSCHD = dc(:,5); % scheduled compounded dc voltage
ANMXR = dc(:,15); % nominal maximum rectifier firing angle
ANMNR = dc(:,16); % nominal minimum rectifier firing angle
GAMMX = dc(:,32); % nominal maximum inverter firing angle
GAMMN = dc(:,33); % nominal minimum inverter firing angle
SETVL = abs(SETVL);
% Convert the voltage on rectifier side and inverter side
% The value is calculated as basekV/VSCHD
% basekV is the bus base voltage, VSCHD is the scheduled compounded
% voltage
dcline = zeros(ndc, c.LOSS1); % initiate the hvdc data format
indr = dc(:,13); % rectifier end bus number
indi = dc(:,30); % inverter end bus number
dcind = [indr indi]; 
% bus nominal voltage
Vr = bus(e2i(indr), VM);
Vi = bus(e2i(indi), VM);
%% Calculate the real power input at the from end
PMW = zeros(ndc, 1);
for i = 1:ndc
    if MDC(i) == 1
        PMW(i) = SETVL(i); % SETVL is the desired real power demand
    elseif MDC(i) == 2;
        PMW(i) = SETVL(i)*VSCHD(i)/1000; % SETVL is the current in amps (need devide 1000 to convert to MW)
    else PMW(i) = 0;
    end
end
%% calculate reactive power limits
[Qrmin,Qrmax] = psse_convert_hvdc_Qlims(ANMXR,ANMNR,PMW);    %% rectifier end
[Qimin,Qimax] = psse_convert_hvdc_Qlims(GAMMX,GAMMN,PMW);    %% inverter end
%% calculate the loss coefficient (Only consider the l1)
% l1 = P'.*RDC;

%% conclude all info
status = ones(ndc, 1);
status(MDC==0) = 0;     %% set status of blocked HVDC lines to zero
% dcline(:,[1 2 3 4 5 8 9 10 11 12 13 14 15]) = [indr,indi,status,PMW, PMW, Vr, Vi,0.85*PMW, 1.15*PMW, Qrmin, Qrmax, Qimin, Qimax];
dcline(:, [c.F_BUS c.T_BUS c.BR_STATUS c.PF c.PT c.VF c.VT ...
            c.PMIN c.PMAX c.QMINF c.QMAXF c.QMINT c.QMAXT]) = ...
    [indr indi status PMW PMW Vr Vi 0.85*PMW 1.15*PMW Qrmin Qrmax Qimin Qimax];


function [Qmin, Qmax] = psse_convert_hvdc_Qlims(alphamax,alphamin,P)
%PSSE_CONVERT_HVDC_QLIMS calculate HVDC line reactive power limits
%
%   [Qmin, Qmax] = psse_convert_hvdc_Qlims(alphamax,alphamin,P)
%
% Inputs:
%       alphamax :  maximum firing angle
%       alphamin :  minimum steady-state rectifier firing angle
%       P :         real power demand
% Outputs:
%       Qmin :  lower limit of reactive power
%       Qmax :  upper limit of reactive power 
%
% Note:
%   This function calculates the reactive power at the rectifier or inverter
%   end. It is assumed the maximum overlap angle is 60 degree (see
%   Kimbark's book). The maximum reactive power is calculated with the
%   power factor:
%       pf = acosd(0.5*(cosd(alphamax(i))+cosd(60))),
%   where, 60 is the maximum delta angle.

len = length(alphamax);
phi = zeros(size(alphamax));
Qmin = phi;
Qmax = phi;
for i = 1:len
    %% minimum reactive power calculated under assumption of no overlap angle
    %% i.e. power factor equals to tan(alpha)
    Qmin(i) = P(i)*tand(alphamin(i));

    %% maximum reactive power calculated when overlap angle reaches max
    %% value (60 deg). I.e.
    %%      cos(phi) = 1/2*(cos(alpha)+cos(delta))
    %%      Q = P*tan(phi)
    phi(i) = acosd(0.5*(cosd(alphamax(i))+cosd(60)));
    Qmax(i) = P(i)*tand(phi(i));
    if Qmin(i)<0
        Qmin(i) = -Qmin(i);
    end
    if Qmax(i)<0
        Qmax(i) = -Qmax(i);
    end
end
