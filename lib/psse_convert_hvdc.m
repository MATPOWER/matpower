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
%   $Id$
%   by Yujia Zhu, PSERC ASU
%   and Ray Zimmerman, PSERC Cornell
%   Based on mpdcin.m and mpqhvdccal.m, written by:
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
c = idx_dcline;

nb = size(bus, 1);
ndc = size(dc, 1);
if ~ndc
	dcline = [];
	return;
end

%% v1: (http://www.ee.washington.edu/research/pstca/formats/pti.txt)
%%  I,MDC,RDC,SETVL,VSCHD,VCMOD,RCOMP,DELTI,METER (1-9)
%%  IPR,NBR,ALFMAX,ALFMN,RCR,XCR,EBASR,TRR,TAPR,TPMXR,TPMNR,TSTPR (13-24)
%%  IPI,NBI,GAMMX,GAMMN,RCI,XCI,EBASI,TRI,TAPI,TPMXI,TPMNI,TSTPI (30-41)
%% v29-30:
%%  I,MDC,RDC,SETVL,VSCHD,VCMOD,RCOMP,DELTI,METER,DCVMIN,CCCITMX,CCCACC (1-12)
%%  IPR,NBR,ALFMX,ALFMN,RCR,XCR,EBASR,TRR,TAPR,TMXR,TMNR,STPR,ICR,IFR,ITR,IDR,XCAPR (13-29)
%%  IPI,NBI,GAMMX,GAMMN,RCI,XCI,EBASI,TRI,TAPI,TMXI,TMNI,STPI,ICI,IFI,ITI,IDI,XCAPI (30-46)
%% v31-33:
%%  'NAME',MDC,RDC,SETVL,VSCHD,VCMOD,RCOMP,DELTI,METER,DCVMIN,CCCITMX,CCCACC (1-12)
%%  IPR,NBR,ANMXR,ANMNR,RCR,XCR,EBASR,TRR,TAPR,TMXR,TMNR,STPR,ICR,IFR,ITR,IDR,XCAPR (13-29)
%%  IPI,NBI,ANMXI,ANMNI,RCI,XCI,EBASI,TRI,TAPI,TMXI,TMNI,STPI,ICI,IFI,ITI,IDI,XCAPI (30-46)

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
nindr = interp1q(bus(:, BUS_I),(1:nb)',indr);
nindi = interp1q(bus(:, BUS_I),(1:nb)',indi);
% bus nominal voltage
Vr = bus(nindr, VM);
Vi = bus(nindi, VM);
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

%% Conclude all info
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
