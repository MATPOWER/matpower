function fname_out = save2psse(fname, mpc, rawver)
%SAVE2PSSE  Saves a MATPOWER case to PSS/E RAW format.
%   SAVE2PSSE(FNAME, MPC)
%
%   FNAME = SAVE2PSSE(FNAME, ...)
%
%   Saves a MATPOWER case struct MPC as a PSS/E RAW file. The FNAME
%   parameter is a string containing the name of the file to be created
%   or overwritten. If FNAME does not include a file extension, '.raw'
%   will be added. Optionally returns the, possibly updated, filename.
%   Currently exports to RAW format Rev 33.

%   To do:
%   SAVE2PSSE(FNAME, CASESTRUCT, VERSION) (not yet implmented)

%   MATPOWER
%   Copyright (c) 2017 Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
c = idx_dcline;

%% verify valid filename
[pathstr, fcn_name, extension] = fileparts(fname);
if isempty(extension)
    extension = '.raw';
end
fname = fullfile(pathstr, [fcn_name extension]);

%% open file
[fd, msg] = fopen(fname, 'wt');     %% print it to a text file
if fd == -1
    error(['save2psse: ', msg]);
end

%% create map of external bus numbers to bus indices
i2e = mpc.bus(:, BUS_I);
e2i = sparse(max(i2e), 1);
e2i(i2e) = (1:size(mpc.bus, 1))';

%% Case Identification Data
fprintf(fd, '0, %d, 33, 0, 0, 60  / %s - MATPOWER %s\n', mpc.baseMVA, datestr(clock), mpver);
fprintf(fd, '\n');
fprintf(fd, '\n');

%% Bus Data
nb = size(mpc.bus, 1);  %% number of buses
if isfield(mpc, 'bus_name')
    maxlen = max(cellfun(@length, mpc.bus_name));
    tmp = cell2mat(cellfun( @(s)sprintf('%-*s', max(maxlen, 12), s), ...
                            mpc.bus_name, 'UniformOutput', 0 ));
    bn = tmp(:, 1:12);
else
    bn = cell2mat(cellfun( @(s)sprintf('BUS %-8d', s), ...
                           num2cell(mpc.bus(:, BUS_I)), 'UniformOutput', 0 ));
end
%%             I,                       'NAME',BASEKV,IDE,AREA,ZONE,OWNER,   VM,    VA,NVHI, NVLO, EVHI, EVLO
fprintf(fd, '%6d, ''%c%c%c%c%c%c%c%c%c%c%c%c'', %9.7g, %d, %4d, %4d, %d, %11.9g, %11.9g, %4g, %4g, %4g, %4g\n', ...
    [ mpc.bus(:, BUS_I) double(bn) ...
      mpc.bus(:, [BASE_KV BUS_TYPE BUS_AREA ZONE]) ones(nb, 1) ...
      mpc.bus(:, [VM VA VMAX VMIN VMAX VMIN]) ...
    ]');
fprintf(fd, '0 / END OF BUS DATA, BEGIN LOAD DATA\n');

%% Load Data
ild = find(mpc.bus(:, PD) ~= 0 | mpc.bus(:, PD) ~= 0);  %% bus indices of fixed loads
nld = length(ild);
if nld
    %%             I, ID,STATUS,AREA,ZONE, PL,   QL,   IP,   IQ,   YP,   YQ,OWNER, SCALE, INTRPT
    fprintf(fd, '%6d, %2d, %d, %4d, %4d, %9.7g, %9.7g, %.9g, %.9g, %.9g, %.9g, %d, %d, %d\n', ...
        [ mpc.bus(ild, BUS_I) ones(nld, 2) ...
          mpc.bus(ild, [BUS_AREA ZONE PD QD]) ...
          zeros(nld, 4) ones(nld, 2) zeros(nld, 1) ...
        ]');
end
idl = find(isload(mpc.gen));        %% gen indices of dispatchable loads
ndl = length(idl);
ib = e2i(mpc.gen(idl, GEN_BUS));    %% bus indices corresponding to dispatchable loads
if ndl
    %%             I, ID,STATUS,AREA,ZONE, PL,   QL,   IP,   IQ,   YP,   YQ,OWNER, SCALE, INTRPT
    fprintf(fd, '%6d, %2d, %d, %4d, %4d, %9.7g, %9.7g, %.9g, %.9g, %.9g, %.9g, %d, %d, %d\n', ...
        [ mpc.bus(ib, BUS_I) ones(nld, 2) mpc.bus(ib, [BUS_AREA ZONE]) ...
          -mpc.gen(idl, [PMIN QMIN]) ...
          zeros(nld, 4) ones(nld, 3) ...
        ]');
end
%%----- TODO: add correct ID value for multiple loads at a bus
fprintf(fd, '0 / END OF LOAD DATA, BEGIN FIXED SHUNT DATA\n');

%% Fixed Bus Shunt Data
ifs = find(mpc.bus(:, GS) ~= 0 | mpc.bus(:, BS) ~= 0);
nfs = length(ifs);
if nfs
    %%             I, ID,STATUS,GL, BL
    fprintf(fd, '%6d, %d, %d, %7g, %7g\n', ...
        [ mpc.bus(ifs, BUS_I) ones(nfs, 2) mpc.bus(ifs, [GS BS]) ]');
end
fprintf(fd, '0 / END OF FIXED SHUNT DATA, BEGIN GENERATOR DATA\n');

%% Generator Data
ig = find(~isload(mpc.gen));    %% gen indices (of real generators)
ng = length(ig);
wind = zeros(ng, 1);
if isfield(mpc, 'gentype')
    iw = find(cellfun(@(s)strcmp(s,'wind'), mpc.genfuel));
    wind(iw) = 1;
end
if isfield(mpc, 'genfuel')
    iw = find(cellfun(@(s)strcmp(s(1),'W'), mpc.gentype));
    wind(iw) = 1;
end
if ng
    %%             I,  ID,    PG,    QG,    QT,    QB,    VS,IREG,MBASE,ZR,ZX, RT, XT,GTAP,STAT,RMPCT,PT,PB,O1,F1,...,O4,F4,WMOD,WPF
    fprintf(fd, '%6d, %2d, %9.7g, %9.7g, %9.7g, %9.7g, %8.7g, %d, %7g, %g, %g, %g, %g, %g, %d, %g, %9.7g, %9.7g, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n', ...
        [ mpc.gen(ig, GEN_BUS) ones(ng, 1) mpc.gen(ig, [PG QG QMAX QMIN VG]) ...
          zeros(ng, 1) mpc.gen(ig, MBASE) zeros(ng, 1) ones(ng, 1) ...
          zeros(ng, 2) ones(ng, 1) (mpc.gen(ig, GEN_STATUS) > 0) ...
          100*ones(ng, 1) mpc.gen(ig, [PMAX PMIN]) ones(ng, 2) ...
          zeros(ng, 1) ones(ng, 1) zeros(ng, 1) ones(ng, 1) ...
          zeros(ng, 1) ones(ng, 1) wind ones(ng, 1) ...
        ]');
end
%%----- TODO: add correct ID value for multiple generators at a bus
fprintf(fd, '0 / END OF GENERATOR DATA, BEGIN BRANCH DATA\n');

%% Non-Transformer Branch Data
il = find(mpc.branch(:, TAP) == 0 & mpc.branch(:, SHIFT) == 0);    %% branch indices (of non-transformers)
nl = length(il);
if nl
    %%             I,   J,CKT,     R,     X,     B,RATEA,RATEB,RATEC,GI,BI,GJ, BJ, ST,MET,LEN, O1, F1,...,             O4, F4
    fprintf(fd, '%6d, %6d, %d, %8.7g, %8.7g, %8.7g, %7g, %7g, %7g, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n', ...
        [ mpc.branch(il, [F_BUS T_BUS]) ones(nl, 1) ...
          mpc.branch(il, [BR_R BR_X BR_B RATE_A RATE_B RATE_C]) ...
          zeros(nl, 4) mpc.branch(il, [BR_STATUS]) ones(nl, 1) ...
          zeros(nl, 1) ones(nl, 2) zeros(nl, 1) ones(nl, 1) ...
          zeros(nl, 1) ones(nl, 1) zeros(nl, 1) ones(nl, 1) ...
        ]');
end
fprintf(fd, '0 / END OF BRANCH DATA, BEGIN TRANSFORMER DATA\n');

%% Transformer Data
it = find(mpc.branch(:, TAP) ~= 0 | mpc.branch(:, SHIFT) ~= 0);    %% branch indices (of transformers)
nt = length(it);
if nt
    %%                                                                                                                          R1-2,X1-2,SBASE1-2                                                          WINDV2,NOMV2
    %%             I,    J, K,CKT, CW, CZ, CM,MAG1,MAG2,NMETR,      'NAME',STAT, O1, F1,...,             O4, F4,          VECGRP            WINDV1,NOMV1,ANG1,RATA1,RATB1,RATC1,COD1,CONT1,RMA1,RMI1,VMA1,VMI1,NTP1,TAB1,CR1,CX1,CNXA1\nWINDV2,NOMV2
    fprintf(fd, '%6d, %6d, %d, %d, %d, %d, %d, %d, %d, %d, ''            '', %d, %d, %d, %d, %d, %d, %d, %d, %d, ''            ''\n%g, %g, %g\n%5g, %d, %5g, %7g, %7g, %7g, %d, %d, %g, %g, %d, %d, %g, %g, %g\n%5g, %g\n', ...
        [ mpc.branch(it, [F_BUS T_BUS]) zeros(nt, 1) ones(nt, 4) ...
          zeros(nt, 2) 2*ones(nt, 1) mpc.branch(it, BR_STATUS) ones(nt, 2) ...
          zeros(nt, 1) ones(nt, 1) zeros(nt, 1) ones(nt, 1) ...
          zeros(nt, 1) ones(nt, 1) ...
          mpc.branch(it, [BR_R BR_X]) mpc.baseMVA*ones(nt, 1) ...
          mpc.branch(it, TAP) zeros(nt, 1) mpc.branch(it, SHIFT) ...
          mpc.branch(it, [RATE_A RATE_B RATE_C]) zeros(nt, 2) ...
          1.1*ones(nt, 1) 0.9*ones(nt, 1) 33*ones(nt, 1) zeros(nt, 4) ...
          ones(nt, 1) zeros(nt, 1) ...
        ]');
end
fprintf(fd, '0 / END OF TRANSFORMER DATA, BEGIN AREA DATA\n');

%% Area Interchange Data
%% I, ISW, PDES, PTOL, 'ARNAME'
fprintf(fd, '0 / END OF AREA DATA, BEGIN TWO-TERMINAL DC DATA\n');

%% Two-Terminal DC Transmission Line Data
if isfield(mpc, 'dcline')
    ndc = size(mpc.dcline, 1);
    if ndc
        dn = cell2mat(cellfun(  @(s)sprintf('DCLINE %-5d', s), ...
                                num2cell((1:ndc)'), 'UniformOutput', 0 ));
        %%                                                                                             IPR,NBR,ANMXR,ANMNR,RCR,XCR,EBASR,TRR,TAPR,TMXR,TMNR,    STPR,ICR,IFR,ITR,IDR,XCAPR
        %%                                'NAME',MDC, RDC, SETVL,VSCHD,VCMOD,RCOMP,DELTI,METER,DCVMIN,CCCITMX,CCCACC                                                                       IPI,NBI,ANMXI,ANMNI,RCI,XCI,EBASI,TRI,TAPI,TMXI,TMNI,STPI,ICI,IFI,ITI,IDI,XCAPI
        fprintf(fd, '''%c%c%c%c%c%c%c%c%c%c%c%c'', 1, 0.0, %9.7g, 500, 0.0, 0.0, 0.0, I, 0.0, 20, 1.0\n%6d, 1, 25.0, 15.0, 0, 0, %8.7g, 1.0, 1.0, 1.5, 0.51, 0.00625, 0, 0, 0, ''1'', 0.0\n%6d, 1, 25.0, 15.0, 0, 0, %8.7g, 1.0, 1.0, 1.5, 0.51, 0.00625, 0, 0, 0, ''1'', 0.0\n', ...
            [ double(dn) mpc.dcline(:, c.PF) ...
              mpc.dcline(:, c.F_BUS) ...
              mpc.bus(e2i(mpc.dcline(:, c.F_BUS)), BASE_KV) ...
              mpc.dcline(:, c.T_BUS) ...
              mpc.bus(e2i(mpc.dcline(:, c.T_BUS)), BASE_KV) ...
            ]');
    end
end
%%----- TODO: attempt some values to approximate losses
fprintf(fd, '0 / END OF TWO-TERMINAL DC DATA, BEGIN VOLTAGE SOURCE CONVERTER DATA\n');

%% Voltage Source Converter (VSC) DC Transmission Line Data
%% 'NAME', MDC, RDC, O1, F1, ... O4, F4
%% IBUS,TYPE,MODE,DCSET,ACSET,ALOSS,BLOSS,MINLOSS,SMAX,IMAX,PWF,MAXQ,MINQ,REMOT,RMPCT
%% IBUS,TYPE,MODE,DCSET,ACSET,ALOSS,BLOSS,MINLOSS,SMAX,IMAX,PWF,MAXQ,MINQ,REMOT,RMPCT
fprintf(fd, '0 / END OF VOLTAGE SOURCE CONVERTER DATA, BEGIN IMPEDANCE CORRECTION DATA\n');

%% Transformer Impedance Correction Tables
%% I, T1, F1, T2, F2, T3, F3, ... T11, F11
fprintf(fd, '0 / END OF IMPEDANCE CORRECTION DATA, BEGIN MULTI-TERMINAL DC DATA\n');

%% Multi-Terminal DC Transmission Line Data
%% 'NAME', NCONV, NDCBS, NDCLN, MDC, VCONV, VCMOD, VCONVN
%% IB,N,ANGMX,ANGMN,RC,XC,EBAS,TR,TAP,TPMX,TPMN,TSTP,SETVL,DCPF,MARG,CNVCOD
%% IDC, IB, AREA, ZONE, 'DCNAME', IDC2, RGRND, OWNER
%% IDC, JDC, DCCKT, MET, RDC, LDC
fprintf(fd, '0 / END OF MULTI-TERMINAL DC DATA, BEGIN MULTI-SECTION LINE DATA\n');

%% Multi-Section Line Grouping Data
%% I, J, ID, MET, DUM1, DUM2, ... DUM9
fprintf(fd, '0 / END OF MULTI-SECTION LINE DATA, BEGIN ZONE DATA\n');

%% Zone Data
%% I, 'ZONAME'
fprintf(fd, '0 / END OF ZONE DATA, BEGIN INTER-AREA TRANSFER DATA\n');

%% Interarea Transfer Data
%% ARFROM, ARTO, TRID, PTRAN
fprintf(fd, '0 / END OF INTER-AREA TRANSFER DATA, BEGIN OWNER DATA\n');

%% Owner Data
%% I, 'OWNAME'
fprintf(fd, '0 / END OF OWNER DATA, BEGIN FACTS CONTROL DEVICE DATA\n');

%% FACTS Device Data
%% 'NAME',I,J,MODE,PDES,QDES,VSET,SHMX,TRMX,VTMN,VTMX,VSMX,IMX,LINX,RMPCT,OWNER,SET1,SET2,VSREF,REMOT,'MNAME'
fprintf(fd, '0 / END OF FACTS CONTROL DEVICE DATA, BEGIN SWITCHED SHUNT DATA\n');

%% Switched Shunt Data
%% I, MODSW, ADJM, STAT, VSWHI, VSWLO, SWREM, RMPCT, 'RMIDNT', BINIT, N1, B1, N2, B2, ... N8, B8
fprintf(fd, '0 / END OF SWITCHED SHUNT DATA, BEGIN GNE DEVICE DATA\n');

%% GNE Device Data
fprintf(fd, '0 / END OF GNE DEVICE DATA, BEGIN INDUCTION MACHINE DATA\n');

%% Induction Machine Data
fprintf(fd, '0 / END OF INDUCTION MACHINE DATA\n');

%% End of Data Indicator
fprintf(fd, 'Q\n');

%% close file
if fd ~= 1
    fclose(fd);
end

if nargout > 0
    fname_out = fname;
end
