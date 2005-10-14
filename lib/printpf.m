function printpf(baseMVA, bus, gen, branch, f, success, et, fd, mpopt)
%PRINTPF   Prints power flow results.
%   printpf(baseMVA, bus, gen, branch, f, success, et, fd, mpopt) prints
%   powerflow results to fd (a file descriptor which defaults to STDOUT).
%   mpopt is a MATPOWER options vector (see 'help mpoption' for details).
%   Uses default options if this parameter is not given. The objective
%   function value is given in f and the elapsed time (seconds to compute
%   opf) in et.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialization -----
%% default arguments
if nargin < 9
    mpopt = mpoption;   %% use default options
    if nargin < 8
        fd = 1;         %% print to stdio by default
    end
end
if isempty(f)
    isOPF = 0;      %% have only simple PF data
else
    isOPF = 1;      %% have OPF data
end

%% options
dc              = mpopt(10);        %% use DC formulation?
OUT_ALL         = mpopt(32);
OUT_ANY         = OUT_ALL == 1;     %% set to true if any pretty output is to be generated
OUT_SYS_SUM     = OUT_ALL == 1 | (OUT_ALL == -1 & mpopt(33));
OUT_AREA_SUM    = OUT_ALL == 1 | (OUT_ALL == -1 & mpopt(34));
OUT_BUS         = OUT_ALL == 1 | (OUT_ALL == -1 & mpopt(35));
OUT_BRANCH      = OUT_ALL == 1 | (OUT_ALL == -1 & mpopt(36));
OUT_GEN         = OUT_ALL == 1 | (OUT_ALL == -1 & mpopt(37));
OUT_ANY         = OUT_ANY | (OUT_ALL == -1 & (OUT_SYS_SUM | OUT_AREA_SUM | OUT_BUS | OUT_BRANCH | OUT_GEN));
if OUT_ALL == -1
    OUT_ALL_LIM = mpopt(38);
elseif OUT_ALL == 1
    OUT_ALL_LIM = 2;
else
    OUT_ALL_LIM = 0;
end
OUT_ANY         = OUT_ANY | OUT_ALL_LIM >= 1;
if OUT_ALL_LIM == -1
    OUT_V_LIM       = mpopt(39);
    OUT_LINE_LIM    = mpopt(40);
    OUT_PG_LIM      = mpopt(41);
    OUT_QG_LIM      = mpopt(42);
else
    OUT_V_LIM       = OUT_ALL_LIM;
    OUT_LINE_LIM    = OUT_ALL_LIM;
    OUT_PG_LIM      = OUT_ALL_LIM;
    OUT_QG_LIM      = OUT_ALL_LIM;
end
OUT_ANY         = OUT_ANY | (OUT_ALL_LIM == -1 & (OUT_V_LIM | OUT_LINE_LIM | OUT_PG_LIM | OUT_QG_LIM));
OUT_RAW         = mpopt(43);

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% constant
j = sqrt(-1);

%% internal bus number
e2i = zeros(max(bus(:, BUS_I)), 1);     %% need internal bus numbering for a second
e2i(bus(:, BUS_I)) = [1:size(bus, 1)]';

%% sizes of things
nb = size(bus, 1);      %% number of buses
nl = size(branch, 1);   %% number of branches
ng = size(gen, 1);      %% number of generators

%% zero out some data to make printout consistent for DC case
if dc
    bus(:, [QD, BS])            = zeros(nb, 2);
    gen(:, [QG, QMAX, QMIN])    = zeros(ng, 3);
    branch(:, [BR_R, BR_B])     = zeros(nl, 2);
end

%% parameters
ties = find(bus(e2i(branch(:, F_BUS)), BUS_AREA) ~= bus(e2i(branch(:, T_BUS)), BUS_AREA));
                        %% area inter-ties
tap = ones(nl, 1);                              %% default tap ratio = 1 for lines
xfmr = find(branch(:, TAP));                    %% indices of transformers
tap(xfmr) = branch(xfmr, TAP);                  %% include transformer tap ratios
tap = tap .* exp(-j*pi/180 * branch(:, SHIFT)); %% add phase shifters
nzld = find(bus(:, PD) | bus(:, QD));
sorted_areas = sort(bus(:, BUS_AREA));
s_areas = sorted_areas([1; find(diff(sorted_areas))+1]);    %% area numbers
na = length(s_areas);                           %% number of areas
nzsh = find(bus(:, GS) | bus(:, BS));
allg = find( ~isload(gen) );
ong  = find( gen(:, GEN_STATUS) > 0 & ~isload(gen) );
onld = find( gen(:, GEN_STATUS) > 0 &  isload(gen) );
V = bus(:, VM) .* exp(sqrt(-1) * pi/180 * bus(:, VA));
out = find(branch(:, BR_STATUS) == 0);          %% out-of-service branches
nout = length(out);
if dc
    loss = zeros(nl,1);
else
    loss = baseMVA * abs(V(e2i(branch(:, F_BUS))) ./ tap - V(e2i(branch(:, T_BUS)))) .^ 2 ./ ...
                (branch(:, BR_R) - j * branch(:, BR_X));
end
fchg = abs(V(e2i(branch(:, F_BUS))) ./ tap) .^ 2 .* branch(:, BR_B) * baseMVA / 2;
tchg = abs(V(e2i(branch(:, T_BUS)))       ) .^ 2 .* branch(:, BR_B) * baseMVA / 2;
loss(out) = zeros(nout, 1);
fchg(out) = zeros(nout, 1);
tchg(out) = zeros(nout, 1);

%%----- print the stuff -----
if OUT_ANY
    %% convergence & elapsed time
    if success
        fprintf(fd, '\nConverged in %.2f seconds', et);
    else
        fprintf(fd, '\nDid not converge (%.2f seconds)\n', et);
    end
    
    %% objective function value
    if isOPF
        fprintf(fd, '\nObjective Function Value = %.2f $/hr', f);
    end
end
if OUT_SYS_SUM
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n|     System Summary                                                           |');
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n\nHow many?                How much?              P (MW)            Q (MVAr)');
    fprintf(fd, '\n---------------------    -------------------  -------------  -----------------');
    fprintf(fd, '\nBuses         %6d     Total Gen Capacity   %7.1f       %7.1f to %.1f', nb, sum(gen(allg, PMAX)), sum(gen(allg, QMIN)), sum(gen(allg, QMAX)));
    fprintf(fd, '\nGenerators     %5d     On-line Capacity     %7.1f       %7.1f to %.1f', length(allg), sum(gen(ong, PMAX)), sum(gen(ong, QMIN)), sum(gen(ong, QMAX)));
    fprintf(fd, '\nCommitted Gens %5d     Generation (actual)  %7.1f           %7.1f', length(ong), sum(gen(ong, PG)), sum(gen(ong, QG)));
    fprintf(fd, '\nLoads          %5d     Load                 %7.1f           %7.1f', length(nzld)+length(onld), sum(bus(nzld, PD))-sum(gen(onld, PG)), sum(bus(nzld, QD))-sum(gen(onld, QG)));
    fprintf(fd, '\n  Fixed        %5d       Fixed              %7.1f           %7.1f', length(nzld), sum(bus(nzld, PD)), sum(bus(nzld, QD)));
    fprintf(fd, '\n  Dispatchable %5d       Dispatchable       %7.1f of %-7.1f%7.1f', length(onld), -sum(gen(onld, PG)), -sum(gen(onld, PMIN)), -sum(gen(onld, QG)));
    fprintf(fd, '\nShunts         %5d     Shunt (inj)          %7.1f           %7.1f', length(nzsh), ...
        -sum(bus(nzsh, VM) .^ 2 .* bus(nzsh, GS)), sum(bus(nzsh, VM) .^ 2 .* bus(nzsh, BS)) );
    fprintf(fd, '\nBranches       %5d     Losses (I^2 * Z)     %8.2f          %8.2f', nl, sum(real(loss)), sum(imag(loss)) );
    fprintf(fd, '\nTransformers   %5d     Branch Charging (inj)     -            %7.1f', length(xfmr), sum(fchg) + sum(tchg) );
    fprintf(fd, '\nInter-ties     %5d     Total Inter-tie Flow %7.1f           %7.1f', length(ties), sum(abs(branch(ties, PF)-branch(ties, PT))) / 2, sum(abs(branch(ties, QF)-branch(ties, QT))) / 2);
    fprintf(fd, '\nAreas          %5d', length(s_areas));
    fprintf(fd, '\n');
    fprintf(fd, '\n                          Minimum                      Maximum');
    fprintf(fd, '\n                 -------------------------  --------------------------------');
    [minv, mini] = min(bus(:, VM));
    [maxv, maxi] = max(bus(:, VM));
    fprintf(fd, '\nVoltage Magnitude %7.3f p.u. @ bus %-4d     %7.3f p.u. @ bus %-4d', minv, bus(mini, BUS_I), maxv, bus(maxi, BUS_I));
    [minv, mini] = min(bus(:, VA));
    [maxv, maxi] = max(bus(:, VA));
    fprintf(fd, '\nVoltage Angle   %8.2f deg   @ bus %-4d   %8.2f deg   @ bus %-4d', minv, bus(mini, BUS_I), maxv, bus(maxi, BUS_I));
    if ~dc
        [maxv, maxi] = max(real(loss));
        fprintf(fd, '\nP Losses (I^2*R)             -              %8.2f MW    @ line %d-%d', maxv, branch(maxi, F_BUS), branch(maxi, T_BUS));
        [maxv, maxi] = max(imag(loss));
        fprintf(fd, '\nQ Losses (I^2*X)             -              %8.2f MVAr  @ line %d-%d', maxv, branch(maxi, F_BUS), branch(maxi, T_BUS));
    end
    if isOPF
        [minv, mini] = min(bus(:, LAM_P));
        [maxv, maxi] = max(bus(:, LAM_P));
        fprintf(fd, '\nLambda P        %8.2f $/MWh @ bus %-4d   %8.2f $/MWh @ bus %-4d', minv, bus(mini, BUS_I), maxv, bus(maxi, BUS_I));
        [minv, mini] = min(bus(:, LAM_Q));
        [maxv, maxi] = max(bus(:, LAM_Q));
        fprintf(fd, '\nLambda Q        %8.2f $/MWh @ bus %-4d   %8.2f $/MWh @ bus %-4d', minv, bus(mini, BUS_I), maxv, bus(maxi, BUS_I));
    end
    fprintf(fd, '\n');
end

if OUT_AREA_SUM
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n|     Area Summary                                                             |');
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\nArea  # of      # of Gens        # of Loads         # of    # of   # of   # of');
    fprintf(fd, '\n Num  Buses   Total  Online   Total  Fixed  Disp    Shunt   Brchs  Xfmrs   Ties');
    fprintf(fd, '\n----  -----   -----  ------   -----  -----  -----   -----   -----  -----  -----');
    for i=1:length(s_areas)
        a = s_areas(i);
        ib = find(bus(:, BUS_AREA) == a);
        ig = find(bus(e2i(gen(:, GEN_BUS)), BUS_AREA) == a & ~isload(gen));
        igon = find(bus(e2i(gen(:, GEN_BUS)), BUS_AREA) == a & gen(:, GEN_STATUS) > 0 & ~isload(gen));
        ildon = find(bus(e2i(gen(:, GEN_BUS)), BUS_AREA) == a & gen(:, GEN_STATUS) > 0 & isload(gen));
        inzld = find(bus(:, BUS_AREA) == a & (bus(:, PD) | bus(:, QD)));
        inzsh = find(bus(:, BUS_AREA) == a & (bus(:, GS) | bus(:, BS)));
        ibrch = find(bus(e2i(branch(:, F_BUS)), BUS_AREA) == a & bus(e2i(branch(:, T_BUS)), BUS_AREA) == a);
        in_tie = find(bus(e2i(branch(:, F_BUS)), BUS_AREA) == a & bus(e2i(branch(:, T_BUS)), BUS_AREA) ~= a);
        out_tie = find(bus(e2i(branch(:, F_BUS)), BUS_AREA) ~= a & bus(e2i(branch(:, T_BUS)), BUS_AREA) == a);
        if length(xfmr)
            nxfmr = length(find(bus(e2i(branch(xfmr, F_BUS)), BUS_AREA) == a & bus(e2i(branch(xfmr, T_BUS)), BUS_AREA) == a));
        else
            nxfmr = 0;
        end
        fprintf(fd, '\n%3d  %6d   %5d  %5d   %5d  %5d  %5d   %5d   %5d  %5d  %5d', ...
            a, length(ib), length(ig), length(igon), ...
            length(inzld)+length(ildon), length(inzld), length(ildon), ...
            length(inzsh), length(ibrch), nxfmr, length(in_tie)+length(out_tie));
    end
    fprintf(fd, '\n----  -----   -----  ------   -----  -----  -----   -----   -----  -----  -----');
    fprintf(fd, '\nTot: %6d   %5d  %5d   %5d  %5d  %5d   %5d   %5d  %5d  %5d', ...
        nb, length(allg), length(ong), length(nzld)+length(onld), ...
        length(nzld), length(onld), length(nzsh), nl, length(xfmr), length(ties));
    fprintf(fd, '\n');
    fprintf(fd, '\nArea      Total Gen Capacity           On-line Gen Capacity         Generation');
    fprintf(fd, '\n Num     MW           MVAr            MW           MVAr             MW    MVAr');
    fprintf(fd, '\n----   ------  ------------------   ------  ------------------    ------  ------');
    for i=1:length(s_areas)
        a = s_areas(i);
        ig = find(bus(e2i(gen(:, GEN_BUS)), BUS_AREA) == a & ~isload(gen));
        igon = find(bus(e2i(gen(:, GEN_BUS)), BUS_AREA) == a & gen(:, GEN_STATUS) > 0 & ~isload(gen));
        fprintf(fd, '\n%3d   %7.1f  %7.1f to %-7.1f  %7.1f  %7.1f to %-7.1f   %7.1f %7.1f', ...
            a, sum(gen(ig, PMAX)), sum(gen(ig, QMIN)), sum(gen(ig, QMAX)), ...
            sum(gen(igon, PMAX)), sum(gen(igon, QMIN)), sum(gen(igon, QMAX)), ...
            sum(gen(igon, PG)), sum(gen(igon, QG)) );
    end
    fprintf(fd, '\n----   ------  ------------------   ------  ------------------    ------  ------');
    fprintf(fd, '\nTot:  %7.1f  %7.1f to %-7.1f  %7.1f  %7.1f to %-7.1f   %7.1f %7.1f', ...
            sum(gen(allg, PMAX)), sum(gen(allg, QMIN)), sum(gen(allg, QMAX)), ...
            sum(gen(ong, PMAX)), sum(gen(ong, QMIN)), sum(gen(ong, QMAX)), ...
            sum(gen(ong, PG)), sum(gen(ong, QG)) );
    fprintf(fd, '\n');
    fprintf(fd, '\nArea    Disp Load Cap       Disp Load         Fixed Load        Total Load');
    fprintf(fd, '\n Num      MW     MVAr       MW     MVAr       MW     MVAr       MW     MVAr');
    fprintf(fd, '\n----    ------  ------    ------  ------    ------  ------    ------  ------');
    Qlim = (gen(:, QMIN) == 0) .* gen(:, QMAX) + (gen(:, QMAX) == 0) .* gen(:, QMIN);
    for i=1:length(s_areas)
        a = s_areas(i);
        ildon = find(bus(e2i(gen(:, GEN_BUS)), BUS_AREA) == a & gen(:, GEN_STATUS) > 0 & isload(gen));
        inzld = find(bus(:, BUS_AREA) == a & (bus(:, PD) | bus(:, QD)));
        fprintf(fd, '\n%3d    %7.1f %7.1f   %7.1f %7.1f   %7.1f %7.1f   %7.1f %7.1f', ...
            a, -sum(gen(ildon, PMIN)), ...
            -sum(Qlim(ildon)), ...
            -sum(gen(ildon, PG)), -sum(gen(ildon, QG)), ...
            sum(bus(inzld, PD)), sum(bus(inzld, QD)), ...
            -sum(gen(ildon, PG)) + sum(bus(inzld, PD)), ...
            -sum(gen(ildon, QG)) + sum(bus(inzld, QD)) );
    end
    fprintf(fd, '\n----    ------  ------    ------  ------    ------  ------    ------  ------');
    fprintf(fd, '\nTot:   %7.1f %7.1f   %7.1f %7.1f   %7.1f %7.1f   %7.1f %7.1f', ...
            -sum(gen(onld, PMIN)), ...
            -sum(Qlim(onld)), ...
            -sum(gen(onld, PG)), -sum(gen(onld, QG)), ...
            sum(bus(nzld, PD)), sum(bus(nzld, QD)), ...
            -sum(gen(onld, PG)) + sum(bus(nzld, PD)), ...
            -sum(gen(onld, QG)) + sum(bus(nzld, QD)) );
    fprintf(fd, '\n');
    fprintf(fd, '\nArea      Shunt Inj        Branch      Series Losses      Net Export');
    fprintf(fd, '\n Num      MW     MVAr     Charging      MW     MVAr       MW     MVAr');
    fprintf(fd, '\n----    ------  ------    --------    ------  ------    ------  ------');
    for i=1:length(s_areas)
        a = s_areas(i);
        inzsh = find(bus(:, BUS_AREA) == a & (bus(:, GS) | bus(:, BS)));
        ibrch = find(bus(e2i(branch(:, F_BUS)), BUS_AREA) == a & bus(e2i(branch(:, T_BUS)), BUS_AREA) == a & branch(:, BR_STATUS));
        in_tie = find(bus(e2i(branch(:, F_BUS)), BUS_AREA) ~= a & bus(e2i(branch(:, T_BUS)), BUS_AREA) == a & branch(:, BR_STATUS));
        out_tie = find(bus(e2i(branch(:, F_BUS)), BUS_AREA) == a & bus(e2i(branch(:, T_BUS)), BUS_AREA) ~= a & branch(:, BR_STATUS));
        fprintf(fd, '\n%3d    %7.1f %7.1f    %7.1f    %7.2f %7.2f   %7.1f %7.1f', ...
            a, -sum(bus(inzsh, VM) .^ 2 .* bus(inzsh, GS)), ...
            sum(bus(inzsh, VM) .^ 2 .* bus(inzsh, BS)), ...
            sum(fchg(ibrch)) + sum(tchg(ibrch)) + sum(fchg(out_tie)) + sum(tchg(in_tie)), ...
            sum(real(loss(ibrch))) + sum(real(loss([in_tie; out_tie]))) / 2, ...
            sum(imag(loss(ibrch))) + sum(imag(loss([in_tie; out_tie]))) / 2, ...
            sum(branch(in_tie, PT))+sum(branch(out_tie, PF)) - sum(real(loss([in_tie; out_tie]))) / 2, ...
            sum(branch(in_tie, QT))+sum(branch(out_tie, QF)) - sum(imag(loss([in_tie; out_tie]))) / 2  );
    end
    fprintf(fd, '\n----    ------  ------    --------    ------  ------    ------  ------');
    fprintf(fd, '\nTot:   %7.1f %7.1f    %7.1f    %7.2f %7.2f       -       -', ...
        -sum(bus(nzsh, VM) .^ 2 .* bus(nzsh, GS)), ...
        sum(bus(nzsh, VM) .^ 2 .* bus(nzsh, BS)), ...
        sum(fchg) + sum(tchg), sum(real(loss)), sum(imag(loss)) );
    fprintf(fd, '\n');
end

%% generator data
if OUT_GEN
    if isOPF
        genlamP = bus(e2i(gen(:, GEN_BUS)), LAM_P);
        genlamQ = bus(e2i(gen(:, GEN_BUS)), LAM_Q);
    end
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n|     Generator Data                                                           |');
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n Gen   Bus   Status     Pg        Qg   ');
    if isOPF, fprintf(fd, '   Lambda ($/MVA-hr)'); end
    fprintf(fd, '\n  #     #              (MW)     (MVAr) ');
    if isOPF, fprintf(fd, '     P         Q    '); end
    fprintf(fd, '\n----  -----  ------  --------  --------');
    if isOPF, fprintf(fd, '  --------  --------'); end
    for k = 1:length(ong)
        i = ong(k);
        fprintf(fd, '\n%3d %6d     %2d ', i, gen(i, GEN_BUS), gen(i, GEN_STATUS));
        if gen(i, GEN_STATUS) > 0 & (gen(i, PG) | gen(i, QG))
            fprintf(fd, '%10.2f%10.2f', gen(i, PG), gen(i, QG));
        else
            fprintf(fd, '       -         -  ');
        end
        if isOPF, fprintf(fd, '%10.2f%10.2f', genlamP(i), genlamQ(i)); end
    end
    fprintf(fd, '\n                     --------  --------');
    fprintf(fd, '\n            Total: %9.2f%10.2f', sum(gen(ong, PG)), sum(gen(ong, QG)));
    fprintf(fd, '\n');
    if length(onld) > 0
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n|     Dispatchable Load Data                                                   |');
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n Gen   Bus   Status     Pd        Qd   ');
        if isOPF, fprintf(fd, '   Lambda ($/MVA-hr)'); end
        fprintf(fd, '\n  #     #              (MW)     (MVAr) ');
        if isOPF, fprintf(fd, '     P         Q    '); end
        fprintf(fd, '\n----  -----  ------  --------  --------');
        if isOPF, fprintf(fd, '  --------  --------'); end
        for k = 1:length(onld)
            i = onld(k);
            fprintf(fd, '\n%3d %6d     %2d ', i, gen(i, GEN_BUS), gen(i, GEN_STATUS));
            if gen(i, GEN_STATUS) > 0 & (gen(i, PG) | gen(i, QG))
                fprintf(fd, '%10.2f%10.2f', -gen(i, PG), -gen(i, QG));
            else
                fprintf(fd, '       -         -  ');
            end
            if isOPF, fprintf(fd, '%10.2f%10.2f', genlamP(i), genlamQ(i)); end
        end
        fprintf(fd, '\n                     --------  --------');
        fprintf(fd, '\n            Total: %9.2f%10.2f', -sum(gen(onld, PG)), -sum(gen(onld, QG)));
        fprintf(fd, '\n');
    end
end
        
%% bus data
if OUT_BUS
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n|     Bus Data                                                                 |');
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n Bus      Voltage          Generation             Load        ');
    if isOPF, fprintf(fd, '  Lambda($/MVA-hr)'); end
    fprintf(fd, '\n  #   Mag(pu) Ang(deg)   P (MW)   Q (MVAr)   P (MW)   Q (MVAr)');
    if isOPF, fprintf(fd, '     P        Q   '); end
    fprintf(fd, '\n----- ------- --------  --------  --------  --------  --------');
    if isOPF, fprintf(fd, '  -------  -------'); end
    for i = 1:nb
        fprintf(fd, '\n%5d%7.3f%9.3f', bus(i, [BUS_I, VM, VA]));
        g  = find(gen(:, GEN_STATUS) > 0 & gen(:, GEN_BUS) == bus(i, BUS_I) & ...
                    ~isload(gen));
        ld = find(gen(:, GEN_STATUS) > 0 & gen(:, GEN_BUS) == bus(i, BUS_I) & ...
                    isload(gen));
        if ~isempty(g)
            fprintf(fd, '%10.2f%10.2f', sum(gen(g, PG)), sum(gen(g, QG)));
        else
            fprintf(fd, '       -         -  ');
        end
        if bus(i, PD) | bus(i, QD) | ~isempty(ld)
            if ~isempty(ld)
                fprintf(fd, '%10.2f*%9.2f*', bus(i, PD) - sum(gen(ld, PG)), ...
                                             bus(i, QD) - sum(gen(ld, QG)));
            else
                fprintf(fd, '%10.2f%10.2f ', bus(i, [PD, QD]));
            end
        else
            fprintf(fd, '       -         -   ');
        end
        if isOPF
            fprintf(fd, '%9.3f', bus(i, LAM_P));
            if abs(bus(i, LAM_Q)) > 1e-6
                fprintf(fd, '%8.3f', bus(i, LAM_Q));
            else
                fprintf(fd, '     -');
            end
        end
    end
    fprintf(fd, '\n                        --------  --------  --------  --------');
    fprintf(fd, '\n               Total: %9.2f %9.2f %9.2f %9.2f', ...
        sum(gen(ong, PG)), sum(gen(ong, QG)), ...
        sum(bus(nzld, PD)) - sum(gen(onld, PG)), ...
        sum(bus(nzld, QD)) - sum(gen(onld, QG)));
    fprintf(fd, '\n');
end

%% branch data
if OUT_BRANCH
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n|     Branch Data                                                              |');
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\nBrnch   From   To    From Bus Injection   To Bus Injection     Loss (I^2 * Z)  ');
    fprintf(fd, '\n  #     Bus    Bus    P (MW)   Q (MVAr)   P (MW)   Q (MVAr)   P (MW)   Q (MVAr)');
    fprintf(fd, '\n-----  -----  -----  --------  --------  --------  --------  --------  --------');
    fprintf(fd, '\n%4d%7d%7d%10.2f%10.2f%10.2f%10.2f%10.3f%10.2f', ...
            [   [1:nl]', branch(:, [F_BUS, T_BUS]), ...
                branch(:, [PF, QF]), branch(:, [PT, QT]), ...
                real(loss), imag(loss) ...
            ]');
    fprintf(fd, '\n                                                             --------  --------');
    fprintf(fd, '\n                                                    Total:%10.3f%10.2f', ...
            sum(real(loss)), sum(imag(loss)));
    fprintf(fd, '\n');
end
    
%%-----  constraint data  -----
if isOPF
    ctol = mpopt(16);   %% constraint violation tolerance
    %% voltage constraints
    if OUT_V_LIM == 2 | (OUT_V_LIM == 1 & ...
                         (any(bus(:, VM) < bus(:, VMIN) + ctol) | ...
                          any(bus(:, VM) > bus(:, VMAX) - ctol) | ...
                          any(bus(:, MU_VMIN) > 1e-6) | ...
                          any(bus(:, MU_VMAX) > 1e-6)))
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n|     Voltage Constraints                                                      |');
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\nBus #  Vmin mu    Vmin    |V|   Vmax    Vmax mu');
        fprintf(fd, '\n-----  --------   -----  -----  -----   --------');        
        for i = 1:nb
            if OUT_V_LIM == 2 | (OUT_V_LIM == 1 & ...
                         (bus(i, VM) < bus(i, VMIN) + ctol | ...
                          bus(i, VM) > bus(i, VMAX) - ctol | ...
                          bus(i, MU_VMIN) > 1e-6 | bus(i, MU_VMAX) > 1e-6))
                fprintf(fd, '\n%5d', bus(i, BUS_I));
                if bus(i, VM) < bus(i, VMIN) + ctol | bus(i, MU_VMIN) > 1e-6
                    fprintf(fd, '%10.3f', bus(i, MU_VMIN));                    
                else
                    fprintf(fd, '      -   ');
                end
                fprintf(fd, '%8.3f%7.3f%7.3f', bus(i, [VMIN, VM, VMAX]));
                if bus(i, VM) > bus(i, VMAX) - ctol | bus(i, MU_VMAX) > 1e-6
                    fprintf(fd, '%10.3f', bus(i, MU_VMAX));
                else
                    fprintf(fd, '      -    ');
                end
            end
        end
        fprintf(fd, '\n');
    end
        
    %% generator P constraints
    if OUT_PG_LIM == 2 | OUT_QG_LIM == 2 | ...
            (OUT_PG_LIM == 1 & (any(gen(ong, PG) < gen(ong, PMIN) + ctol) | ...
                                any(gen(ong, PG) > gen(ong, PMAX) - ctol) | ...
                                any(gen(ong, MU_PMIN) > 1e-6) | ...
                                any(gen(ong, MU_PMAX) > 1e-6))) | ...
            (OUT_QG_LIM == 1 & (any(gen(ong, QG) < gen(ong, QMIN) + ctol) | ...
                                any(gen(ong, QG) > gen(ong, QMAX) - ctol) | ...
                                any(gen(ong, MU_QMIN) > 1e-6) | ...
                                any(gen(ong, MU_QMAX) > 1e-6)))
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n|     Generation Constraints                                                   |');
        fprintf(fd, '\n================================================================================');
    end
    if OUT_PG_LIM == 2 | (OUT_PG_LIM == 1 & ...
                             (any(gen(ong, PG) < gen(ong, PMIN) + ctol) | ...
                              any(gen(ong, PG) > gen(ong, PMAX) - ctol) | ...
                              any(gen(ong, MU_PMIN) > 1e-6) | ...
                              any(gen(ong, MU_PMAX) > 1e-6)))
        fprintf(fd, '\n Gen   Bus                Active Power Limits');
        fprintf(fd, '\n  #     #    Pmin mu    Pmin       Pg       Pmax    Pmax mu');
        fprintf(fd, '\n----  -----  -------  --------  --------  --------  -------');
        for k = 1:length(ong)
            i = ong(k);
            if OUT_PG_LIM == 2 | (OUT_PG_LIM == 1 & ...
                        (gen(i, PG) < gen(i, PMIN) + ctol | ...
                         gen(i, PG) > gen(i, PMAX) - ctol | ...
                         gen(i, MU_PMIN) > 1e-6 | gen(i, MU_PMAX) > 1e-6))
                fprintf(fd, '\n%4d%6d ', i, gen(i, GEN_BUS));
                if gen(i, PG) < gen(i, PMIN) + ctol | gen(i, MU_PMIN) > 1e-6
                    fprintf(fd, '%8.3f', gen(i, MU_PMIN));
                else
                    fprintf(fd, '     -  ');
                end
                if gen(i, PG)
                    fprintf(fd, '%10.2f%10.2f%10.2f', gen(i, [PMIN, PG, PMAX]));
                else
                    fprintf(fd, '%10.2f       -  %10.2f', gen(i, [PMIN, PMAX]));
                end
                if gen(i, PG) > gen(i, PMAX) - ctol | gen(i, MU_PMAX) > 1e-6
                    fprintf(fd, '%9.3f', gen(i, MU_PMAX));
                else
                    fprintf(fd, '      -  ');
                end
            end
        end
        fprintf(fd, '\n');
    end
        
    %% generator Q constraints
    if OUT_QG_LIM == 2 | (OUT_QG_LIM == 1 & ...
                             (any(gen(ong, QG) < gen(ong, QMIN) + ctol) | ...
                              any(gen(ong, QG) > gen(ong, QMAX) - ctol) | ...
                              any(gen(ong, MU_QMIN) > 1e-6) | ...
                              any(gen(ong, MU_QMAX) > 1e-6)))
        fprintf(fd, '\nGen  Bus              Reactive Power Limits');
        fprintf(fd, '\n #    #   Qmin mu    Qmin       Qg       Qmax    Qmax mu');
        fprintf(fd, '\n---  ---  -------  --------  --------  --------  -------');
        for k = 1:length(ong)
            i = ong(k);
            if OUT_QG_LIM == 2 | (OUT_QG_LIM == 1 & ...
                        (gen(i, QG) < gen(i, QMIN) + ctol | ...
                         gen(i, QG) > gen(i, QMAX) - ctol | ...
                         gen(i, MU_QMIN) > 1e-6 | gen(i, MU_QMAX) > 1e-6))
                fprintf(fd, '\n%3d%5d', i, gen(i, GEN_BUS));
                if gen(i, QG) < gen(i, QMIN) + ctol | gen(i, MU_QMIN) > 1e-6
                    fprintf(fd, '%8.3f', gen(i, MU_QMIN));
                else
                    fprintf(fd, '     -  ');
                end
                if gen(i, QG)
                    fprintf(fd, '%10.2f%10.2f%10.2f', gen(i, [QMIN, QG, QMAX]));
                else
                    fprintf(fd, '%10.2f       -  %10.2f', gen(i, [QMIN, QMAX]));
                end
                if gen(i, QG) > gen(i, QMAX) - ctol | gen(i, MU_QMAX) > 1e-6
                    fprintf(fd, '%9.3f', gen(i, MU_QMAX));
                else
                    fprintf(fd, '      -  ');
                end
            end
        end
        fprintf(fd, '\n');
    end
        
    %% dispatchable load P constraints
    if OUT_PG_LIM == 2 | OUT_QG_LIM == 2 | ...
            (OUT_PG_LIM == 1 & (any(gen(onld, PG) < gen(onld, PMIN) + ctol) | ...
                                any(gen(onld, PG) > gen(onld, PMAX) - ctol) | ...
                                any(gen(onld, MU_PMIN) > 1e-6) | ...
                                any(gen(onld, MU_PMAX) > 1e-6))) | ...
            (OUT_QG_LIM == 1 & (any(gen(onld, QG) < gen(onld, QMIN) + ctol) | ...
                                any(gen(onld, QG) > gen(onld, QMAX) - ctol) | ...
                                any(gen(onld, MU_QMIN) > 1e-6) | ...
                                any(gen(onld, MU_QMAX) > 1e-6)))
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n|     Dispatchable Load Constraints                                            |');
        fprintf(fd, '\n================================================================================');
    end
    if OUT_PG_LIM == 2 | (OUT_PG_LIM == 1 & ...
                             (any(gen(onld, PG) < gen(onld, PMIN) + ctol) | ...
                              any(gen(onld, PG) > gen(onld, PMAX) - ctol) | ...
                              any(gen(onld, MU_PMIN) > 1e-6) | ...
                              any(gen(onld, MU_PMAX) > 1e-6)))
        fprintf(fd, '\nGen  Bus               Active Power Limits');
        fprintf(fd, '\n #    #   Pmin mu    Pmin       Pg       Pmax    Pmax mu');
        fprintf(fd, '\n---  ---  -------  --------  --------  --------  -------');
        for k = 1:length(onld)
            i = onld(k);
            if OUT_PG_LIM == 2 | (OUT_PG_LIM == 1 & ...
                        (gen(i, PG) < gen(i, PMIN) + ctol | ...
                         gen(i, PG) > gen(i, PMAX) - ctol | ...
                         gen(i, MU_PMIN) > 1e-6 | gen(i, MU_PMAX) > 1e-6))
                fprintf(fd, '\n%3d%5d', i, gen(i, GEN_BUS));
                if gen(i, PG) < gen(i, PMIN) + ctol | gen(i, MU_PMIN) > 1e-6
                    fprintf(fd, '%8.3f', gen(i, MU_PMIN));
                else
                    fprintf(fd, '     -  ');
                end
                if gen(i, PG)
                    fprintf(fd, '%10.2f%10.2f%10.2f', gen(i, [PMIN, PG, PMAX]));
                else
                    fprintf(fd, '%10.2f       -  %10.2f', gen(i, [PMIN, PMAX]));
                end
                if gen(i, PG) > gen(i, PMAX) - ctol | gen(i, MU_PMAX) > 1e-6
                    fprintf(fd, '%9.3f', gen(i, MU_PMAX));
                else
                    fprintf(fd, '      -  ');
                end
            end
        end
        fprintf(fd, '\n');
    end
        
    %% dispatchable load Q constraints
    if OUT_QG_LIM == 2 | (OUT_QG_LIM == 1 & ...
                             (any(gen(onld, QG) < gen(onld, QMIN) + ctol) | ...
                              any(gen(onld, QG) > gen(onld, QMAX) - ctol) | ...
                              any(gen(onld, MU_QMIN) > 1e-6) | ...
                              any(gen(onld, MU_QMAX) > 1e-6)))
        fprintf(fd, '\nGen  Bus              Reactive Power Limits');
        fprintf(fd, '\n #    #   Qmin mu    Qmin       Qg       Qmax    Qmax mu');
        fprintf(fd, '\n---  ---  -------  --------  --------  --------  -------');
        for k = 1:length(onld)
            i = onld(k);
            if OUT_QG_LIM == 2 | (OUT_QG_LIM == 1 & ...
                        (gen(i, QG) < gen(i, QMIN) + ctol | ...
                         gen(i, QG) > gen(i, QMAX) - ctol | ...
                         gen(i, MU_QMIN) > 1e-6 | gen(i, MU_QMAX) > 1e-6))
                fprintf(fd, '\n%3d%5d', i, gen(i, GEN_BUS));
                if gen(i, QG) < gen(i, QMIN) + ctol | gen(i, MU_QMIN) > 1e-6
                    fprintf(fd, '%8.3f', gen(i, MU_QMIN));
                else
                    fprintf(fd, '     -  ');
                end
                if gen(i, QG)
                    fprintf(fd, '%10.2f%10.2f%10.2f', gen(i, [QMIN, QG, QMAX]));
                else
                    fprintf(fd, '%10.2f       -  %10.2f', gen(i, [QMIN, QMAX]));
                end
                if gen(i, QG) > gen(i, QMAX) - ctol | gen(i, MU_QMAX) > 1e-6
                    fprintf(fd, '%9.3f', gen(i, MU_QMAX));
                else
                    fprintf(fd, '      -  ');
                end
            end
        end
        fprintf(fd, '\n');
    end
        
    %% line flow constraints
    if mpopt(24) == 1   %% P limit
        Sf = branch(:, PF);
        St = branch(:, PT);
        str = '\n  #     Bus    Pf  mu     Pf      |Pmax|      Pt      Pt  mu   Bus';
    else                %% |S| limit
        Sf = abs(branch(:, PF) + j * branch(:, QF));
        St = abs(branch(:, PT) + j * branch(:, QT));
        str = '\n  #     Bus   |Sf| mu    |Sf|     |Smax|     |St|    |St| mu   Bus';
    end
    if OUT_LINE_LIM == 2 | (OUT_LINE_LIM == 1 & ...
                        (any(abs(Sf) > branch(:, RATE_A) - ctol) | ...
                         any(abs(St) > branch(:, RATE_A) - ctol) | ...
                         any(branch(:, MU_SF) > 1e-6) | ...
                         any(branch(:, MU_ST) > 1e-6)))
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n|     Branch Flow Constraints                                                  |');
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\nBrnch   From     "From" End        Limit       "To" End        To');
        fprintf(fd, str);
        fprintf(fd, '\n-----  -----  -------  --------  --------  --------  -------  -----');
        for i = 1:nl
            if OUT_LINE_LIM == 2 | (OUT_LINE_LIM == 1 & ...
                         (Sf(i) > branch(i, RATE_A) - ctol | ...
                          St(i) > branch(i, RATE_A) - ctol | ...
                          branch(i, MU_SF) > 1e-6 | branch(i, MU_ST) > 1e-6))
                fprintf(fd, '\n%4d%7d', i, branch(i, F_BUS));
                if Sf(i) > branch(i, RATE_A) - ctol | branch(i, MU_SF) > 1e-6
                    fprintf(fd, '%10.3f', branch(i, MU_SF));
                else
                    fprintf(fd, '      -   ');
                end
                fprintf(fd, '%9.2f%10.2f%10.2f', ...
                    [Sf(i), branch(i, RATE_A), St(i)]);
                if St(i) > branch(i, RATE_A) - ctol | branch(i, MU_ST) > 1e-6
                    fprintf(fd, '%10.3f', branch(i, MU_ST));
                else
                    fprintf(fd, '      -   ');
                end
                fprintf(fd, '%6d', branch(i, T_BUS));
            end
        end
        fprintf(fd, '\n');
    end
end

%% print raw data for Perl database interface
if OUT_RAW
    fprintf(fd, '----------  raw PB::Soln data below  ----------\n');
    fprintf(fd, 'bus\n');
    if isOPF
        fprintf(fd, '%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\n', ...
                    bus(:, [BUS_I, BUS_TYPE, VM, VA, LAM_P, LAM_Q, MU_VMAX, MU_VMIN])');
    
        fprintf(fd, 'branch\n');
        fprintf(fd, '%d\t%g\t%g\t%g\t%g\t%g\t%g\n', ...
                    [[1:nl]' branch(:, [PF, QF, PT, QT, MU_SF, MU_ST])]');
    
        fprintf(fd, 'gen\n');
        fprintf(fd, '%d\t%g\t%g\t%g\t%d\t%g\t%g\t%g\t%g\n', ...
                    [[1:ng]' gen(:, [PG, QG, VG, GEN_STATUS, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN])]');
    else
        fprintf(fd, '%d\t%d\t%f\t%f\t%d\t%d\t%d\t%d\n', ...
                    [bus(:, [BUS_I, BUS_TYPE, VM, VA]) zeros(nb, 4)]');
    
        fprintf(fd, 'branch\n');
        fprintf(fd, '%d\t%f\t%f\t%f\t%f\t%d\t%d\n', ...
                    [[1:nl]' branch(:, [PF, QF, PT, QT]) zeros(nl, 2)]');
    
        fprintf(fd, 'gen\n');
        fprintf(fd, '%d\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\n', ...
                    [[1:ng]' gen(:, [PG, QG, VG, GEN_STATUS]) zeros(ng, 4)]');
    end
    fprintf(fd, 'end\n');
    fprintf(fd, '----------  raw PB::Soln data above  ----------\n');
end

return;
