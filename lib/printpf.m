function printpf(baseMVA, bus, gen, branch, f, success, et, fd, mpopt)
%PRINTPF   Prints power flow results.
%   PRINTPF(RESULTS, FD, MPOPT)
%   PRINTPF(BASEMVA, BUS, GEN, BRANCH, F, SUCCESS, ET, FD, MPOPT)
%
%   Prints power flow and optimal power flow results to FD (a file
%   descriptor which defaults to STDOUT), with the details of what
%   gets printed controlled by the optional MPOPT argument, which is a
%   MATPOWER options struct (see MPOPTION for details).
%
%   The data can either be supplied in a single RESULTS struct, or
%   in the individual arguments: BASEMVA, BUS, GEN, BRANCH, F, SUCCESS
%   and ET, where F is the OPF objective function value, SUCCESS is
%   true if the solution converged and false otherwise, and ET is the
%   elapsed time for the computation in seconds. If F is given, it is
%   assumed that the output is from an OPF run, otherwise it is assumed
%   to be a simple power flow run.
%
%   Examples:
%       mpopt = mpoptions('out.gen', 1, 'out.bus', 0, 'out.branch', 0);
%       [fd, msg] = fopen(fname, 'at');
%       results = runopf(mpc);
%       printpf(results);
%       printpf(results, fd);
%       printpf(results, fd, mpopt);
%       printpf(baseMVA, bus, gen, branch, f, success, et);
%       printpf(baseMVA, bus, gen, branch, f, success, et, fd);
%       printpf(baseMVA, bus, gen, branch, f, success, et, fd, mpopt);
%       fclose(fd);

%   MATPOWER
%   Copyright (c) 1996-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialization -----
%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% default arguments
if isstruct(baseMVA)
    have_results_struct = 1;
    results = baseMVA;
    if nargin < 3 || isempty(gen)
        mpopt = mpoption;   %% use default options
    else
        mpopt = gen;
    end
    if mpopt.out.all == 0
        return;         %% nothin' to see here, bail out now
    end
    if nargin < 2 || isempty(bus)
        fd = 1;         %% print to stdio by default
    else
        fd = bus;
    end
    [baseMVA, bus, gen, branch, success, et] = ...
        deal(results.baseMVA, results.bus, results.gen, results.branch, ...
            results.success, results.et);
    if isfield(results, 'f') && ~isempty(results.f)
        f = results.f;
    else
        f = [];
    end
else
    have_results_struct = 0;
    if nargin < 9
        mpopt = mpoption;   %% use default options
        if nargin < 8
            fd = 1;         %% print to stdio by default
        end
    end
    if mpopt.out.all == 0
        return;         %% nothin' to see here, bail out now
    end
end
isOPF = ~isempty(f);    %% FALSE -> only simple PF data, TRUE -> OPF data

%% options
isDC            = strcmp(upper(mpopt.model), 'DC');

SUPPRESS        = mpopt.out.suppress_detail;
if SUPPRESS == -1
    if size(bus, 1) > 500
        SUPPRESS = 1;
    else
        SUPPRESS = 0;
    end
end
OUT_ALL         = mpopt.out.all;
OUT_FORCE       = mpopt.out.force;
OUT_ANY         = OUT_ALL == 1;     %% set to true if any pretty output is to be generated
OUT_SYS_SUM     = OUT_ALL == 1 || (OUT_ALL == -1 && mpopt.out.sys_sum);
OUT_AREA_SUM    = OUT_ALL == 1 || (OUT_ALL == -1 && ~SUPPRESS && mpopt.out.area_sum);
OUT_BUS         = OUT_ALL == 1 || (OUT_ALL == -1 && ~SUPPRESS && mpopt.out.bus);
OUT_BRANCH      = OUT_ALL == 1 || (OUT_ALL == -1 && ~SUPPRESS && mpopt.out.branch);
OUT_GEN         = OUT_ALL == 1 || (OUT_ALL == -1 && ~SUPPRESS && mpopt.out.gen);
OUT_ANY         = OUT_ANY || (OUT_ALL == -1 && ...
                    (OUT_SYS_SUM || OUT_AREA_SUM || OUT_BUS || ...
                    OUT_BRANCH || OUT_GEN));
if OUT_ALL == -1
    OUT_ALL_LIM = ~SUPPRESS * mpopt.out.lim.all;
elseif OUT_ALL == 1
    OUT_ALL_LIM = 2;
else
    OUT_ALL_LIM = 0;
end
OUT_ANY         = OUT_ANY || OUT_ALL_LIM >= 1;
if OUT_ALL_LIM == -1
    OUT_V_LIM       = ~SUPPRESS * mpopt.out.lim.v;
    OUT_LINE_LIM    = ~SUPPRESS * mpopt.out.lim.line;
    OUT_PG_LIM      = ~SUPPRESS * mpopt.out.lim.pg;
    OUT_QG_LIM      = ~SUPPRESS * mpopt.out.lim.qg;
else
    OUT_V_LIM       = OUT_ALL_LIM;
    OUT_LINE_LIM    = OUT_ALL_LIM;
    OUT_PG_LIM      = OUT_ALL_LIM;
    OUT_QG_LIM      = OUT_ALL_LIM;
end
OUT_ANY         = OUT_ANY || (OUT_ALL_LIM == -1 && (OUT_V_LIM || OUT_LINE_LIM || OUT_PG_LIM || OUT_QG_LIM));

%%----- print the stuff -----
if OUT_ANY
    ptol = 1e-4;        %% tolerance for displaying shadow prices
    if isOPF && ~isDC && strcmp(upper(mpopt.opf.ac.solver), 'SDPOPF')
        isSDP = 1;
        ptol = 0.1;     %% tolerance for displaying shadow prices
        if have_results_struct && isfield(results, 'mineigratio') && ~isempty(results.mineigratio)
            mineigratio = results.mineigratio;
        else
            mineigratio = [];
        end
        if have_results_struct && isfield(results, 'zero_eval') && ~isempty(results.zero_eval)
            zero_eval = results.zero_eval;
        else
            zero_eval = [];
        end
    else
        isSDP = 0;
    end

    %% create map of external bus numbers to bus indices
    i2e = bus(:, BUS_I);
    e2i = sparse(max(i2e), 1);
    e2i(i2e) = (1:size(bus, 1))';

    %% sizes of things
    nb = size(bus, 1);      %% number of buses
    nl = size(branch, 1);   %% number of branches
    ng = size(gen, 1);      %% number of generators

    %% zero out some data to make printout consistent for DC case
    if isDC
        bus(:, [QD, BS])            = zeros(nb, 2);
        gen(:, [QG, QMAX, QMIN])    = zeros(ng, 3);
        branch(:, [BR_R, BR_B])     = zeros(nl, 2);
    end

    %% parameters
    ties = find(bus(e2i(branch(:, F_BUS)), BUS_AREA) ~= bus(e2i(branch(:, T_BUS)), BUS_AREA));
                            %% area inter-ties
    xfmr = find(branch(:, TAP));                    %% indices of transformers
    nzld = find((bus(:, PD) | bus(:, QD)) & bus(:, BUS_TYPE) ~= NONE);
    sorted_areas = sort(bus(:, BUS_AREA));
    s_areas = sorted_areas([1; find(diff(sorted_areas))+1]);    %% area numbers
    nzsh = find((bus(:, GS) | bus(:, BS)) & bus(:, BUS_TYPE) ~= NONE);
    allg = find( ~isload(gen) );
    alld = find(  isload(gen) );
    ong  = find( gen(:, GEN_STATUS) > 0 & ~isload(gen) );
    onld = find( gen(:, GEN_STATUS) > 0 &  isload(gen) );
    V = bus(:, VM) .* exp(sqrt(-1) * pi/180 * bus(:, VA));
    opt = struct('type', 'FIXED');
    [Pdf, Qdf] = total_load(bus, gen, 'bus', struct('type', 'FIXED'), mpopt);
    [Pdd, Qdd] = total_load(bus, gen, 'bus', struct('type', 'DISPATCHABLE'), mpopt);
    if isDC
        loss = zeros(nl, 1);
        fchg = loss;
        tchg = loss;
    else
        [loss, fchg, tchg] = get_losses(baseMVA, bus, branch);
    end

    %% convergence & elapsed time
    if success
        if isSDP
            fprintf(fd, '\nSolution satisfies rank and consistency conditions, %.2f seconds.\nmineigratio = %0.5g, zero_eval = %0.5g', et, mineigratio, zero_eval);
        else
            fprintf(fd, '\nConverged in %.2f seconds', et);
        end
    else
        if isSDP
            fprintf(fd, '\n>>>>>  Solution does NOT satisfy rank and/or consistency conditions (%.2f seconds).  <<<<<\nmineigratio = %0.5g, zero_eval = %0.5g\n', et, mineigratio, zero_eval);
        else
            fprintf(fd, '\n>>>>>  Did NOT converge (%.2f seconds)  <<<<<\n', et);
        end
    end

    %% objective function value
    if isOPF && (success || OUT_FORCE)
        fprintf(fd, '\nObjective Function Value = %.2f $/hr', f);
    end
end
if OUT_SYS_SUM && (success || OUT_FORCE)
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n|     System Summary                                                           |');
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n\nHow many?                How much?              P (MW)            Q (MVAr)');
    fprintf(fd, '\n---------------------    -------------------  -------------  -----------------');
    fprintf(fd, '\nBuses         %6d     Total Gen Capacity   %7.1f       %7.1f to %.1f', nb, sum(gen(allg, PMAX)), sum(gen(allg, QMIN)), sum(gen(allg, QMAX)));
    fprintf(fd, '\nGenerators     %5d     On-line Capacity     %7.1f       %7.1f to %.1f', length(allg), sum(gen(ong, PMAX)), sum(gen(ong, QMIN)), sum(gen(ong, QMAX)));
    fprintf(fd, '\nCommitted Gens %5d     Generation (actual)  %7.1f           %7.1f', length(ong), sum(gen(ong, PG)), sum(gen(ong, QG)));
    fprintf(fd, '\nLoads          %5d     Load                 %7.1f           %7.1f', length(nzld)+length(onld), sum(Pdf(nzld))-sum(gen(onld, PG)), sum(Qdf(nzld))-sum(gen(onld, QG)));
    fprintf(fd, '\n  Fixed        %5d       Fixed              %7.1f           %7.1f', length(nzld), sum(Pdf(nzld)), sum(Qdf(nzld)));
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
    if ~isDC
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

if OUT_AREA_SUM && (success || OUT_FORCE)
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
        inzld = find(bus(:, BUS_AREA) == a & (Pdf | Qdf) & bus(:, BUS_TYPE) ~= NONE);
        inzsh = find(bus(:, BUS_AREA) == a & (bus(:, GS) | bus(:, BS)) & bus(:, BUS_TYPE) ~= NONE);
        ibrch = find(bus(e2i(branch(:, F_BUS)), BUS_AREA) == a & bus(e2i(branch(:, T_BUS)), BUS_AREA) == a);
        in_tie = find(bus(e2i(branch(:, F_BUS)), BUS_AREA) == a & bus(e2i(branch(:, T_BUS)), BUS_AREA) ~= a);
        out_tie = find(bus(e2i(branch(:, F_BUS)), BUS_AREA) ~= a & bus(e2i(branch(:, T_BUS)), BUS_AREA) == a);
        if isempty(xfmr)
            nxfmr = 0;
        else
            nxfmr = length(find(bus(e2i(branch(xfmr, F_BUS)), BUS_AREA) == a & bus(e2i(branch(xfmr, T_BUS)), BUS_AREA) == a));
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
        inzld = find(bus(:, BUS_AREA) == a & (Pdf | Qdf));
        fprintf(fd, '\n%3d    %7.1f %7.1f   %7.1f %7.1f   %7.1f %7.1f   %7.1f %7.1f', ...
            a, -sum(gen(ildon, PMIN)), ...
            -sum(Qlim(ildon)), ...
            -sum(gen(ildon, PG)), -sum(gen(ildon, QG)), ...
            sum(Pdf(inzld)), sum(Qdf(inzld)), ...
            -sum(gen(ildon, PG)) + sum(Pdf(inzld)), ...
            -sum(gen(ildon, QG)) + sum(Qdf(inzld)) );
    end
    fprintf(fd, '\n----    ------  ------    ------  ------    ------  ------    ------  ------');
    fprintf(fd, '\nTot:   %7.1f %7.1f   %7.1f %7.1f   %7.1f %7.1f   %7.1f %7.1f', ...
            -sum(gen(onld, PMIN)), ...
            -sum(Qlim(onld)), ...
            -sum(gen(onld, PG)), -sum(gen(onld, QG)), ...
            sum(Pdf(nzld)), sum(Qdf(nzld)), ...
            -sum(gen(onld, PG)) + sum(Pdf(nzld)), ...
            -sum(gen(onld, QG)) + sum(Qdf(nzld)) );
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
if OUT_GEN && (success || OUT_FORCE)
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
    for k = 1:length(allg)
        i = allg(k);
        fprintf(fd, '\n%3d %6d     %2d ', i, gen(i, GEN_BUS), gen(i, GEN_STATUS));
        if gen(i, GEN_STATUS) > 0 && (gen(i, PG) || gen(i, QG))
            fprintf(fd, '%10.2f%10.2f', gen(i, PG), gen(i, QG));
        else
            fprintf(fd, '       -         -  ');
        end
        if isOPF, fprintf(fd, '%10.2f%10.2f', genlamP(i), genlamQ(i)); end
    end
    fprintf(fd, '\n                     --------  --------');
    fprintf(fd, '\n            Total: %9.2f%10.2f', sum(gen(ong, PG)), sum(gen(ong, QG)));
    fprintf(fd, '\n');
    if ~isempty(alld)
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n|     Dispatchable Load Data                                                   |');
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n Gen   Bus   Status     Pd        Qd   ');
        if isOPF, fprintf(fd, '   Lambda ($/MVA-hr)'); end
        fprintf(fd, '\n  #     #              (MW)     (MVAr) ');
        if isOPF, fprintf(fd, '     P         Q    '); end
        fprintf(fd, '\n----  -----  ------  --------  --------');
        if isOPF, fprintf(fd, '  --------  --------'); end
        for k = 1:length(alld)
            i = alld(k);
            fprintf(fd, '\n%3d %6d     %2d ', i, gen(i, GEN_BUS), gen(i, GEN_STATUS));
            if gen(i, GEN_STATUS) > 0 && (gen(i, PG) || gen(i, QG))
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
if OUT_BUS && (success || OUT_FORCE)
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
        if bus(i, BUS_TYPE) == REF
            fprintf(fd, '*');
        elseif bus(i, BUS_TYPE) == NONE
            fprintf(fd, 'x');
        else
            fprintf(fd, ' ');
        end
        g  = find(gen(:, GEN_STATUS) > 0 & gen(:, GEN_BUS) == bus(i, BUS_I) & ...
                    ~isload(gen));
        ld = find(gen(:, GEN_STATUS) > 0 & gen(:, GEN_BUS) == bus(i, BUS_I) & ...
                    isload(gen));
        if ~isempty(g)
            fprintf(fd, '%9.2f%10.2f', sum(gen(g, PG)), sum(gen(g, QG)));
        else
            fprintf(fd, '      -         -  ');
        end
        if Pdf(i) || Qdf(i) || ~isempty(ld)
            if ~isempty(ld)
                fprintf(fd, '%10.2f*%9.2f*', Pdf(i) - sum(gen(ld, PG)), ...
                                             Qdf(i) - sum(gen(ld, QG)));
            else
                fprintf(fd, '%10.2f%10.2f ', [ Pdf(i) Qdf(i) ]);
            end
        else
            fprintf(fd, '       -         -   ');
        end
        if isOPF
            fprintf(fd, '%9.3f', bus(i, LAM_P));
            if abs(bus(i, LAM_Q)) > ptol
                fprintf(fd, '%8.3f', bus(i, LAM_Q));
            else
                fprintf(fd, '     -');
            end
        end
    end
    fprintf(fd, '\n                        --------  --------  --------  --------');
    fprintf(fd, '\n               Total: %9.2f %9.2f %9.2f %9.2f', ...
        sum(gen(ong, PG)), sum(gen(ong, QG)), ...
        sum(Pdf(nzld)) - sum(gen(onld, PG)), ...
        sum(Qdf(nzld)) - sum(gen(onld, QG)));
    fprintf(fd, '\n');
end

%% branch data
if OUT_BRANCH && (success || OUT_FORCE)
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n|     Branch Data                                                              |');
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\nBrnch   From   To    From Bus Injection   To Bus Injection     Loss (I^2 * Z)  ');
    fprintf(fd, '\n  #     Bus    Bus    P (MW)   Q (MVAr)   P (MW)   Q (MVAr)   P (MW)   Q (MVAr)');
    fprintf(fd, '\n-----  -----  -----  --------  --------  --------  --------  --------  --------');
    fprintf(fd, '\n%4d%7d%7d%10.2f%10.2f%10.2f%10.2f%10.3f%10.2f', ...
            [   (1:nl)', branch(:, [F_BUS, T_BUS]), ...
                branch(:, [PF, QF]), branch(:, [PT, QT]), ...
                real(loss), imag(loss) ...
            ]');
    fprintf(fd, '\n                                                             --------  --------');
    fprintf(fd, '\n                                                    Total:%10.3f%10.2f', ...
            sum(real(loss)), sum(imag(loss)));
    fprintf(fd, '\n');
end

%%-----  constraint data  -----
if OUT_ANY && isOPF && (success || OUT_FORCE)
    ctol = mpopt.opf.violation; %% constraint violation tolerance
    %% voltage constraints
    if ~isDC && (OUT_V_LIM == 2 || (OUT_V_LIM == 1 && ...
                         (any(bus(:, VM) < bus(:, VMIN) + ctol) || ...
                          any(bus(:, VM) > bus(:, VMAX) - ctol) || ...
                          any(bus(:, MU_VMIN) > ptol) || ...
                          any(bus(:, MU_VMAX) > ptol))))
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n|     Voltage Constraints                                                      |');
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\nBus #  Vmin mu    Vmin    |V|   Vmax    Vmax mu');
        fprintf(fd, '\n-----  --------   -----  -----  -----   --------');
        for i = 1:nb
            if OUT_V_LIM == 2 || (OUT_V_LIM == 1 && ...
                         (bus(i, VM) < bus(i, VMIN) + ctol || ...
                          bus(i, VM) > bus(i, VMAX) - ctol || ...
                          bus(i, MU_VMIN) > ptol || bus(i, MU_VMAX) > ptol))
                fprintf(fd, '\n%5d', bus(i, BUS_I));
                if bus(i, VM) < bus(i, VMIN) + ctol || bus(i, MU_VMIN) > ptol
                    fprintf(fd, '%10.3f', bus(i, MU_VMIN));
                else
                    fprintf(fd, '      -   ');
                end
                fprintf(fd, '%8.3f%7.3f%7.3f', bus(i, [VMIN, VM, VMAX]));
                if bus(i, VM) > bus(i, VMAX) - ctol || bus(i, MU_VMAX) > ptol
                    fprintf(fd, '%10.3f', bus(i, MU_VMAX));
                else
                    fprintf(fd, '      -    ');
                end
            end
        end
        fprintf(fd, '\n');
    end

    %% generator P constraints
    if OUT_PG_LIM == 2 || ...
            (OUT_PG_LIM == 1 && (any(gen(ong, PG) < gen(ong, PMIN) + ctol) || ...
                                any(gen(ong, PG) > gen(ong, PMAX) - ctol) || ...
                                any(gen(ong, MU_PMIN) > ptol) || ...
                                any(gen(ong, MU_PMAX) > ptol))) || ...
            (~isDC && (OUT_QG_LIM == 2 || ...
            (OUT_QG_LIM == 1 && (any(gen(ong, QG) < gen(ong, QMIN) + ctol) || ...
                                any(gen(ong, QG) > gen(ong, QMAX) - ctol) || ...
                                any(gen(ong, MU_QMIN) > ptol) || ...
                                any(gen(ong, MU_QMAX) > ptol)))))
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n|     Generation Constraints                                                   |');
        fprintf(fd, '\n================================================================================');
    end
    if OUT_PG_LIM == 2 || (OUT_PG_LIM == 1 && ...
                             (any(gen(ong, PG) < gen(ong, PMIN) + ctol) || ...
                              any(gen(ong, PG) > gen(ong, PMAX) - ctol) || ...
                              any(gen(ong, MU_PMIN) > ptol) || ...
                              any(gen(ong, MU_PMAX) > ptol)))
        fprintf(fd, '\n Gen   Bus                  Active Power Limits');
        fprintf(fd, '\n  #     #     Pmin mu     Pmin       Pg       Pmax    Pmax mu');
        fprintf(fd, '\n----  -----   -------   --------  --------  --------  -------');
        for k = 1:length(ong)
            i = ong(k);
            if OUT_PG_LIM == 2 || (OUT_PG_LIM == 1 && ...
                        (gen(i, PG) < gen(i, PMIN) + ctol || ...
                         gen(i, PG) > gen(i, PMAX) - ctol || ...
                         gen(i, MU_PMIN) > ptol || gen(i, MU_PMAX) > ptol))
                fprintf(fd, '\n%4d%6d ', i, gen(i, GEN_BUS));
                if gen(i, PG) < gen(i, PMIN) + ctol || gen(i, MU_PMIN) > ptol
                    fprintf(fd, '%10.3f', gen(i, MU_PMIN));
                else
                    fprintf(fd, '      -   ');
                end
                if gen(i, PG)
                    fprintf(fd, '%10.2f%10.2f%10.2f', gen(i, [PMIN, PG, PMAX]));
                else
                    fprintf(fd, '%10.2f       -  %10.2f', gen(i, [PMIN, PMAX]));
                end
                if gen(i, PG) > gen(i, PMAX) - ctol || gen(i, MU_PMAX) > ptol
                    fprintf(fd, '%10.3f', gen(i, MU_PMAX));
                else
                    fprintf(fd, '      -   ');
                end
            end
        end
        fprintf(fd, '\n');
    end

    %% generator Q constraints
    if ~isDC && (OUT_QG_LIM == 2 || (OUT_QG_LIM == 1 && ...
                             (any(gen(ong, QG) < gen(ong, QMIN) + ctol) || ...
                              any(gen(ong, QG) > gen(ong, QMAX) - ctol) || ...
                              any(gen(ong, MU_QMIN) > ptol) || ...
                              any(gen(ong, MU_QMAX) > ptol))))
        fprintf(fd, '\n Gen   Bus                 Reactive Power Limits');
        fprintf(fd, '\n  #     #     Qmin mu     Qmin       Qg       Qmax    Qmax mu');
        fprintf(fd, '\n----  -----   -------   --------  --------  --------  -------');
        for k = 1:length(ong)
            i = ong(k);
            if OUT_QG_LIM == 2 || (OUT_QG_LIM == 1 && ...
                        (gen(i, QG) < gen(i, QMIN) + ctol || ...
                         gen(i, QG) > gen(i, QMAX) - ctol || ...
                         gen(i, MU_QMIN) > ptol || gen(i, MU_QMAX) > ptol))
                fprintf(fd, '\n%4d%6d ', i, gen(i, GEN_BUS));
                if gen(i, QG) < gen(i, QMIN) + ctol || gen(i, MU_QMIN) > ptol
                    fprintf(fd, '%10.3f', gen(i, MU_QMIN));
                else
                    fprintf(fd, '      -   ');
                end
                if gen(i, QG)
                    fprintf(fd, '%10.2f%10.2f%10.2f', gen(i, [QMIN, QG, QMAX]));
                else
                    fprintf(fd, '%10.2f       -  %10.2f', gen(i, [QMIN, QMAX]));
                end
                if gen(i, QG) > gen(i, QMAX) - ctol || gen(i, MU_QMAX) > ptol
                    fprintf(fd, '%10.3f', gen(i, MU_QMAX));
                else
                    fprintf(fd, '      -   ');
                end
            end
        end
        fprintf(fd, '\n');
    end

    %% dispatchable load P constraints
    if ~isempty(onld) && (OUT_PG_LIM == 2 || ...
            (OUT_PG_LIM == 1 && (any(gen(onld, PG) < gen(onld, PMIN) + ctol) || ...
                                any(gen(onld, PG) > gen(onld, PMAX) - ctol) || ...
                                any(gen(onld, MU_PMIN) > ptol) || ...
                                any(gen(onld, MU_PMAX) > ptol))) || ...
            (~isDC && (OUT_QG_LIM == 2 || ...
            (OUT_QG_LIM == 1 && (any(gen(onld, QG) < gen(onld, QMIN) + ctol) || ...
                                any(gen(onld, QG) > gen(onld, QMAX) - ctol) || ...
                                any(gen(onld, MU_QMIN) > ptol) || ...
                                any(gen(onld, MU_QMAX) > ptol))))))
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n|     Dispatchable Load Constraints                                            |');
        fprintf(fd, '\n================================================================================');
    end
    if ~isempty(onld) && (OUT_PG_LIM == 2 || (OUT_PG_LIM == 1 && ...
                             (any(gen(onld, PG) < gen(onld, PMIN) + ctol) || ...
                              any(gen(onld, PG) > gen(onld, PMAX) - ctol) || ...
                              any(gen(onld, MU_PMIN) > ptol) || ...
                              any(gen(onld, MU_PMAX) > ptol))))
        fprintf(fd, '\n Gen   Bus                  Active Power Limits');
        fprintf(fd, '\n  #     #     Pmin mu     Pmin       Pg       Pmax    Pmax mu');
        fprintf(fd, '\n----  -----   -------   --------  --------  --------  -------');
        for k = 1:length(onld)
            i = onld(k);
            if OUT_PG_LIM == 2 || (OUT_PG_LIM == 1 && ...
                        (gen(i, PG) < gen(i, PMIN) + ctol || ...
                         gen(i, PG) > gen(i, PMAX) - ctol || ...
                         gen(i, MU_PMIN) > ptol || gen(i, MU_PMAX) > ptol))
                fprintf(fd, '\n%4d%6d ', i, gen(i, GEN_BUS));
                if gen(i, PG) < gen(i, PMIN) + ctol || gen(i, MU_PMIN) > ptol
                    fprintf(fd, '%10.3f', gen(i, MU_PMIN));
                else
                    fprintf(fd, '      -   ');
                end
                if gen(i, PG)
                    fprintf(fd, '%10.2f%10.2f%10.2f', gen(i, [PMIN, PG, PMAX]));
                else
                    fprintf(fd, '%10.2f       -  %10.2f', gen(i, [PMIN, PMAX]));
                end
                if gen(i, PG) > gen(i, PMAX) - ctol || gen(i, MU_PMAX) > ptol
                    fprintf(fd, '%10.3f', gen(i, MU_PMAX));
                else
                    fprintf(fd, '      -   ');
                end
            end
        end
        fprintf(fd, '\n');
    end

    %% dispatchable load Q constraints
    if ~isDC && ~isempty(onld) && (OUT_QG_LIM == 2 || (OUT_QG_LIM == 1 && ...
                             (any(gen(onld, QG) < gen(onld, QMIN) + ctol) || ...
                              any(gen(onld, QG) > gen(onld, QMAX) - ctol) || ...
                              any(gen(onld, MU_QMIN) > ptol) || ...
                              any(gen(onld, MU_QMAX) > ptol))))
        fprintf(fd, '\n Gen   Bus                 Reactive Power Limits');
        fprintf(fd, '\n  #     #     Qmin mu     Qmin       Qg       Qmax    Qmax mu');
        fprintf(fd, '\n----  -----   -------   --------  --------  --------  -------');
        for k = 1:length(onld)
            i = onld(k);
            if OUT_QG_LIM == 2 || (OUT_QG_LIM == 1 && ...
                        (gen(i, QG) < gen(i, QMIN) + ctol || ...
                         gen(i, QG) > gen(i, QMAX) - ctol || ...
                         gen(i, MU_QMIN) > ptol || gen(i, MU_QMAX) > ptol))
                fprintf(fd, '\n%4d%6d ', i, gen(i, GEN_BUS));
                if gen(i, QG) < gen(i, QMIN) + ctol || gen(i, MU_QMIN) > ptol
                    fprintf(fd, '%10.3f', gen(i, MU_QMIN));
                else
                    fprintf(fd, '      -   ');
                end
                if gen(i, QG)
                    fprintf(fd, '%10.2f%10.2f%10.2f', gen(i, [QMIN, QG, QMAX]));
                else
                    fprintf(fd, '%10.2f       -  %10.2f', gen(i, [QMIN, QMAX]));
                end
                if gen(i, QG) > gen(i, QMAX) - ctol || gen(i, MU_QMAX) > ptol
                    fprintf(fd, '%10.3f', gen(i, MU_QMAX));
                else
                    fprintf(fd, '      -   ');
                end
            end
        end
        fprintf(fd, '\n');
    end

    %% line flow constraints
    lim_type = upper(mpopt.opf.flow_lim(1));
    if isDC || lim_type == 'P' || lim_type == '2'   %% |P| limit
        Ff = branch(:, PF);
        Ft = branch(:, PT);
        str = '\n  #     Bus    Pf  mu     Pf      |Pmax|      Pt      Pt  mu   Bus';
    elseif lim_type == 'I'                          %% |I| limit
        Ff = abs( (branch(:, PF) + 1j * branch(:, QF)) ./ V(e2i(branch(:, F_BUS))) );
        Ft = abs( (branch(:, PT) + 1j * branch(:, QT)) ./ V(e2i(branch(:, T_BUS))) );
        str = '\n  #     Bus   |If| mu    |If|     |Imax|     |It|    |It| mu   Bus';
    else                                            %% |S| limit
        Ff = abs(branch(:, PF) + 1j * branch(:, QF));
        Ft = abs(branch(:, PT) + 1j * branch(:, QT));
        str = '\n  #     Bus   |Sf| mu    |Sf|     |Smax|     |St|    |St| mu   Bus';
    end
    if any(branch(:, RATE_A) ~= 0) && (OUT_LINE_LIM == 2 || (OUT_LINE_LIM == 1 && ...
                        (any(abs(Ff) > branch(:, RATE_A) - ctol) || ...
                         any(abs(Ft) > branch(:, RATE_A) - ctol) || ...
                         any(branch(:, MU_SF) > ptol) || ...
                         any(branch(:, MU_ST) > ptol))))
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n|     Branch Flow Constraints                                                  |');
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\nBrnch   From     "From" End        Limit       "To" End        To');
        fprintf(fd, str);
        fprintf(fd, '\n-----  -----  -------  --------  --------  --------  -------  -----');
        for i = 1:nl
            if branch(i, RATE_A) ~= 0 && (OUT_LINE_LIM == 2 || (OUT_LINE_LIM == 1 && ...
                   (abs(Ff(i)) > branch(i, RATE_A) - ctol || ...
                    abs(Ft(i)) > branch(i, RATE_A) - ctol || ...
                    branch(i, MU_SF) > ptol || branch(i, MU_ST) > ptol)))
                fprintf(fd, '\n%4d%7d', i, branch(i, F_BUS));
                if Ff(i) > branch(i, RATE_A) - ctol || branch(i, MU_SF) > ptol
                    fprintf(fd, '%10.3f', branch(i, MU_SF));
                else
                    fprintf(fd, '      -   ');
                end
                fprintf(fd, '%9.2f%10.2f%10.2f', ...
                    [Ff(i), branch(i, RATE_A), Ft(i)]);
                if Ft(i) > branch(i, RATE_A) - ctol || branch(i, MU_ST) > ptol
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

%% execute userfcn callbacks for 'printpf' stage
if have_results_struct && isfield(results, 'userfcn') && (success || OUT_FORCE)
    if ~isOPF   %% turn off option for all constraints if it isn't an OPF
        mpopt = mpoption(mpopt, 'out.lim.all', 0);
    end
    run_userfcn(results.userfcn, 'printpf', results, fd, mpopt);
end
if OUT_ANY && ~success
    if OUT_FORCE
        if isSDP
            fprintf(fd, '\n>>>>>  Solution does NOT satisfy rank and/or consistency conditions (%.2f seconds).  <<<<<\nmineigratio = %0.5g, zero_eval = %0.5g\n', et, mineigratio, zero_eval);
        else
            fprintf(fd, '\n>>>>>  Did NOT converge (%.2f seconds)  <<<<<\n', et);
        end
    end
    fprintf('\n');
end
