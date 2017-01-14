function [groupss, isolated] = case_info(mpc, fd)
%CASE_INFO Prints information about islands in a network.
%   CASE_INFO(MPC)
%   CASE_INFO(MPC, FD)
%   [GROUPS, ISOLATED] = CASE_INFO(...)
%
%   Prints out detailed information about a MATPOWER case. Optionally prints
%   to an open file, whose file identifier, as returned by FOPEN, is
%   specified in the optional second parameter FD. Optional return arguments
%   include GROUPS and ISOLATED buses, as returned by FIND_ISLANDS.

%   TO DO: Add checking of bus types (isolated bus not marked as such).
%          Check for infeasible limits PMAX < PMIN, QMAX < QMIN, VMAX < VMIN,
%               etc.
%          Warn about islands without reference buses.
%          Warn about PV and ref buses without generators.
%          Separate branch charging injections from series losses.
%          Report islands that are connected by DC lines.

%   MATPOWER
%   Copyright (c) 2012-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
c = idx_dcline;

if nargin < 2
    fd = 1;     %% print to stdio by default
end

tic;

%% load case if necessary
mpc = loadcase(mpc);

%% sizes of things
nb  = size(mpc.bus, 1);     %% number of buses
nl  = size(mpc.branch, 1);  %% number of branches
if isfield(mpc, 'dcline')   %% number of DC lines
    ndc = size(mpc.dcline, 1);
else
    ndc = 0;
end
ng  = size(mpc.gen, 1);     %% number of dispatchable injections

%% max bus number
mb  = max(abs([mpc.bus(:, BUS_I); mpc.gen(:, GEN_BUS); ...
            mpc.branch(:, F_BUS); mpc.branch(:, T_BUS)]));
nbase = 10^(fix(log10(mb))+1);  %% negative bus numbers
                                %% are subtracted from nbase
                                %% to get something positive

%% check for negative or zero bus numbers in bus matrix
bus_i = mpc.bus(:, BUS_I);
nonpos_bus = find(bus_i <= 0);
bus_i(nonpos_bus) = -bus_i(nonpos_bus) + nbase;

%% create external to internal bus map, with range for
%% converted non-positive bus numbers
e2i = sparse(bus_i, ones(nb, 1), 1:nb, nbase+mb, 1);

%% find unknown bus numbers from gen, branch or dcline matrices
unknown_gbus = unknown_buses(e2i, nbase, mpc.gen(:, GEN_BUS));
unknown_fbus = unknown_buses(e2i, nbase, mpc.branch(:, F_BUS));
unknown_tbus = unknown_buses(e2i, nbase, mpc.branch(:, T_BUS));
if ndc
    unknown_fbusdc = unknown_buses(e2i, nbase, mpc.dcline(:, c.F_BUS));
    unknown_tbusdc = unknown_buses(e2i, nbase, mpc.dcline(:, c.T_BUS));
end

if length(nonpos_bus)
    fprintf(fd, 'Bad bus numbers:              %d\n', length(nonpos_bus));
    for k = 1:length(nonpos_bus)
        s = sprintf('bus(%d, BUS_I)', nonpos_bus(k));
        fprintf(fd, '%24s = %d\n', s, mpc.bus(nonpos_bus(k), BUS_I));
    end
end
if ~isempty(unknown_gbus)
    fprintf(fd, 'Unknown generator buses:      %d\n', length(unknown_gbus));
    for k = 1:length(unknown_gbus)
        s = sprintf('gen(%d, GEN_BUS)', unknown_gbus(k));
        fprintf(fd, '%24s = %d\n', s, mpc.gen(unknown_gbus(k), GEN_BUS));
    end

    %% remove them
    mpc.gen(unknown_gbus, :) = [];
    ng  = size(mpc.gen, 1);     %% number of dispatchable injections
end
if ~isempty(unknown_fbus)
    fprintf(fd, 'Unknown branch "from" buses:  %d\n', length(unknown_fbus));
    for k = 1:length(unknown_fbus)
        s = sprintf('branch(%d, F_BUS)', unknown_fbus(k));
        fprintf(fd, '%24s = %d\n', s, mpc.branch(unknown_fbus(k), F_BUS));
    end
end
if ~isempty(unknown_tbus)
    fprintf(fd, 'Unknown branch "to" buses:    %d\n', length(unknown_tbus));
    for k = 1:length(unknown_tbus)
        s = sprintf('branch(%d, T_BUS)', unknown_tbus(k));
        fprintf(fd, '%24s = %d\n', s, mpc.branch(unknown_tbus(k), T_BUS));
    end
end
%% remove branches connected to unknown buses
if ~isempty(unknown_fbus) || ~isempty(unknown_tbus)
    tmp = unique([unknown_fbus; unknown_tbus]);
    mpc.branch(tmp, :) = [];    %% remove them
    nl  = size(mpc.branch, 1);  %% number of branches
end
if ndc
    if ~isempty(unknown_fbusdc)
        fprintf(fd, 'Unknown DC line "from" buses: %d\n', length(unknown_fbusdc));
        for k = 1:length(unknown_fbusdc)
            s = sprintf('dcline(%d, c.F_BUS)', unknown_fbusdc(k));
            fprintf(fd, '%24s = %d\n', s, mpc.dcline(unknown_fbusdc(k), c.F_BUS));
        end
    end
    if ~isempty(unknown_tbusdc)
        fprintf(fd, 'Unknown DC line "to" buses:   %d\n', length(unknown_tbusdc));
        for k = 1:length(unknown_tbusdc)
            s = sprintf('dcline(%d, c.T_BUS)', unknown_tbusdc(k));
            fprintf(fd, '%24s = %d\n', s, mpc.dcline(unknown_tbusdc(k), c.T_BUS));
        end
    end
    %% remove branches connected to unknown buses
    if ~isempty(unknown_fbusdc) || ~isempty(unknown_tbusdc)
        tmp = unique([unknown_fbusdc; unknown_tbusdc]);
        mpc.dcline(tmp, :) = [];    %% remove them
        ndc  = size(mpc.dcline, 1); %% number of DC lines
    end
end

if length(nonpos_bus) == 0  %% no bad bus numbers

%% build connectivity structures
C_on = sparse(1:nl, e2i(mpc.branch(:, F_BUS)), -mpc.branch(:, BR_STATUS), nl, nb) + ...
    sparse(1:nl, e2i(mpc.branch(:, T_BUS)),  mpc.branch(:, BR_STATUS), nl, nb);
C = sparse(1:nl, e2i(mpc.branch(:, F_BUS)), -1, nl, nb) + ...
    sparse(1:nl, e2i(mpc.branch(:, T_BUS)),  1, nl, nb);
if ndc
    Cdc_on = sparse(1:ndc, e2i(mpc.dcline(:, c.F_BUS)), -mpc.dcline(:, c.BR_STATUS), ndc, nb) + ...
        sparse(1:ndc, e2i(mpc.dcline(:, c.T_BUS)),  mpc.dcline(:, c.BR_STATUS), ndc, nb);
    Cdc = sparse(1:ndc, e2i(mpc.dcline(:, c.F_BUS)), -1, ndc, nb) + ...
        sparse(1:ndc, e2i(mpc.dcline(:, c.T_BUS)),  1, ndc, nb);
end
Cg_on = sparse(1:ng, e2i(mpc.gen(:, GEN_BUS)), mpc.gen(:, GEN_STATUS), ng, nb);
Cg = sparse(1:ng, e2i(mpc.gen(:, GEN_BUS)), 1, ng, nb);

%% check for islands
fprintf(fd, 'Checking connectivity ... ');
[groups, isolated] = connected_components(C_on);

ngr = length(groups);   %% number of islands
nis = length(isolated); %% number of isolated buses
have_isolated = nis > 0;
if ngr == 1
    if have_isolated
        if nis == 1, s = ''; else, s = 'es'; end
        fprintf(fd, 'single connected network, plus %d isolated bus%s\n', nis, s);
    else
        fprintf(fd, 'single fully connected network\n');
    end
else
    if nis == 1, s = ''; else, s = 'es'; end
    fprintf(fd, '%d connected groups, %d isolated bus%s\n', ngr, nis, s);
end

%% collect info on groups
bron  = mpc.branch(:, BR_STATUS)  > 0;
broff = mpc.branch(:, BR_STATUS) <= 0;
if ndc
    dcon  = mpc.dcline(:, c.BR_STATUS)  > 0;
    dcoff = mpc.dcline(:, c.BR_STATUS) <= 0;
end
gon   = mpc.gen(:, GEN_STATUS)  > 0;
goff  = mpc.gen(:, GEN_STATUS) <= 0;

%% initialize data structures
d0 = struct( ...
    'nb', 0, ...
    'nl', 0, ...
    'nl_on', 0, ...
    'nl_off', 0, ...
    'nlt', 0, ...
...%     'nlt_on', 0, ... %% (always zero)
    'nlt_off', 0, ...
    'ndc', 0, ...
    'ndc_on', 0, ...
    'ndc_off', 0, ...
    'ndct', 0, ...
    'ndct_on', 0, ...
    'ndct_off', 0, ...
    'ndc_all', 0, ...
    'ng', 0, ...
    'ng_on', 0, ...
    'ng_off', 0, ...
    'nsh', 0, ...
    'nld', 0, ...
    'nld_on', 0, ...
    'nld_off', 0, ...
    'nfld', 0, ...
    'ndld', 0, ...
    'ndld_on', 0, ...
    'ndld_off', 0, ...
    'Pmax', 0, ...
    'Qmax', 0, ...
    'Pmax_on', 0, ...
    'Qmax_on', 0, ...
    'Pmax_off', 0, ...
    'Qmax_off', 0, ...
    'Pmin', 0, ...
    'Qmin', 0, ...
    'Pmin_on', 0, ...
    'Qmin_on', 0, ...
    'Pmin_off', 0, ...
    'Qmin_off', 0, ...
    'Pg', 0, ...
    'Qg', 0, ...
    'Ps', 0, ...
    'Qs', 0, ...
    'Ploss', 0, ...
    'Qloss', 0, ...
    'Pd', 0, ...
    'Qd', 0, ...
    'Pd_fixed', 0, ...
    'Qd_fixed', 0, ...
    'Pd_disp_cap', 0, ...
    'Qd_disp_cap', 0, ...
    'Pd_disp_cap_on', 0, ...
    'Qd_disp_cap_on', 0, ...
    'Pd_disp_cap_off', 0, ...
    'Qd_disp_cap_off', 0, ...
    'Pd_disp', 0, ...
    'Qd_disp', 0, ...
    'Pd_curtailed', 0, ...
    'Qd_curtailed', 0, ...
    'Pd_cap', 0, ...
    'Qd_cap', 0, ...
    'Pd_cap_on', 0, ...
    'Qd_cap_on', 0, ...
    'Pd_cap_off', 0, ...
    'Qd_cap_off', 0, ...
    'Pdc', 0, ...
    'Pmaxdc', 0, ...
    'Pmaxdc_on', 0, ...
    'Pmaxdc_off', 0, ...
    'Pmindc', 0, ...
    'Pmindc_on', 0, ...
    'Pmindc_off', 0 ...
);
d = d0;
total = d0;
allrefs = find(mpc.bus(:, BUS_TYPE) == REF);
refs = {};
nrefs = 0;
ibr_tie_all = [];
idc_tie_all = [];
idc_tie_all_on  = [];
idc_tie_all_off = [];

%% gather data
fields = fieldnames(d);
for k = 1:ngr+have_isolated
    %% initialize d(k)
    for f = 1:length(fields)
        ff = fields{f};
        d(k).(ff) = d0.(ff);
    end
    if k > ngr
        b = isolated';
        ibr = [];
        idc = [];
    else
        b = groups{k};
        %% branches with both ends in group
        ibr = find(sum(abs(C(:, b)), 2) & ~sum(C(:, b), 2));
        if ndc
            idc = find(sum(abs(Cdc(:, b)), 2) & ~sum(Cdc(:, b), 2));
        else
            idc = [];
        end
    end
    %% branches with one end in group
    ibr_tie = find(sum(C(:, b), 2));
    if ndc
        idc_tie = find(sum(Cdc(:, b), 2));
    else
        idc_tie = [];
    end
    refs{k} = b(find(mpc.bus(b, BUS_TYPE) == REF));
    nrefs = nrefs + length(refs{k});

    ibr_on  = ibr(find(bron(ibr)));
    ibr_off = ibr(find(broff(ibr)));
%     ibr_tie_on  = ibr_tie(find(bron(ibr_tie))); %% (always empty)
    ibr_tie_off = ibr_tie(find(broff(ibr_tie)));
    ibr_tie_all = unique([ ibr_tie_all; ibr_tie_off ]);
    if ndc
        idc_on  = idc(find(dcon(idc)));
        idc_off = idc(find(dcoff(idc)));
        idc_tie_on  = idc_tie(find(dcon(idc_tie)));
        idc_tie_off = idc_tie(find(dcoff(idc_tie)));
        idc_tie_all = unique([ idc_tie_all; idc_tie ]);
        idc_tie_all_on  = unique([ idc_tie_all_on;  idc_tie_on ]);
        idc_tie_all_off = unique([ idc_tie_all_off; idc_tie_off ]);
    end
    ig      = find(sum(abs(Cg(:, b)), 2));
    if isempty(ig)
        ig_on   = [];
        ig_off  = [];
        idld_on = [];
        idld_off= [];
    else
        ig_on   = ig(find(gon(ig)  & ~isload(mpc.gen(ig, :))));
        ig_off  = ig(find(goff(ig) & ~isload(mpc.gen(ig, :))));
        idld_on = ig(find(gon(ig)  &  isload(mpc.gen(ig, :))));
        idld_off= ig(find(goff(ig) &  isload(mpc.gen(ig, :))));
    end

    d(k).nb = length(b);        %% # of buses
    d(k).nl = length(ibr);      %% # of branches
    d(k).nl_on  = length(ibr_on);   %% # of in-service branches
    d(k).nl_off = length(ibr_off);  %% # of out-of-service branches
    d(k).nlt = length(ibr_tie);     %% # of tie-lines
%     d(k).nlt_on  = length(ibr_tie_on);  %% # of in-service tie-lines (always zero)
    d(k).nlt_off = length(ibr_tie_off); %% # of out-of-service tie-lines
    if ndc
        d(k).ndc = length(idc);      %% # of dc lines
        d(k).ndc_on  = length(idc_on);   %% # of in-service dc lines
        d(k).ndc_off = length(idc_off);  %% # of out-of-service dc lines
        d(k).ndct = length(idc_tie);      %% # of dc tie-lines
        d(k).ndct_on  = length(idc_tie_on);   %% # of in-service dc tie-lines
        d(k).ndct_off = length(idc_tie_off);  %% # of out-of-service dc tie-lines
        d(k).ndc_all = d(k).ndc + d(k).ndct;
    end
    d(k).ng = length(ig_on) + length(ig_off);   %% # of gen
    d(k).ng_on  = length(ig_on);    %% # of in-service gens
    d(k).ng_off = length(ig_off);   %% # of out-of-service gens
    d(k).nsh = length(find(mpc.bus(b, GS) ~= 0 | mpc.bus(b, BS) ~= 0)); %% # of shunt elements
    d(k).nfld = length(find(mpc.bus(b, PD) ~= 0 | mpc.bus(b, QD) ~= 0));    %% # of fixed loads
    d(k).ndld = length(idld_on) + length(idld_off); %% # of disp loads
    d(k).ndld_on  = length(idld_on);    %% # of in-service disp loads
    d(k).ndld_off = length(idld_off);   %% # of out-of-service disp loads
    d(k).nld = d(k).nfld + d(k).ndld;   %% # of fixed + disp loads
    d(k).nld_on  = d(k).nfld + d(k).ndld_on;    %% # of in-service fixed + disp loads
    d(k).nld_off = d(k).ndld_off;       %% # of out-of-service fixed + disp loads

    d(k).Pmax_on    = sum(mpc.gen(ig_on,  PMAX));
    d(k).Pmax_off   = sum(mpc.gen(ig_off, PMAX));
    d(k).Pmax       = d(k).Pmax_on + d(k).Pmax_off;
    d(k).Pmin_on    = sum(mpc.gen(ig_on,  PMIN));
    d(k).Pmin_off   = sum(mpc.gen(ig_off, PMIN));
    d(k).Pmin       = d(k).Pmin_on + d(k).Pmin_off;
    d(k).Pg         = sum(mpc.gen(ig_on,  PG));

    d(k).Qmax_on    = sum(mpc.gen(ig_on,  QMAX));
    d(k).Qmax_off   = sum(mpc.gen(ig_off, QMAX));
    d(k).Qmax       = d(k).Qmax_on + d(k).Qmax_off;
    d(k).Qmin_on    = sum(mpc.gen(ig_on,  QMIN));
    d(k).Qmin_off   = sum(mpc.gen(ig_off, QMIN));
    d(k).Qmin       = d(k).Qmin_on + d(k).Qmin_off;
    d(k).Qg         = sum(mpc.gen(ig_on,  QG));

    d(k).Ps         = sum(-mpc.bus(b, VM) .^ 2 .* mpc.bus(b, GS));
    d(k).Qs         = sum( mpc.bus(b, VM) .^ 2 .* mpc.bus(b, BS));
    if size(mpc.branch, 2) > PF
        d(k).Ploss      = sum(mpc.branch(ibr_on, PF) + mpc.branch(ibr_on, PT));
        d(k).Qloss      = sum(mpc.branch(ibr_on, QF) + mpc.branch(ibr_on, QT));
    end

    d(k).Pd_fixed   = sum(mpc.bus(b, PD));
    d(k).Qd_fixed   = sum(mpc.bus(b, QD));
    d(k).Pd_disp_cap_on     = sum(-mpc.gen(idld_on,  PMIN));
    d(k).Qd_disp_cap_on     = sum(-mpc.gen(idld_on,  QMIN));
    d(k).Pd_disp_cap_off    = sum(-mpc.gen(idld_off, PMIN));
    d(k).Qd_disp_cap_off    = sum(-mpc.gen(idld_off, QMIN));
    d(k).Pd_disp_cap        = d(k).Pd_disp_cap_on + d(k).Pd_disp_cap_off;
    d(k).Qd_disp_cap        = d(k).Qd_disp_cap_on + d(k).Qd_disp_cap_off;
    d(k).Pd_disp            = sum(-mpc.gen(idld_on,  PG));
    d(k).Qd_disp            = sum(-mpc.gen(idld_on,  QG));
    d(k).Pd_curtailed       = d(k).Pd_disp_cap_on - d(k).Pd_disp;
    d(k).Qd_curtailed       = d(k).Qd_disp_cap_on - d(k).Qd_disp;
    d(k).Pd                 = d(k).Pd_fixed + d(k).Pd_disp;
    d(k).Qd                 = d(k).Qd_fixed + d(k).Qd_disp;
    d(k).Pd_cap             = d(k).Pd_fixed + d(k).Pd_disp_cap;
    d(k).Qd_cap             = d(k).Qd_fixed + d(k).Qd_disp_cap;
    d(k).Pd_cap_on          = d(k).Pd_fixed + d(k).Pd_disp_cap_on;
    d(k).Qd_cap_on          = d(k).Qd_fixed + d(k).Qd_disp_cap_on;
    d(k).Pd_cap_off         = d(k).Pd_disp_cap_off;
    d(k).Qd_cap_off         = d(k).Qd_disp_cap_off;

    if ndc
        f = find(sum(Cdc(:, b), 2) < 0);
        t = find(sum(Cdc(:, b), 2) > 0);
        d(k).Pdc    = sum(mpc.dcline(f, c.PF)) - sum(mpc.dcline(t, c.PT));
        d(k).Pmaxdc = sum(mpc.dcline(f, c.PMAX)) - sum(mpc.dcline(t, c.PMAX));
        d(k).Pmindc = sum(mpc.dcline(f, c.PMIN)) - sum(mpc.dcline(t, c.PMIN));
        f = find(sum(Cdc(:, b), 2) < 0 & dcon);
        t = find(sum(Cdc(:, b), 2) > 0 & dcon);
        d(k).Pmaxdc_on = sum(mpc.dcline(f, c.PMAX)) - sum(mpc.dcline(t, c.PMAX));
        d(k).Pmindc_on = sum(mpc.dcline(f, c.PMIN)) - sum(mpc.dcline(t, c.PMIN));
        f = find(sum(Cdc(:, b), 2) < 0 & dcoff);
        t = find(sum(Cdc(:, b), 2) > 0 & dcoff);
        d(k).Pmaxdc_off = sum(mpc.dcline(f, c.PMAX)) - sum(mpc.dcline(t, c.PMAX));
        d(k).Pmindc_off = sum(mpc.dcline(f, c.PMIN)) - sum(mpc.dcline(t, c.PMIN));
    end

    %% accumulate totals
    for f = 1:length(fields)
        ff = fields{f};
        total.(ff) = total.(ff) + d(k).(ff);
    end

    total.nl = nl;
    total.nl_on  = length(find(bron));
    total.nl_off = length(find(broff));
    total.ndc = ndc;
    if ndc
        total.ndc_on  = length(find(dcon));
        total.ndc_off = length(find(dcoff));
    else
        total.ndc_on  = 0;
        total.ndc_off = 0;
    end
end
total.nlt       = length(ibr_tie_all);
% total.nlt_on  = 0;
total.nlt_off   = length(ibr_tie_all);
if ndc
    total.ndct      = length(idc_tie_all);
    total.ndct_on   = length(idc_tie_all_on);
    total.ndct_off  = length(idc_tie_all_off);
end

%% print summary
et = toc;
fprintf(fd, 'Elapsed time is %f seconds.\n', et);
fprintf(fd, '================================================================================\n');
pages = ceil((ngr + have_isolated + 1) / 5);
for page = 1:pages
    if page == 1
        if ngr == 1 && ~have_isolated
            islands = [];
        else
            islands = 1:min(4, ngr+have_isolated);
        end
    else
        fprintf(fd, '--------------------------------------------------------------------------------\n');
        islands = (5*(page-1)):min(5*page-1, ngr+have_isolated);
    end

    %% header row 1
    fprintf(fd, '%-20s', '');
    if page == 1
        fprintf(fd, '    Full    ');
    end
    for k = islands
        if k > ngr
            fprintf(fd, '  Isolated  ');
        else
            fprintf(fd, '   Island   ');
        end
    end
    fprintf(fd, '\n');

    %% header row 2
    fprintf(fd, '%-20s', '');
    if page == 1
        fprintf(fd, '   System   ');
    end
    for k = islands
        if k > ngr
            fprintf(fd, '    Buses   ');
        else
            fprintf(fd, '  %5d     ', k);
        end
    end
    fprintf(fd, '\n');

    %% header row 3
    fprintf(fd, '%-20s', 'Number of:');
    if page == 1
        fprintf(fd, ' ---------- ');
    end
    for k = islands
        fprintf(fd, ' ---------- ');
    end
    fprintf(fd, '\n');

    p = struct('page', page, 'islands', islands, 'total', total, 'd', d);

    print_row(fd, p, ' %8d   ', '  buses',          'nb');
    print_row(fd, p, ' %8d   ', '  loads',          'nld');
    print_row(fd, p, ' %8d   ', '    on',           'nld_on');
    print_row(fd, p, ' %8d   ', '    off',          'nld_off');
    print_row(fd, p, ' %8d   ', '    fixed',        'nfld');
    print_row(fd, p, ' %8d   ', '    dispatchable', 'ndld');
    print_row(fd, p, ' %8d   ', '      on',         'ndld_on');
    print_row(fd, p, ' %8d   ', '      off',        'ndld_off');
    print_row(fd, p, ' %8d   ', '  generators',     'ng');
    print_row(fd, p, ' %8d   ', '    on',           'ng_on');
    print_row(fd, p, ' %8d   ', '    off',          'ng_off');
    print_row(fd, p, ' %8d   ', '  shunt elements', 'nsh');
    print_row(fd, p, ' %8d   ', '  branches',       'nl');
    print_row(fd, p, ' %8d   ', '    on',           'nl_on');
    print_row(fd, p, ' %8d   ', '    off',          'nl_off');
%     print_row(fd, p, ' %8d   ', '    ties',         'nlt');     %% (always same as nlt_off)
    print_row(fd, p, ' %8d   ', '    ties (off)',   'nlt_off');
    if ndc
        print_row(fd, p, ' %8d   ', '  DC lines',     'ndc_all');
        print_row(fd, p, ' %8d   ', '    within',     'ndc');
        print_row(fd, p, ' %8d   ', '      on',       'ndc_on');
        print_row(fd, p, ' %8d   ', '      off',      'ndc_off');
        print_row(fd, p, ' %8d   ', '    ties',       'ndct');
        print_row(fd, p, ' %8d   ', '      on',       'ndct_on');
        print_row(fd, p, ' %8d   ', '      off',      'ndct_off');
    end

    fprintf(fd, '\n%-20s\n', 'Load');
    fprintf(fd,   '%-20s\n', '  active (MW)');
    print_row(fd, p, '%11.1f ', '    dispatched',       'Pd');
    print_row(fd, p, '%11.1f ', '      fixed',          'Pd_fixed');
    print_row(fd, p, '%11.1f ', '      dispatchable',   'Pd_disp');
    print_row(fd, p, '%11.1f ', '    curtailed',        'Pd_curtailed');
    print_row(fd, p, '%11.1f ', '    nominal',          'Pd_cap');
    print_row(fd, p, '%11.1f ', '      on',             'Pd_cap_on');
    print_row(fd, p, '%11.1f ', '      off',            'Pd_cap_off');
    print_row(fd, p, '%11.1f ', '      fixed',          'Pd_fixed');
    print_row(fd, p, '%11.1f ', '      dispatchable',   'Pd_disp_cap');
    print_row(fd, p, '%11.1f ', '        on',           'Pd_disp_cap_on');
    print_row(fd, p, '%11.1f ', '        off',          'Pd_disp_cap_off');
    fprintf(fd,   '%-20s\n', '  reactive (MVAr)');
    print_row(fd, p, '%11.1f ', '    dispatched',       'Qd');
    print_row(fd, p, '%11.1f ', '      fixed',          'Qd_fixed');
    print_row(fd, p, '%11.1f ', '      dispatchable',   'Qd_disp');
    print_row(fd, p, '%11.1f ', '    curtailed',        'Qd_curtailed');
    print_row(fd, p, '%11.1f ', '    nominal',          'Qd_cap');
    print_row(fd, p, '%11.1f ', '      on',             'Qd_cap_on');
    print_row(fd, p, '%11.1f ', '      off',            'Qd_cap_off');
    print_row(fd, p, '%11.1f ', '      fixed',          'Qd_fixed');
    print_row(fd, p, '%11.1f ', '      dispatchable',   'Qd_disp_cap');
    print_row(fd, p, '%11.1f ', '        on',           'Qd_disp_cap_on');
    print_row(fd, p, '%11.1f ', '        off',          'Qd_disp_cap_off');

    fprintf(fd, '\n%-20s\n', 'Generation');
    fprintf(fd,   '%-20s\n', '  active (MW)');
    print_row(fd, p, '%11.1f ', '    dispatched',       'Pg');
    print_row(fd, p, '%11.1f ', '    max capacity',     'Pmax');
    print_row(fd, p, '%11.1f ', '      on',             'Pmax_on');
    print_row(fd, p, '%11.1f ', '      off',            'Pmax_off');
    print_row(fd, p, '%11.1f ', '    min capacity',     'Pmin');
    print_row(fd, p, '%11.1f ', '      on',             'Pmin_on');
    print_row(fd, p, '%11.1f ', '      off',            'Pmin_off');
    fprintf(fd,   '%-20s\n', '  reactive (MVAr)');
    print_row(fd, p, '%11.1f ', '    dispatched',       'Qg');
    print_row(fd, p, '%11.1f ', '    max capacity',     'Qmax');
    print_row(fd, p, '%11.1f ', '      on',             'Qmax_on');
    print_row(fd, p, '%11.1f ', '      off',            'Qmax_off');
    print_row(fd, p, '%11.1f ', '    min capacity',     'Qmin');
    print_row(fd, p, '%11.1f ', '      on',             'Qmin_on');
    print_row(fd, p, '%11.1f ', '      off',            'Qmin_off');

    fprintf(fd, '\n%-20s\n', 'Shunt Injections');
    print_row(fd, p, '%11.1f ', '    active (MW)',      'Ps');
    print_row(fd, p, '%11.1f ', '    reactive (MVAr)',  'Qs');

    fprintf(fd, '\n%-20s\n', 'Branch Losses');
    print_row(fd, p, '%11.1f ', '    active (MW)',      'Ploss');
    print_row(fd, p, '%11.1f ', '    reactive (MVAr)',  'Qloss');

    fprintf(fd, '\n%-20s\n', 'DC line');
    fprintf(fd,   '%-20s\n', '  export (MW)');
    print_row(fd, p, '%11.1f ', '    dispatch',      'Pdc');
    print_row(fd, p, '%11.1f ', '    max capacity',  'Pmaxdc');
    print_row(fd, p, '%11.1f ', '      on',          'Pmaxdc_on');
    print_row(fd, p, '%11.1f ', '      off',         'Pmaxdc_off');
    print_row(fd, p, '%11.1f ', '    min capacity',  'Pmindc');
    print_row(fd, p, '%11.1f ', '      on',          'Pmindc_on');
    print_row(fd, p, '%11.1f ', '      off',         'Pmindc_off');

    fprintf(fd, '\n%-20s\n', 'Reference Buses');

    fprintf(fd, '%-20s', '  num of ref buses');
    if page == 1
        fprintf(fd, ' %8d   ', nrefs);
    end
    for k = islands
        fprintf(fd, ' %8d   ', length(refs{k}));
    end
    fprintf(fd, '\n');

    for j = 1:nrefs
        if j == 1
            fprintf(fd, '%-20s', '  ref bus numbers');
        else
            fprintf(fd, '%-20s', '');
        end
        if page == 1
            fprintf(fd, ' %8d   ', mpc.bus(allrefs(j), BUS_I));
        end
        for k = islands
            if j <= length(refs{k})
                fprintf(fd, ' %8d   ', mpc.bus(refs{k}(j), BUS_I));
            else
                fprintf(fd, ' %8s   ', '');
            end
        end
        fprintf(fd, '\n');
    end

    if page ~= pages
        fprintf(fd, '\n\n');
    end
end

end     %% no bad bus numbers

%% assign output arguments as requested
if nargout > 0
    groupss = groups;
end


function print_row(fd, p, template, name, field)
templatez = sprintf('%%%ds', length(sprintf(template, 0)));
fprintf(fd, '%-20s', name);
if p.page == 1
    if p.total.(field) == 0
        fprintf(fd, templatez, '-   ');
    else
        fprintf(fd, template, p.total.(field));
    end
end
for k = p.islands
    if p.d(k).(field) == 0
        fprintf(fd, templatez, '-   ');
    else
        fprintf(fd, template, p.d(k).(field));
    end
end
fprintf(fd, '\n');


function unknown = unknown_buses(e2i, nbase, bus_list)
nonpos = find(bus_list <= 0);
bus_list(nonpos) = -bus_list(nonpos) + nbase;
unknown = find(e2i(bus_list) == 0);
