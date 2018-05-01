function mpc = toggle_dcline(mpc, on_off)
%TOGGLE_DCLINE Enable, disable or check status of DC line modeling.
%   MPC = TOGGLE_DCLINE(MPC, 'on')
%   MPC = TOGGLE_DCLINE(MPC, 'off')
%   T_F = TOGGLE_DCLINE(MPC, 'status')
%
%   Enables, disables or checks the status of a set of OPF userfcn
%   callbacks to implement DC lines as a pair of linked generators.
%   While it uses the OPF extension mechanism, this implementation
%   works for simple power flow as well as OPF problems.
%
%   These callbacks expect to find a 'dcline' field in the input MPC,
%   where MPC.dcline is an ndc x 17 matrix with columns as defined
%   in IDX_DCLINE, where ndc is the number of DC lines.
%
%   The 'int2ext' callback also packages up flow results and stores them
%   in appropriate columns of MPC.dcline.
%
%   NOTE: Because of the way this extension modifies the number of
%   rows in the gen and gencost matrices, caution must be taken
%   when using it with other extensions that deal with generators.
%
%   Examples:
%       mpc = loadcase('t_case9_dcline');
%       mpc = toggle_dcline(mpc, 'on');
%       results1 = runpf(mpc);
%       results2 = runopf(mpc);
%
%   See also IDX_DCLINE, ADD_USERFCN, REMOVE_USERFCN, RUN_USERFCN.

%   MATPOWER
%   Copyright (c) 2011-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if strcmp(upper(on_off), 'ON')
    %% define named indices into data matrices
    c = idx_dcline;

    %% check for proper input data
    if ~isfield(mpc, 'dcline') || size(mpc.dcline, 2) < c.LOSS1
        error('toggle_dcline: case must contain a ''dcline'' field, an ndc x %d matrix.', c.LOSS1);
    end
    if isfield(mpc, 'dclinecost') && size(mpc.dcline, 1) ~= size(mpc.dclinecost, 1)
        error('toggle_dcline: number of rows in ''dcline'' field (%d) and ''dclinecost'' field (%d) do not match.', ...
            size(mpc.dcline, 1), size(mpc.dclinecost, 1));
    end
    l0 = mpc.dcline(:, c.LOSS0);
    l1 = mpc.dcline(:, c.LOSS1);
    k = find( (l0 + l1 .* mpc.dcline(:, c.PMIN) < 0) | ...
              (l0 + l1 .* mpc.dcline(:, c.PMAX) < 0) );
    if ~isempty(k)
        warning('toggle_dcline: loss can be negative for DC line from bus %d to %d\n', ...
            [mpc.dcline(k, c.F_BUS:c.T_BUS)]');
    end

    %% add callback functions
    %% note: assumes all necessary data included in 1st arg (mpc, om, results)
    %%       so, no additional explicit args are needed
    mpc = add_userfcn(mpc, 'ext2int', @userfcn_dcline_ext2int);
    mpc = add_userfcn(mpc, 'formulation', @userfcn_dcline_formulation);
    mpc = add_userfcn(mpc, 'int2ext', @userfcn_dcline_int2ext);
    mpc = add_userfcn(mpc, 'printpf', @userfcn_dcline_printpf);
    mpc = add_userfcn(mpc, 'savecase', @userfcn_dcline_savecase);
    mpc.userfcn.status.dcline = 1;
elseif strcmp(upper(on_off), 'OFF')
    mpc = remove_userfcn(mpc, 'savecase', @userfcn_dcline_savecase);
    mpc = remove_userfcn(mpc, 'printpf', @userfcn_dcline_printpf);
    mpc = remove_userfcn(mpc, 'int2ext', @userfcn_dcline_int2ext);
    mpc = remove_userfcn(mpc, 'formulation', @userfcn_dcline_formulation);
    mpc = remove_userfcn(mpc, 'ext2int', @userfcn_dcline_ext2int);
    mpc.userfcn.status.dcline = 0;
elseif strcmp(upper(on_off), 'STATUS')
    if isfield(mpc, 'userfcn') && isfield(mpc.userfcn, 'status') && ...
            isfield(mpc.userfcn.status, 'dcline')
        mpc = mpc.userfcn.status.dcline;
    else
        mpc = 0;
    end
else
    error('toggle_dcline: 2nd argument must be ''on'', ''off'' or ''status''');
end


%%-----  ext2int  ------------------------------------------------------
function mpc = userfcn_dcline_ext2int(mpc, mpopt, args)
%
%   mpc = userfcn_dcline_ext2int(mpc, mpopt, args)
%
%   This is the 'ext2int' stage userfcn callback that prepares the input
%   data for the formulation stage. It expects to find a 'dcline' field
%   in mpc as described above. The optional args are not currently used.
%   It adds two dummy generators for each in-service DC line, with the
%   appropriate upper and lower generation bounds and corresponding
%   entries in gencost. It also expands columns of any A and N matrices
%   accordingly, if present.

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
c = idx_dcline;

%% initialize some things
if isfield(mpc, 'dclinecost')
    havecost = 1;
else
    havecost = 0;
end

%% save version with external indexing
mpc.order.ext.dcline = mpc.dcline;              %% external indexing
if havecost
    mpc.order.ext.dclinecost = mpc.dclinecost;  %% external indexing
end

%% work with only in-service DC lines
mpc.order.dcline.status.on  = find(mpc.dcline(:, c.BR_STATUS) >  0);
mpc.order.dcline.status.off = find(mpc.dcline(:, c.BR_STATUS) <= 0);

%% remove out-of-service DC lines
dc = mpc.dcline(mpc.order.dcline.status.on, :); %% only in-service DC lines
if havecost
    dcc = mpc.dclinecost(mpc.order.dcline.status.on, :);    %% only in-service DC lines
    mpc.dclinecost = dcc;
end
ndc = size(dc, 1);          %% number of in-service DC lines
o = mpc.order;

%%-----  convert stuff to internal indexing  -----
dc(:, c.F_BUS) = o.bus.e2i(dc(:, c.F_BUS));
dc(:, c.T_BUS) = o.bus.e2i(dc(:, c.T_BUS));
mpc.dcline = dc;

%%-----  create gens to represent DC line terminals  -----
%% ensure consistency of initial values of PF, PT and losses
%% (for simple power flow cases)
dc(:, c.PT) = dc(:, c.PF) - (dc(:, c.LOSS0) + dc(:, c.LOSS1) .* dc(:, c.PF));

%% create gens
fg = zeros(ndc, size(mpc.gen, 2));
fg(:, MBASE)        = 100;
fg(:, GEN_STATUS)   =  dc(:, c.BR_STATUS);  %% status (should be all 1's)
fg(:, PMIN)         = -Inf;
fg(:, PMAX)         =  Inf;
tg = fg;
fg(:, GEN_BUS)      =  dc(:, c.F_BUS);      %% from bus
tg(:, GEN_BUS)      =  dc(:, c.T_BUS);      %% to bus
fg(:, PG)           = -dc(:, c.PF);         %% flow (extracted at "from")
tg(:, PG)           =  dc(:, c.PT);         %% flow (injected at "to")
fg(:, QG)           =  dc(:, c.QF);         %% VAr injection at "from"
tg(:, QG)           =  dc(:, c.QT);         %% VAr injection at "to"
fg(:, VG)           =  dc(:, c.VF);         %% voltage set-point at "from"
tg(:, VG)           =  dc(:, c.VT);         %% voltage set-point at "to"
k = find(dc(:, c.PMIN) >= 0);           %% min positive direction flow
if ~isempty(k)                              %% contrain at "from" end
    fg(k, PMAX)     = -dc(k, c.PMIN);       %% "from" extraction lower lim
end
k = find(dc(:, c.PMAX) >= 0);           %% max positive direction flow
if ~isempty(k)                              %% contrain at "from" end
    fg(k, PMIN)     = -dc(k, c.PMAX);       %% "from" extraction upper lim
end
k = find(dc(:, c.PMIN) < 0);            %% max negative direction flow
if ~isempty(k)                              %% contrain at "to" end
    tg(k, PMIN)     =  dc(k, c.PMIN);       %% "to" injection lower lim
end
k = find(dc(:, c.PMAX) < 0);            %% min negative direction flow
if ~isempty(k)                              %% contrain at "to" end
    tg(k, PMAX)     =  dc(k, c.PMAX);       %% "to" injection upper lim
end
fg(:, QMIN)         =  dc(:, c.QMINF);      %% "from" VAr injection lower lim
fg(:, QMAX)         =  dc(:, c.QMAXF);      %% "from" VAr injection upper lim
tg(:, QMIN)         =  dc(:, c.QMINT);      %%  "to"  VAr injection lower lim
tg(:, QMAX)         =  dc(:, c.QMAXT);      %%  "to"  VAr injection upper lim

%% fudge PMAX a bit if necessary to avoid triggering
%% dispatchable load constant power factor constraints
fg(isload(fg), PMAX) = -1e-6;
tg(isload(tg), PMAX) = -1e-6;

%% set all terminal buses to PV (except ref bus)
refbus = find(mpc.bus(:, BUS_TYPE) == REF);
mpc.bus(dc(:, c.F_BUS), BUS_TYPE) = PV;
mpc.bus(dc(:, c.T_BUS), BUS_TYPE) = PV;
mpc.bus(refbus, BUS_TYPE) = REF;

%% expand A and N, if present
nb = size(mpc.bus, 1);
ng = size(mpc.gen, 1);
if isfield(mpc, 'A') && ~isempty(mpc.A)
    [mA, nA] = size(mpc.A);
    if nA >= 2*nb + 2*ng    %% assume AC dimensions
        mpc.A = [   mpc.A(:, 1:2*nb+ng)      sparse(mA, 2*ndc) ...
                    mpc.A(:, 2*nb+ng+(1:ng)) sparse(mA, 2*ndc) ...
                    mpc.A(:, (2*nb+2*ng+1:nA)) ];
    else                    %% assume DC dimensions
        mpc.A = [   mpc.A(:, 1:nb+ng)   sparse(mA, 2*ndc) ...
                    mpc.A(:, (nb+ng+1:nA)) ];
    end
end
if isfield(mpc, 'N') && ~isempty(mpc.N)
    [mN, nN] = size(mpc.N);
    if nN >= 2*nb + 2*ng    %% assume AC dimensions
        mpc.N = [   mpc.N(:, 1:2*nb+ng)      sparse(mN, 2*ndc) ...
                    mpc.N(:, 2*nb+ng+(1:ng)) sparse(mN, 2*ndc) ...
                    mpc.N(:, (2*nb+2*ng+1:nN)) ];
    else                    %% assume DC dimensions
        mpc.N = [   mpc.N(:, 1:nb+ng)   sparse(mN, 2*ndc) ...
                    mpc.N(:, (nb+ng+1:nN)) ];
    end
end

%% append dummy gens
mpc.gen = [mpc.gen; fg; tg];

%% gencost
if isfield(mpc, 'gencost') && ~isempty(mpc.gencost)
    ngcc = size(mpc.gencost, 2);    %% dimensions of gencost
    [pc, qc] = pqcost(mpc.gencost, ng); %% split out re/active costs
    if havecost                 %% user has provided costs
        ndccc = size(dcc, 2);       %% number of dclinecost columns
        ccc = max([ngcc; ndccc]);   %% number of columns in new gencost
        if ccc > ngcc               %% right zero-pad gencost
            pc = [pc zeros(ng, ccc-ngcc)];
            if ~isempty(qc)         %% pad Qg costs, too
                qc = [qc zeros(ng, ccc-ngcc)];
            end
            ngcc = ccc;
        end

        %% flip function across vertical axis and append to gencost
        %% (PF for DC line = -PG for dummy gen at "from" bus)
        for k = 1:ndc
            if dcc(k, MODEL) == POLYNOMIAL
                nc = dcc(k, NCOST);
                temp = dcc(k, NCOST+(1:nc));
                %% flip sign on coefficients of odd terms
                %% (every other starting with linear term,
                %%  that is, the next to last one)
                temp((nc-1):-2:1) = -temp((nc-1):-2:1);
            else  %% dcc(k, MODEL) == PW_LINEAR
                nc = dcc(k, NCOST);
                temp = dcc(k, NCOST+(1:2*nc));
                %% switch sign on horizontal coordinate
                xx = -temp(1:2:2*nc);
                yy =  temp(2:2:2*nc);
                temp(1:2:2*nc) = xx(end:-1:1);
                temp(2:2:2*nc) = yy(end:-1:1);
            end
            padding = zeros(1, ngcc-NCOST-length(temp));
            gck = [dcc(k, 1:NCOST) temp padding];
            
            %% append to gencost
            pc = [pc; gck];
        end
        %% define zero cost for "to" end gen
        dcgc = ones(ndc, 1) * [2 0 0 2 zeros(1, ngcc-4)];
    else
        %% use zero cost as default on "from" end gen
        dcgc = ones(ndc, 1) * [2 0 0 2 zeros(1, ngcc-4)];
        pc = [pc; dcgc];
    end
    %% always use zero cost on "to" end gen
    pc = [pc; dcgc];

    %% reassemble gencost
    if ~isempty(qc)
        qc = [qc; dcgc; dcgc];  %% always use zero cost Qg, both ends
        mpc.gencost = [pc; qc];
    else
        mpc.gencost = pc;
    end
end


%%-----  formulation  --------------------------------------------------
function om = userfcn_dcline_formulation(om, mpopt, args)
%
%   om = userfcn_dcline_formulation(om, mpopt, args)
%
%   This is the 'formulation' stage userfcn callback that defines the
%   user constraints for the dummy generators representing DC lines.
%   It expects to find a 'dcline' field in the mpc stored in om, as
%   described above. By the time it is passed to this callback,
%   MPC.dcline should contain only in-service lines and the from and
%   two bus columns should be converted to internal indexing. The
%   optional args are not currently used.
%
%   If Pf, Pt and Ploss are the flow at the "from" end, flow at the
%   "to" end and loss respectively, and L0 and L1 are the linear loss
%   coefficients, the the relationships between them is given by:
%       Pf - Ploss = Pt
%       Ploss = L0 + L1 * Pf
%   If Pgf and Pgt represent the injections of the dummy generators
%   representing the DC line injections into the network, then
%   Pgf = -Pf and Pgt = Pt, and we can combine all of the above to
%   get the following constraint on Pgf ang Pgt:
%       -Pgf - (L0 - L1 * Pgf) = Pgt
%   which can be written:
%       -L0 <= (1 - L1) * Pgf + Pgt <= -L0

%% define named indices into data matrices
c = idx_dcline;

%% initialize some things
mpc = om.get_mpc();
dc = mpc.dcline;
ndc = size(dc, 1);              %% number of in-service DC lines
ng  = size(mpc.gen, 1) - 2*ndc; %% number of original gens/disp loads

%% constraints
nL0 = -dc(:, c.LOSS0) / mpc.baseMVA;
L1  =  dc(:, c.LOSS1);
Adc = [sparse(ndc, ng) spdiags(1-L1, 0, ndc, ndc) speye(ndc, ndc)];

%% add them to the model
om.add_lin_constraint('dcline', Adc, nL0, nL0, {'Pg'});


%%-----  int2ext  ------------------------------------------------------
function results = userfcn_dcline_int2ext(results, mpopt, args)
%
%   results = userfcn_dcline_int2ext(results, mpopt, args)
%
%   This is the 'int2ext' stage userfcn callback that converts everything
%   back to external indexing and packages up the results. It expects to
%   find a 'dcline' field in the results struct as described for mpc
%   above. It also expects that the last 2*ndc entries in the gen and
%   gencost matrices correspond to the in-service DC lines (where ndc is
%   the number of rows in MPC.dcline. These extra rows are removed from
%   gen and gencost and the flow is taken from the PG of these gens and
%   placed in the flow column of the appropiate dcline row. Corresponding
%   columns are also removed from any A and N matrices, if present. The
%   optional args are not currently used.

%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
c = idx_dcline;

%% initialize some things
o = results.order;
k = find(o.ext.dcline(:, c.BR_STATUS));
ndc = length(k);                    %% number of in-service DC lines
ng  = size(results.gen, 1) - 2*ndc; %% number of original gens/disp loads
nb  = size(results.bus, 1);

%% extract dummy gens
fg = results.gen(ng    +(1:ndc), :);
tg = results.gen(ng+ndc+(1:ndc), :);

%% remove dummy gens
results.gen     = results.gen(1:ng, :);
if isfield(results, 'gencost') && ~isempty(results.gencost)
    results.gencost = results.gencost(1:ng, :);
end

%% delete corresponding columns from A and N, if present
if isfield(results, 'A') && ~isempty(results.A)
    [mA, nA] = size(results.A);
    if nA >= 2*nb + 2*ng + 4*ndc    %% assume AC dimensions
        results.A = results.A(:, [1:2*nb+ng 2*nb+ng+2*ndc+(1:ng) 2*nb+2*ng+4*ndc+1:nA]);
    else                            %% assume DC dimensions
        results.A = results.A(:, [1:nb+ng nb+ng+2*ndc+1:nA]);
    end
end
if isfield(results, 'N') && ~isempty(results.N)
    [mN, nN] = size(results.N);
    if nN >= 2*nb + 2*ng + 4*ndc    %% assume AC dimensions
        results.N = results.N(:, [1:2*nb+ng 2*nb+ng+2*ndc+(1:ng) 2*nb+2*ng+4*ndc+1:nN]);
    else                            %% assume DC dimensions
        results.N = results.N(:, [1:nb+ng nb+ng+2*ndc+1:nN]);
    end
end

%% get the solved flows
results.dcline(:, c.PF) = -fg(:, PG);
results.dcline(:, c.PT) =  tg(:, PG);
results.dcline(:, c.QF) =  fg(:, QG);
results.dcline(:, c.QT) =  tg(:, QG);
results.dcline(:, c.VF) =  fg(:, VG);
results.dcline(:, c.VT) =  tg(:, VG);
if size(fg, 2) >= MU_QMIN
    results.dcline(:, c.MU_PMIN ) = fg(:, MU_PMAX) + tg(:, MU_PMIN);
    results.dcline(:, c.MU_PMAX ) = fg(:, MU_PMIN) + tg(:, MU_PMAX);
    results.dcline(:, c.MU_QMINF) = fg(:, MU_QMIN);
    results.dcline(:, c.MU_QMAXF) = fg(:, MU_QMAX);
    results.dcline(:, c.MU_QMINT) = tg(:, MU_QMIN);
    results.dcline(:, c.MU_QMAXT) = tg(:, MU_QMAX);
end

%%-----  convert stuff back to external indexing  -----
results.order.int.dcline = results.dcline;  %% save internal version
%% copy results to external version
o.ext.dcline(k, c.PF:c.VT) = results.dcline(:, c.PF:c.VT);
if size(results.dcline, 2) == c.MU_QMAXT
    o.ext.dcline(k, c.MU_PMIN:c.MU_QMAXT) = results.dcline(:, c.MU_PMIN:c.MU_QMAXT);
end
results.dcline = o.ext.dcline;              %% use external version


%%-----  printpf  ------------------------------------------------------
function results = userfcn_dcline_printpf(results, fd, mpopt, args)
%
%   results = userfcn_dcline_printpf(results, fd, mpopt, args)
%
%   This is the 'printpf' stage userfcn callback that pretty-prints the
%   results. It expects a results struct, a file descriptor and a MATPOWER
%   options struct. The optional args are not currently used.

%% define named indices into data matrices
c = idx_dcline;

%% options
SUPPRESS        = mpopt.out.suppress_detail;
if SUPPRESS == -1
    if size(results.bus, 1) > 500
        SUPPRESS = 1;
    else
        SUPPRESS = 0;
    end
end
OUT_ALL = mpopt.out.all;
OUT_BRANCH      = OUT_ALL == 1 || (OUT_ALL == -1 && ~SUPPRESS && mpopt.out.branch);
if OUT_ALL == -1
    OUT_ALL_LIM = ~SUPPRESS * mpopt.out.lim.all;
elseif OUT_ALL == 1
    OUT_ALL_LIM = 2;
else
    OUT_ALL_LIM = 0;
end
if OUT_ALL_LIM == -1
    OUT_LINE_LIM    = ~SUPPRESS * mpopt.out.lim.line;
else
    OUT_LINE_LIM    = OUT_ALL_LIM;
end
ctol = mpopt.opf.violation; %% constraint violation tolerance
ptol = 1e-4;                %% tolerance for displaying shadow prices

%%-----  print results  -----
dc = results.dcline;
ndc = size(dc, 1);
kk = find(dc(:, c.BR_STATUS) ~= 0);
if OUT_BRANCH
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n|     DC Line Data                                                             |');
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n Line    From     To        Power Flow           Loss     Reactive Inj (MVAr)');
    fprintf(fd, '\n   #      Bus     Bus   From (MW)   To (MW)      (MW)       From        To   ');
    fprintf(fd, '\n------  ------  ------  ---------  ---------  ---------  ---------  ---------');
    loss = 0;
    for k = 1:ndc
        if dc(k, c.BR_STATUS)   %% status on
            fprintf(fd, '\n%5d%8d%8d%11.2f%11.2f%11.2f%11.2f%11.2f', ...
                        k, dc(k, c.F_BUS:c.T_BUS), dc(k, c.PF:c.PT), ...
                        dc(k, c.PF) - dc(k, c.PT), dc(k, c.QF:c.QT) );
            loss = loss + dc(k, c.PF) - dc(k, c.PT);
        else
            fprintf(fd, '\n%5d%8d%8d%11s%11s%11s%11s%11s', ...
                        k, dc(k, c.F_BUS:c.T_BUS), '-  ', '-  ', '-  ', '-  ', '-  ');
        end
    end
    fprintf(fd, '\n                                              ---------');
    fprintf(fd, '\n                                     Total:%11.2f\n', loss);
end

if OUT_LINE_LIM == 2 || (OUT_LINE_LIM == 1 && ...
        (any(dc(kk, c.PF) > dc(kk, c.PMAX) - ctol) || ...
         any(dc(kk, c.MU_PMIN) > ptol) || ...
         any(dc(kk, c.MU_PMAX) > ptol)))
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n|     DC Line Constraints                                                      |');
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n Line    From     To          Minimum        Actual Flow       Maximum');
    fprintf(fd, '\n   #      Bus     Bus    Pmin mu     Pmin       (MW)       Pmax      Pmax mu ');
    fprintf(fd, '\n------  ------  ------  ---------  ---------  ---------  ---------  ---------');
    for k = 1:ndc
        if OUT_LINE_LIM == 2 || (OUT_LINE_LIM == 1 && ...
                (dc(k, c.PF) > dc(k, c.PMAX) - ctol || ...
                 dc(k, c.MU_PMIN) > ptol || ...
                 dc(k, c.MU_PMAX) > ptol))
            if dc(k, c.BR_STATUS)   %% status on
                fprintf(fd, '\n%5d%8d%8d', k, dc(k, c.F_BUS:c.T_BUS) );
                if dc(k, c.MU_PMIN) > ptol
                    fprintf(fd, '%11.3f', dc(k, c.MU_PMIN) );
                else
                    fprintf(fd, '%11s', '-  ' );
                end
                fprintf(fd, '%11.2f%11.2f%11.2f', ...
                            dc(k, c.PMIN), dc(k, c.PF), dc(k, c.PMAX) );
                if dc(k, c.MU_PMAX) > ptol
                    fprintf(fd, '%11.3f', dc(k, c.MU_PMAX) );
                else
                    fprintf(fd, '%11s', '-  ' );
                end
            else
                fprintf(fd, '\n%5d%8d%8d%11s%11s%11s%11s%11s', ...
                            k, dc(k, c.F_BUS:c.T_BUS), '-  ', '-  ', '-  ', '-  ', '-  ');
            end
        end
    end
    fprintf(fd, '\n');
end


%%-----  savecase  -----------------------------------------------------
function mpc = userfcn_dcline_savecase(mpc, fd, prefix, args)
%
%   mpc = userfcn_dcline_savecase(mpc, fd, prefix, args)
%
%   This is the 'savecase' stage userfcn callback that prints the M-file
%   code to save the 'dcline' field in the case file. It expects a
%   MATPOWER case struct (mpc), a file descriptor and variable prefix
%   (usually 'mpc.'). The optional args are not currently used.

%% define named indices into data matrices
c = idx_dcline;

%% save it
ncols = size(mpc.dcline, 2);
fprintf(fd, '\n%%%%-----  DC Line Data  -----%%%%\n');
if ncols < c.MU_QMAXT
    fprintf(fd, '%%\tfbus\ttbus\tstatus\tPf\tPt\tQf\tQt\tVf\tVt\tPmin\tPmax\tQminF\tQmaxF\tQminT\tQmaxT\tloss0\tloss1\n');
else
    fprintf(fd, '%%\tfbus\ttbus\tstatus\tPf\tPt\tQf\tQt\tVf\tVt\tPmin\tPmax\tQminF\tQmaxF\tQminT\tQmaxT\tloss0\tloss1\tmuPmin\tmuPmax\tmuQminF\tmuQmaxF\tmuQminT\tmuQmaxT\n');
end
template = '\t%d\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g';
if ncols == c.MU_QMAXT
    template = [template, '\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f'];
end
template = [template, ';\n'];
fprintf(fd, '%sdcline = [\n', prefix);
fprintf(fd, template, mpc.dcline.');
fprintf(fd, '];\n');

%% to do, make it save mpc.dclinecost too (somehow I forgot this)
%% have a look at the saving of gencost in savecase for template
