function mpc = toggle_iflims(mpc, on_off)
%TOGGLE_IFLIMS Enable, disable or check status of set of interface flow limits.
%   MPC = TOGGLE_IFLIMS(MPC, 'on')
%   MPC = TOGGLE_IFLIMS(MPC, 'off')
%   T_F = TOGGLE_IFLIMS(MPC, 'status')
%
%   Enables, disables or checks the status of a set of OPF userfcn
%   callbacks to implement interface flow limits based on a DC flow model.
%
%   These callbacks expect to find an 'if' field in the input MPC, where
%   MPC.if is a struct with the following fields:
%       map     n x 2, defines each interface in terms of a set of 
%               branch indices and directions. Interface I is defined
%               by the set of rows whose 1st col is equal to I. The
%               2nd column is a branch index multiplied by 1 or -1
%               respectively for lines whose orientation is the same
%               as or opposite to that of the interface.
%       lims    nif x 3, defines the DC model flow limits in MW
%               for specified interfaces. The first column is the index
%               of the interface, the 2nd and 3rd columns specify
%               the lower and upper limits on the (DC model) flow
%               across the interface, respectively. Normally, the lower
%               limit is negative, indicating a flow in the opposite
%               direction.
%
%   The 'int2ext' callback also packages up results and stores them in
%   the following output fields of results.if:
%       P       - nif x 1, actual flow across each interface in MW
%       mu.l    - nif x 1, shadow price on lower flow limit, ($/MW)
%       mu.u    - nif x 1, shadow price on upper flow limit, ($/MW)
%
%   See also ADD_USERFCN, REMOVE_USERFCN, RUN_USERFCN, T_CASE30_USERFCNS.

%   MATPOWER
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if strcmp(upper(on_off), 'ON')
    %% check for proper reserve inputs
    if ~isfield(mpc, 'if') || ~isstruct(mpc.if) || ...
            ~isfield(mpc.if, 'map') || ...
            ~isfield(mpc.if, 'lims')
        error('toggle_iflims: case must contain an ''if'' field, a struct defining ''map'' and ''lims''');
    end
    
    %% add callback functions
    %% note: assumes all necessary data included in 1st arg (mpc, om, results)
    %%       so, no additional explicit args are needed
    mpc = add_userfcn(mpc, 'ext2int', @userfcn_iflims_ext2int);
    mpc = add_userfcn(mpc, 'formulation', @userfcn_iflims_formulation);
    mpc = add_userfcn(mpc, 'int2ext', @userfcn_iflims_int2ext);
    mpc = add_userfcn(mpc, 'printpf', @userfcn_iflims_printpf);
    mpc = add_userfcn(mpc, 'savecase', @userfcn_iflims_savecase);
    mpc.userfcn.status.iflims = 1;
elseif strcmp(upper(on_off), 'OFF')
    mpc = remove_userfcn(mpc, 'savecase', @userfcn_iflims_savecase);
    mpc = remove_userfcn(mpc, 'printpf', @userfcn_iflims_printpf);
    mpc = remove_userfcn(mpc, 'int2ext', @userfcn_iflims_int2ext);
    mpc = remove_userfcn(mpc, 'formulation', @userfcn_iflims_formulation);
    mpc = remove_userfcn(mpc, 'ext2int', @userfcn_iflims_ext2int);
    mpc.userfcn.status.iflims = 0;
elseif strcmp(upper(on_off), 'STATUS')
    if isfield(mpc, 'userfcn') && isfield(mpc.userfcn, 'status') && ...
            isfield(mpc.userfcn.status, 'iflims')
        mpc = mpc.userfcn.status.iflims;
    else
        mpc = 0;
    end
else
    error('toggle_iflims: 2nd argument must be ''on'', ''off'' or ''status''');
end


%%-----  ext2int  ------------------------------------------------------
function mpc = userfcn_iflims_ext2int(mpc, mpopt, args)
%
%   mpc = userfcn_iflims_ext2int(mpc, mpopt, args)
%
%   This is the 'ext2int' stage userfcn callback that prepares the input
%   data for the formulation stage. It expects to find an 'if' field in
%   mpc as described above. The optional args are not currently used.

%% initialize some things
ifmap = mpc.if.map;
o = mpc.order;
nl0 = size(o.ext.branch, 1);    %% original number of branches
nl = size(mpc.branch, 1);       %% number of on-line branches

%% save if.map for external indexing
mpc.order.ext.ifmap = ifmap;

%%-----  convert stuff to internal indexing  -----
e2i = zeros(nl0, 1);
e2i(o.branch.status.on) = (1:nl)';  %% ext->int branch index mapping
d = sign(ifmap(:, 2));
br = abs(ifmap(:, 2));
ifmap(:, 2) = d .* e2i(br);
ifmap(ifmap(:, 2) == 0, :) = [];    %% delete branches that are out

mpc.if.map = ifmap;


%%-----  formulation  --------------------------------------------------
function om = userfcn_iflims_formulation(om, mpopt, args)
%
%   om = userfcn_iflims_formulation(om, mpopt, args)
%
%   This is the 'formulation' stage userfcn callback that defines the
%   user costs and constraints for interface flow limits. It expects to
%   find an 'if' field in the mpc stored in om, as described above. The
%   optional args are not currently used.

%% define named indices into data matrices
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% initialize some things
mpc = om.get_mpc();
[baseMVA, bus, branch] = deal(mpc.baseMVA, mpc.bus, mpc.branch);
ifmap = mpc.if.map;
iflims = mpc.if.lims;

%% form B matrices for DC model
[B, Bf, Pbusinj, Pfinj] = makeBdc(baseMVA, bus, branch);
n = size(Bf, 2);                    %% dim of theta

%% form constraints
ifidx = unique(iflims(:, 1));   %% interface number list
nifs = length(ifidx);           %% number of interfaces
Aif = sparse(nifs, n);
lif = zeros(nifs, 1);
uif = zeros(nifs, 1);
for k = 1:nifs
    %% extract branch indices
    br = ifmap(ifmap(:, 1) == ifidx(k), 2);
    if isempty(br)
        error('userfcn_iflims_formulation: interface %d has no in-service branches', k);
    end
    d = sign(br);
    br = abs(br);
    Ak = sparse(1, n);              %% Ak = sum( d(i) * Bf(i, :) )
    bk = 0;                         %% bk = sum( d(i) * Pfinj(i) )
    for i = 1:length(br)
        Ak = Ak + d(i) * Bf(br(i), :);
        bk = bk + d(i) * Pfinj(br(i));
    end
    Aif(k, :) = Ak;
    lif(k) = iflims(k, 2) / baseMVA - bk;
    uif(k) = iflims(k, 3) / baseMVA - bk;
end

%% add interface constraint
om.add_lin_constraint('iflims',  Aif, lif, uif, {'Va'});   %% nifs


%%-----  int2ext  ------------------------------------------------------
function results = userfcn_iflims_int2ext(results, mpopt, args)
%
%   results = userfcn_iflims_int2ext(results, mpopt, args)
%
%   This is the 'int2ext' stage userfcn callback that converts everything
%   back to external indexing and packages up the results. It expects to
%   find an 'if' field in the results struct as described for mpc above.
%   It also expects the results to contain solved branch flows and linear
%   constraints named 'iflims' which are used to populate output fields
%   in results.if. The optional args are not currently used.

%% define named indices into data matrices
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% get internal ifmap
ifmap = results.if.map;
iflims = results.if.lims;

%%-----  convert stuff back to external indexing  -----
results.if.map = results.order.ext.ifmap;

%%-----  results post-processing  -----
ifidx = unique(iflims(:, 1));   %% interface number list
nifs = length(ifidx);           %% number of interfaces
results.if.P = zeros(nifs, 1);
for k = 1:nifs
    %% extract branch indices
    br = ifmap(ifmap(:, 1) ==  ifidx(k), 2);
    d = sign(br);
    br = abs(br);
    results.if.P(k) = sum( d .* results.branch(br, PF) );
end
results.if.mu.l = results.lin.mu.l.iflims / results.baseMVA;
results.if.mu.u = results.lin.mu.u.iflims / results.baseMVA;


%%-----  printpf  ------------------------------------------------------
function results = userfcn_iflims_printpf(results, fd, mpopt, args)
%
%   results = userfcn_iflims_printpf(results, fd, mpopt, args)
%
%   This is the 'printpf' stage userfcn callback that pretty-prints the
%   results. It expects a results struct, a file descriptor and a MATPOWER
%   options struct. The optional args are not currently used.

%% define named indices into data matrices
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%%-----  print results  -----
ptol = 1e-6;        %% tolerance for displaying shadow prices
isOPF           = isfield(results, 'f') && ~isempty(results.f);
SUPPRESS        = mpopt.out.suppress_detail;
if SUPPRESS == -1
    if size(results.bus, 1) > 500
        SUPPRESS = 1;
    else
        SUPPRESS = 0;
    end
end
OUT_ALL         = mpopt.out.all;
OUT_FORCE       = mpopt.out.force;
OUT_BRANCH      = OUT_ALL == 1 || (OUT_ALL == -1 && ~SUPPRESS && mpopt.out.branch);

if isOPF && OUT_BRANCH && (results.success || OUT_FORCE)
    iflims = results.if.lims;
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n|     Interface Flow Limits                                                    |');
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n Interface  Shadow Prc  Lower Lim      Flow      Upper Lim   Shadow Prc');
    fprintf(fd, '\n     #        ($/MW)       (MW)        (MW)        (MW)       ($/MW)   ');
    fprintf(fd, '\n----------  ----------  ----------  ----------  ----------  -----------');
    ifidx = unique(iflims(:, 1));   %% interface number list
    nifs = length(ifidx);           %% number of interfaces
    for k = 1:nifs
        fprintf(fd, '\n%6d ', iflims(k, 1));
        if results.if.mu.l(k) > ptol
            fprintf(fd, '%14.3f', results.if.mu.l(k));
        else
            fprintf(fd, '          -   ');
        end
        fprintf(fd, '%12.2f%12.2f%12.2f', iflims(k, 2), results.if.P(k), iflims(k, 3));
        if results.if.mu.u(k) > ptol
            fprintf(fd, '%13.3f', results.if.mu.u(k));
        else
            fprintf(fd, '         -     ');
        end
    end
    fprintf(fd, '\n');
end


%%-----  savecase  -----------------------------------------------------
function mpc = userfcn_iflims_savecase(mpc, fd, prefix, args)
%
%   mpc = userfcn_iflims_savecase(mpc, fd, prefix, args)
%
%   This is the 'savecase' stage userfcn callback that prints the M-file
%   code to save the 'if' field in the case file. It expects a
%   MATPOWER case struct (mpc), a file descriptor and variable prefix
%   (usually 'mpc.'). The optional args are not currently used.

ifmap = mpc.if.map;
iflims = mpc.if.lims;

fprintf(fd, '\n%%%%-----  Interface Flow Limit Data  -----%%%%\n');
fprintf(fd, '%%%% interface<->branch map data\n');
fprintf(fd, '%%\tifnum\tbranchidx (negative defines opposite direction)\n');
fprintf(fd, '%sif.map = [\n', prefix);
fprintf(fd, '\t%d\t%d;\n', ifmap');
fprintf(fd, '];\n');

fprintf(fd, '\n%%%% interface flow limit data (based on DC model)\n');
fprintf(fd, '%%%% (lower limit should be negative for opposite direction)\n');
fprintf(fd, '%%\tifnum\tlower\tupper\n');
fprintf(fd, '%sif.lims = [\n', prefix);
fprintf(fd, '\t%d\t%g\t%g;\n', iflims');
fprintf(fd, '];\n');

%% save output fields for solved case
if isfield(mpc.if, 'P')
    if exist('serialize', 'file') == 2
        fprintf(fd, '\n%%%% solved values\n');
        fprintf(fd, '%sif.P = %s\n', prefix, serialize(mpc.if.P));
        fprintf(fd, '%sif.mu.l = %s\n', prefix, serialize(mpc.if.mu.l));
        fprintf(fd, '%sif.mu.u = %s\n', prefix, serialize(mpc.if.mu.u));
    else
        url = 'http://www.mathworks.com/matlabcentral/fileexchange/12063';
        warning('MATPOWER:serialize', ...
            'userfcn_iflims_savecase: Cannot save the ''iflims'' output fields without the ''serialize'' function, which is available as a free download from:\n<%s>\n\n', url);
    end
end
