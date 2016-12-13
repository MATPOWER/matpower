function mpc = toggle_softlims(mpc, on_off)
%TOGGLE_SOFTLIMS Relax DC optimal power flow branch limits.
%   MPC = TOGGLE_SOFTLIMS(MPC, 'on')
%   MPC = TOGGLE_SOFTLIMS(MPC, 'off')
%   T_F = TOGGLE_SOFTLIMS(MPC, 'status')
%
%   Enables, disables or checks the status of a set of OPF userfcn
%   callbacks to implement relaxed branch flow limits a DC optimal flow model.
%
%   These callbacks expect to find a 'softlims' field in the input MPC,
%   where MPC.softlims is a struct with the following fields:
%       idx     (optional) n x 1, index vector for branches whos flow
%               limits are to be relaxed, default is to use all on-line
%               branches with non-zero limits specified in RATE_A
%       cost    (optional) n x 1, linear marginal cost per MW of exceeding
%               RATE_A. Can optionally be a scalar, in which case it is
%               applied to all soft limits. Default if not specified is
%               $1000/MW.
%
%   The 'int2ext' callback also packages up results and stores them in
%   the following output fields of results.softlims:
%       overload - nl x 1, amount of overload of each line in MW
%       ovl_cost - nl x 1, total cost of overload in $/hr
%
%   The shadow prices on the soft limit constraints are also returned in the
%   MU_SF and MU_ST columns of the branch matrix.
%   Note: These shadow prices are equal to the corresponding hard limit
%       shadow prices when the soft limits are not violated. When violated,
%       the shadow price on a soft limit constraint is equal to the
%       user-specified soft limit violation cost.
%
%   See also ADD_USERFCN, REMOVE_USERFCN, RUN_USERFCN, T_CASE30_USERFCNS.

%   To do for future versions:
%       - make softlims input field optional, convert all lines if missing
%       Inputs:
%       cost    n x 3, linear marginal cost per MW of exceeding each of
%               RATE_A, RATE_B and RATE_C. Columns 2 and 3 are optional.
%       brkpts  n x npts, allow to specify arbitrary breakpoints at which
%               cost increases, defined as percentages above RATE_A.
%       base_flow   n x 1, arbitrary baseline (other than RATE_A)

%   MATPOWER
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if strcmp(upper(on_off), 'ON')
    %% check for proper softlims inputs
    %% (inputs are optional, defaults handled in ext2int callback)
    
    %% add callback functions
    %% note: assumes all necessary data included in 1st arg (mpc, om, results)
    %%       so, no additional explicit args are needed
    mpc = add_userfcn(mpc, 'ext2int', @userfcn_softlims_ext2int);
    mpc = add_userfcn(mpc, 'formulation', @userfcn_softlims_formulation);
    mpc = add_userfcn(mpc, 'int2ext', @userfcn_softlims_int2ext);
    mpc = add_userfcn(mpc, 'printpf', @userfcn_softlims_printpf);
    mpc = add_userfcn(mpc, 'savecase', @userfcn_softlims_savecase);
    mpc.userfcn.status.softlims = 1;
elseif strcmp(upper(on_off), 'OFF')
    mpc = remove_userfcn(mpc, 'savecase', @userfcn_softlims_savecase);
    mpc = remove_userfcn(mpc, 'printpf', @userfcn_softlims_printpf);
    mpc = remove_userfcn(mpc, 'int2ext', @userfcn_softlims_int2ext);
    mpc = remove_userfcn(mpc, 'formulation', @userfcn_softlims_formulation);
    mpc = remove_userfcn(mpc, 'ext2int', @userfcn_softlims_ext2int);
    mpc.userfcn.status.softlims = 0;
elseif strcmp(upper(on_off), 'STATUS')
    if isfield(mpc, 'userfcn') && isfield(mpc.userfcn, 'status') && ...
            isfield(mpc.userfcn.status, 'softlims')
        mpc = mpc.userfcn.status.softlims;
    else
        mpc = 0;
    end
else
    error('toggle_softlims: 2nd argument must be ''on'', ''off'' or ''status''');
end


%%-----  ext2int  ------------------------------------------------------
function mpc = userfcn_softlims_ext2int(mpc, args)
%
%   mpc = userfcn_softlims_ext2int(mpc, args)
%
%   This is the 'ext2int' stage userfcn callback that prepares the input
%   data for the formulation stage. It expects to find a 'softlims' field in
%   mpc as described above. The optional args are not currently used.

%% define named indices into data matrices
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% check for proper softlims inputs
default_cost = 1000;    %% used if cost is not specified
if isfield(mpc, 'softlims')
    if ~isfield(mpc.softlims, 'idx')
        mpc.softlims.idx = [];
    end
    if ~isfield(mpc.softlims, 'cost')
        mpc.softlims.cost = default_cost;
    end
else
    mpc.softlims.idx = [];
    mpc.softlims.cost = default_cost;
end

%% initialize some things
s = mpc.softlims;
o = mpc.order;
nl0 = size(o.ext.branch, 1);    %% original number of branches
nl = size(mpc.branch, 1);       %% number of on-line branches

%% save softlims struct for external indexing
mpc.order.ext.softlims = s;

%%-----  convert stuff to internal indexing  -----
s = softlims_defaults(s, o.ext.branch);     %% get defaults
e2i = zeros(nl0, 1);
e2i(o.branch.status.on) = (1:nl)';  %% ext->int branch index mapping
s.idx = e2i(s.idx);
k = find(s.idx == 0);   %% find idxes corresponding to off-line branches
s.idx(k) = [];          %% delete them
s.cost(k, :) = [];
k = find(mpc.branch(s.idx, RATE_A) <= 0);   %% find branches w/o limits
s.idx(k) = [];          %% delete them
s.cost(k, :) = [];

%%-----  remove hard limits on branches with soft limits  -----
s.Pfmax = mpc.branch(s.idx, RATE_A) / mpc.baseMVA;  %% save limit first
mpc.branch(s.idx, RATE_A) = 0;                      %% then remove it

mpc.softlims = s;
mpc.order.int.softlims = s;


%%-----  formulation  --------------------------------------------------
function om = userfcn_softlims_formulation(om, args)
%
%   om = userfcn_softlims_formulation(om, args)
%
%   This is the 'formulation' stage userfcn callback that defines the
%   user costs and constraints for interface flow limits. It expects to
%   find a 'softlims' field in the mpc stored in om, as described above. The
%   optional args are not currently used.

%% define named indices into data matrices
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% initialize some things
mpc = get_mpc(om);
[baseMVA, bus, branch] = deal(mpc.baseMVA, mpc.bus, mpc.branch);
s = mpc.softlims;

%% form B matrices for DC model
[B, Bf, Pbusinj, Pfinj] = makeBdc(baseMVA, bus, branch);
n = size(Bf, 2);        %% dim of theta

%% form constraints (flv is flow limit violation variable)
%%   -Bf * Va - flv <=  Pfinj + Pfmax
%%    Bf * Va - flv <= -Pfinj + Pfmax
ns = length(s.idx);     %% number of soft limits
I = speye(ns, ns);
As = [-Bf(s.idx, :) -I; Bf(s.idx, :) -I];
ls = -Inf(2*ns, 1);
us = [   Pfinj(s.idx) + s.Pfmax;
        -Pfinj(s.idx) + s.Pfmax ];

%% costs on flv variable
Cw = s.cost(:, 1) * mpc.baseMVA;

%% add vars, costs, constraints
om = add_vars(om, 'flv', ns, zeros(ns, 1), zeros(ns, 1), Inf(ns, 1));
om = add_costs(om, 'vc', struct('N', I, 'Cw', Cw), {'flv'});
om = add_constraints(om, 'softlims',  As, ls, us, {'Va', 'flv'});   %% 2*ns


%%-----  int2ext  ------------------------------------------------------
function results = userfcn_softlims_int2ext(results, args)
%
%   results = userfcn_softlims_int2ext(results, args)
%
%   This is the 'int2ext' stage userfcn callback that converts everything
%   back to external indexing and packages up the results. It expects to
%   find a 'softlims' field in the results struct as described for mpc above.
%   It also expects the results to contain solved branch flows and linear
%   constraints named 'softlims' which are used to populate output fields
%   in results.softlims. The optional args are not currently used.

%% define named indices into data matrices
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% get internal softlims struct
s = results.softlims;

%%-----  convert stuff back to external indexing  -----
o = results.order;
nl0 = size(o.ext.branch, 1);    %% original number of branches
nl = size(results.branch, 1);   %% number of on-line branches
o.branch.status.on;
results.softlims = results.order.ext.softlims;

%%-----  restore hard limits  -----
results.branch(s.idx, RATE_A) = s.Pfmax * results.baseMVA;

%%-----  results post-processing  -----
%% get shadow prices
n = size(results.lin.mu.u.softlims, 1) / 2;
results.branch(s.idx, MU_ST) = results.lin.mu.u.softlims(1:n) / results.baseMVA;
results.branch(s.idx, MU_SF) = results.lin.mu.u.softlims(n+1:end) / results.baseMVA;
%% get overloads and overload costs
results.softlims.overload = zeros(nl0, 1);
k = find(results.branch(:, RATE_A) & ...
         abs(results.branch(:, PF)) > results.branch(:, RATE_A) );
results.softlims.overload(o.branch.status.on(k)) = ...
        abs(results.branch(k, PF)) - results.branch(k, RATE_A);
results.softlims.ovl_cost = zeros(size(results.softlims.overload));
results.softlims.ovl_cost(o.branch.status.on(s.idx)) = ...
    results.softlims.overload(o.branch.status.on(s.idx)) .* s.cost(:, 1);


%%-----  printpf  ------------------------------------------------------
function results = userfcn_softlims_printpf(results, fd, mpopt, args)
%
%   results = userfcn_softlims_printpf(results, fd, mpopt, args)
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
    s = softlims_defaults(results.softlims, results.branch);   %% get defaults
    k = find(s.overload(s.idx) | sum(results.branch(s.idx, MU_SF:MU_ST), 2) > ptol);

    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n|     Soft Flow Limits                                                         |');
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\nBrnch   From   To     Flow      Limit    Overload     mu');
    fprintf(fd, '\n  #     Bus    Bus    (MW)      (MW)       (MW)     ($/MW)');
    fprintf(fd, '\n-----  -----  -----  --------  --------  --------  ---------');
    fprintf(fd, '\n%4d%7d%7d%10.2f%10.2f%10.2f%11.3f', ...
            [   s.idx(k), results.branch(s.idx(k), [F_BUS, T_BUS]), ...
                results.branch(s.idx(k), [PF, RATE_A]), ...
                s.overload(s.idx(k)), ...
                sum(results.branch(s.idx(k), [MU_SF:MU_ST]), 2) ...
            ]');
    fprintf(fd, '\n                                         --------');
    fprintf(fd, '\n                                Total:%10.2f', ...
            sum(s.overload(s.idx(k))));
    fprintf(fd, '\n');
end


%%-----  savecase  -----------------------------------------------------
function mpc = userfcn_softlims_savecase(mpc, fd, prefix, args)
%
%   mpc = userfcn_softlims_savecase(mpc, fd, mpopt, args)
%
%   This is the 'savecase' stage userfcn callback that prints the M-file
%   code to save the 'softlims' field in the case file. It expects a
%   MATPOWER case struct (mpc), a file descriptor and variable prefix
%   (usually 'mpc.'). The optional args are not currently used.

if isfield(mpc, 'softlims')
    s = mpc.softlims;

    fprintf(fd, '\n%%%%-----  Soft Flow Limit Data  -----%%%%\n');

    if isfield(s, 'idx')
        fprintf(fd, '%%%% branch indexes\n');
        fprintf(fd, '%%\tbranchidx\n');
        if isempty(s.idx)
            fprintf(fd, '%ssoftlims.idx = [];\n\n', prefix);
        else
            fprintf(fd, '%ssoftlims.idx = [\n', prefix);
            fprintf(fd, '\t%d;\n', s.idx);
            fprintf(fd, '];\n\n');
        end
    end

    fprintf(fd, '%%%% violation cost coefficients\n');
    fprintf(fd, '%%\trate_a_cost\n');
    fprintf(fd, '%ssoftlims.cost = [\n', prefix);
    fprintf(fd, '\t%g;\n', s.cost);
    fprintf(fd, '];\n');

    %% save output fields for solved case
    if isfield(mpc.softlims, 'overload')
        fprintf(fd, '\n%%%% overloads\n');
        fprintf(fd, '%%\toverload\n');
        fprintf(fd, '%ssoftlims.overload = [\n', prefix);
        fprintf(fd, '\t%g;\n', s.overload);
        fprintf(fd, '];\n');

        fprintf(fd, '\n%%%% overload costs\n');
        fprintf(fd, '%%\toverload_costs\n');
        fprintf(fd, '%ssoftlims.ovl_cost = [\n', prefix);
        fprintf(fd, '\t%g;\n', s.ovl_cost);
        fprintf(fd, '];\n');
    end
end

%%-----  softlims_defaults  --------------------------------------------
function s = softlims_defaults(s, branch)
%
%   s = softlims_defaults(s, branch)
%
%   Takes a softlims struct that could have an empty 'idx' field or a
%   scalar 'cost' field and fills them out with the defaults, where the
%   default for 'idx' includes all on-line branches with non-zero RATE_A,
%   and the default for the cost is to apply the scalar to each soft limit
%   violation.

%% define named indices into data matrices
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

if isempty(s.idx)
    s.idx = find(branch(:, BR_STATUS) > 0 & branch(:, RATE_A) > 0);
end
if length(s.cost) == 1 && length(s.idx) > 1
    s.cost = s.cost * ones(size(s.idx));
end
