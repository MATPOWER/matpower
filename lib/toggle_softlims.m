function mpc = toggle_softlims(mpc, on_off)
%TOGGLE_SOFTLIMS Relax DC optimal power flow branch limits.
%   MPC = TOGGLE_SOFTLIMS(MPC, 'on')
%   MPC = TOGGLE_SOFTLIMS(MPC, 'off')
%   T_F = TOGGLE_SOFTLIMS(MPC, 'status')
%
%   Enables, disables or checks the status of a set of OPF userfcn
%   callbacks to implement relaxed inequality constraints for an OPF model.
%
%   These callbacks expect to find a 'softlims' field in the input MPC,
%   where MPC.softlims is a struct with fields corresponding to the
%   possible limits, namely:
%       VMIN, VMAX, RATE_A, PMIN, PMAX, QMIN, QMAX, ANGMAX, ANGMIN,
%   Each of these is itself a struct with the following fields, all of which
%   are optional:
%       idx     index of affected buses, branches, or generators
%               For buses these are bus numbers. For all others these are
%               indexes into the respective matrix. The default is to include
%               all online elements for which the constraint in question is
%               not unbounded, except for generators, which also exclude those
%               used to model dispatchable loads (i.e. isload(gen) is true).
%       cost    linear marginal cost of exceeding the original limit
%               Defaults are:
%                   $100,000 $/pu   for VMAX and VMIN
%                      $1000 $/MW   for RATE_A, PMAX, and PMIN
%                      $1000 $/MVAr for QMAX, QMIN
%                      $1000 $/deg  for ANGMAX, ANGMIN
%               The final cost is determined by taking the maximum between
%               the cost above and 2*gen_max_cost, where gen_max_cost is
%               the maximum generator marginal cost of online generators.
%
%       hl_mod  type of modification to hard limit, hl, are:
%                   'remove'  : unbounded limit
%                   'replace' : new hard limit specified in hl_val
%                   'scale'   : new hard limit is hl_val*hl
%                   'shift'   : new hard limit is hl + sign*hl where sign
%                               either 1 or -1 depending on the limit
%                   'none'    : No soft limit for this property
%       hl_val  Value used in conjuction with hl_mod to modify the
%               constraint hard limit. It can be specified as either a
%               scalar, or a vector that MUST be the same length as idx.
%               when...
%                    hl_mod = 'remove': hl_val is ALWAYS Inf
%                    hl_mod = 'replace': hl_val MUST be specified.
%                    hl_mod = 'scale': if hl_val is NOT specified positive
%                                      limits (as well as QMIN and ANGMIN)
%                                      are doubled (hl_val = 2) and negative
%                                      limits are halved (hl_val = 0.5).
%                    hl_mod = 'shift': If hl_val is NOT specified voltage
%                                      limits are shifted by 0.25 pu (in
%                                      appropriate direction) and the rest
%                                      are shifted by 10.
%
%       ub      Internally calculated upperbound on the slack variable that
%               implements the sof limit (handled in the defaults function).
%
%       sav     original limits (handled in the defaults function).
%
%       rval    value to place in mpc structure to eliminate the original 
%               constraint (handled in the default function).
%
%   Defaults are assigned depending on the mpopt.opf.softlims.default
%   option. If it is 0, any limit not specified in the softlims structrue
%   is set to hl_mod = 'none' (i.e ignored). If
%   mpopt.opf.softlims.default = 1 (the default) all limits not included in
%   the softlims structure are set to hl_mod = 'remove' except:
%   softlims
%         .VMIN
%           .hl_mod = 'replace'
%           .hl_val = 0
%         .PMIN
%           .hl_mod = 'replace'
%           .hl_val = 0    for normal generators (PMIN > 0)
%           .hl_val = -Inf for for generators with PMIN < 0 AND PMAX > 0
%
%   The 'int2ext' callback also packages up results and stores them in
%   the following output fields of results.softlims.(prop), where prop is
%   one of the above mentioned limits:
%       overload - amount of overload, i.e. violation of hard-limit.
%       ovl_cost - total cost of overload in $/hr
%
%   The shadow prices on the soft limit constraints are also returned in the
%   relevant columns of the respective matrices (MU_SF, MU_ST for RATE_A,
%   MU_VMAX for VMAX etc.)
%   Note: These shadow prices are equal to the corresponding hard limit
%       shadow prices when the soft limits are not violated. When violated,
%       the shadow price on a soft limit constraint is equal to the
%       user-specified soft limit violation cost.
%
%   If mpopt.opf.softlims.default = 1 (default) the default softlimits are
%   applied to unspecified limits in the softlims structure. If
%   mpopt.opf.softlims.default = 0 the unspecified softlims are ignored.
%   In this case, to initialize the default soft limit a blank structure
%   should be created. For example, to use the default settings for VMAX
%   the command, mpc.softlims.VMAX = struct(), should be issued.
%
%   See also ADD_USERFCN, REMOVE_USERFCN, RUN_USERFCN, T_OPF_SOFTLIMS.

%   To do for future versions:
%       Inputs:
%       cost    n x 3, linear marginal cost per MW of exceeding each of
%               RATE_A, RATE_B and RATE_C. Columns 2 and 3 are optional.
%       brkpts  n x npts, allow to specify arbitrary breakpoints at which
%               cost increases, defined as percentages above RATE_A.
%       base_flow   n x 1, arbitrary baseline (other than RATE_A)

%   MATPOWER
%   Copyright (c) 2009-2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Eran Schweitzer, Arizona State University
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
function mpc = userfcn_softlims_ext2int(mpc, mpopt, args)
%
%   mpc = userfcn_softlims_ext2int(mpc, mpopt, args)
%
%   This is the 'ext2int' stage userfcn callback that prepares the input
%   data for the formulation stage. It expects to find a 'softlims' field in
%   mpc as described above. The optional args are not currently used.

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

% structures used to index into softlims in a loop
lims = softlims_lim2mat();
mat2lims = struct(...
    'bus', {{'VMAX', 'VMIN'}}, ...
    'branch', {{'ANGMAX', 'ANGMIN', 'RATE_A'}}, ...
    'gen', {{'PMAX', 'PMIN', 'QMAX', 'QMIN'}} ...
);

%% check opf.softlims.default
if isempty(mpopt)
    warning('userfcn_softlims_ext2int: Assuming ''mpopt.opf.softlims.default'' = 1, since mpopt was not provided.');
    use_default = 1;
else
    use_default = mpopt.opf.softlims.default;
end 

%% set up softlims defaults
if ~isfield(mpc, 'softlims')
    mpc.softlims = struct();
end
mpc.softlims = softlims_defaults(mpc, mpopt);
mpc.softlims = softlims_init(mpc, mpopt);

%% initialize some things
s = mpc.softlims;
o = mpc.order;
%nl0 = size(o.ext.branch, 1);     %% original number of branches
%nl  = size(mpc.branch, 1);       %% number of on-line branches

%% save softlims struct with external indexing
mpc.order.ext.softlims = s;

%%-----  convert stuff to internal indexing  -----
for mat = {'bus', 'branch', 'gen'}
    mat = mat{1};
    n0  = size(o.ext.(mat), 1);  %% original number
    n   = size(mpc.(mat), 1);    %% on-line number
    e2i = zeros(n0, 1);
    if strcmp(mat,'gen')
        % for generators, the permutation of the generator matrix needs to
        % be accounted for
        e2i(o.gen.status.on(o.gen.i2e)) = (1:n)';
    else
        e2i(o.(mat).status.on) = (1:n)';  %% ext->int index mapping
    end
    for prop = mat2lims.(mat)
        if ~strcmp( s.(prop{:}).hl_mod, 'none')
            s.(prop{:}).idx = e2i(s.(prop{:}).idx);
            k = find(s.(prop{:}).idx == 0); %% find idxes corresponding to off-line elements
            s.(prop{:}).idx(k)     = [];    %% delete them
            s.(prop{:}).cost(k, :) = [];
            s.(prop{:}).sav(k)     = [];
            if ~isscalar(s.(prop{:}).ub)
                s.(prop{:}).ub(k)  = [];
            end
            if ~isscalar(s.(prop{:}).hl_val)
                s.(prop{:}).hl_val(k)  = [];
            end
            if isempty(s.(prop{:}).idx)
                s.(prop{:}).hl_mod = 'none';
            end
        end
    end
end

%%%%-------- remove hard limits on elements with soft limits
for prop = fieldnames(lims).'
    if ~strcmp(s.(prop{:}).hl_mod, 'none')
        mat = lims.(prop{:});  %% mpc matrix
        mpc.(mat)(s.(prop{:}).idx, eval(prop{:})) = s.(prop{:}).rval;
    end
end

mpc.softlims = s;
mpc.order.int.softlims = s;


%%-----  formulation  --------------------------------------------------
function om = userfcn_softlims_formulation(om, mpopt, args)
%
%   om = userfcn_softlims_formulation(om, mpopt, args)
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
mpc = om.get_mpc();
[baseMVA, bus, branch] = deal(mpc.baseMVA, mpc.bus, mpc.branch);
s = mpc.softlims;

%% cheat by sticking mpopt in om temporarily for use by int2ext callback
om.userdata.mpopt = mpopt;

%% add variables, costs, and constraints

%%%%-------- limits that are the same for DC and AC formulation -----
for prop = fieldnames(s).'
    if strcmp(s.(prop{:}).hl_mod, 'none')
        continue
    end
    varname = ['s_', lower(prop{:})];
    cstname = ['cs_', lower(prop{:})];
    ns = length(s.(prop{:}).idx); % number of softlims

    %%%%%% variable and cost
    if ismember(prop{:}, {'VMIN', 'VMAX'})
        om.add_var(varname, ns, zeros(ns, 1), zeros(ns, 1), s.(prop{:}).ub ); %% add variable
        Cw = s.(prop{:}).cost(:,1);  % cost in $/pu
        om.add_quad_cost(cstname, [], Cw, 0, {varname});
    elseif ismember(prop{:}, {'RATE_A', 'PMIN', 'PMAX', 'QMIN', 'QMAX'})
        om.add_var(varname, ns, zeros(ns, 1), zeros(ns, 1), s.(prop{:}).ub / mpc.baseMVA ); %% add variable
        Cw = s.(prop{:}).cost(:,1) * mpc.baseMVA;  % cost in $/MW -> $/pu
        om.add_quad_cost(cstname, [], Cw, 0, {varname});
    elseif ismember(prop{:}, {'ANGMIN', 'ANGMAX'})
        om.add_var(varname, ns, zeros(ns, 1), zeros(ns, 1), s.(prop{:}).ub * pi/180 ); %% add variable
        Cw = s.(prop{:}).cost(:,1) * 180/pi;  % cost in $/deg -> $/rad
        om.add_quad_cost(cstname, [], Cw, 0, {varname});
    else
        error('userfcn_soflims_formulation: woops! property %s is unknown ', prop{:})
    end

    %%%%%% constraints
    if strcmp(prop{:}, 'ANGMIN')
        %%% theta_f - theta_t + s_angmin >= s.ANGMIN.sav
        ns = length(s.ANGMIN.idx);
        Av = sparse([1:ns,1:ns].', [mpc.branch(s.ANGMIN.idx, F_BUS); mpc.branch(s.ANGMIN.idx, T_BUS)], [ones(ns,1);-ones(ns,1)], ns, size(mpc.bus,1));
        As = speye(ns);
        lb = s.ANGMIN.sav * pi/180;
        ub = Inf(ns,1);

        om.add_lin_constraint('soft_angmin', [Av As], lb, ub, {'Va', 's_angmin'});
    end
    if strcmp(prop{:}, 'ANGMAX')
        %%% theta_f - theta_t - s_angmax <= s.ANGMAX.sav
        ns = length(s.ANGMAX.idx);
        Av = sparse([1:ns,1:ns].', [mpc.branch(s.ANGMAX.idx, F_BUS); mpc.branch(s.ANGMAX.idx, T_BUS)], [ones(ns,1);-ones(ns,1)], ns, size(mpc.bus,1));
        As = speye(ns);
        lb = -Inf(ns,1);
        ub = s.ANGMAX.sav * pi/180;

        om.add_lin_constraint('soft_angmax', [Av -As], lb, ub, {'Va', 's_angmax'});
    end
    if strcmp(prop{:}, 'PMIN')
        %%% Pg + s_pmin >= s.PMIN.sav
        ns = length(s.PMIN.idx);
        Av = sparse(1:ns, s.PMIN.idx, 1, ns, size(mpc.gen,1));
        As = speye(ns);
        lb = s.PMIN.sav / mpc.baseMVA;
        ub = Inf(ns,1);

        om.add_lin_constraint('soft_pmin', [Av As], lb, ub, {'Pg', 's_pmin'});
    end
    if strcmp(prop{:}, 'PMAX')
        %%% Pg - s_pmax <= s.PMAX.sav
        ns = length(s.PMAX.idx);
        Av = sparse(1:ns, s.PMAX.idx, 1, ns, size(mpc.gen,1));
        As = speye(ns);
        lb = -Inf(ns,1);
        ub = s.PMAX.sav / mpc.baseMVA;

        om.add_lin_constraint('soft_pmax', [Av -As], lb, ub, {'Pg', 's_pmax'});
    end
end

if strcmp(mpopt.model, 'DC') && ~strcmp(s.RATE_A.hl_mod, 'none')
    ns  = length(s.RATE_A.idx);
    %%% fetch Bf matrix for DC model
    Bf = om.get_userdata('Bf');
    Pfinj = om.get_userdata('Pfinj');

    %%% form constraints
    %%%    Bf * Va - s_rate_a <= -Pfinj + Pfmax
    %%%   -Bf * Va - s_rate_a <=  Pfinj + Pfmax
    I = speye(ns, ns);
    Asf = [ Bf(s.RATE_A.idx, :) -I];
    Ast = [-Bf(s.RATE_A.idx, :) -I];
    lsf = -Inf(ns, 1);
    lst = lsf;
    usf =  -Pfinj(s.RATE_A.idx) + s.RATE_A.sav/mpc.baseMVA ;
    ust =   Pfinj(s.RATE_A.idx) + s.RATE_A.sav/mpc.baseMVA ;

    om.add_lin_constraint('softPf',  Asf, lsf, usf, {'Va', 's_rate_a'});     %% ns
    om.add_lin_constraint('softPt',  Ast, lst, ust, {'Va', 's_rate_a'});     %% ns
else
    %% form AC constraints (see softlims_fcn() below)
    %%%%% voltage limits
    if ~strcmp(s.VMIN.hl_mod, 'none')
        %%% Vm + s_vmin >= s.VMIN.sav
        ns  = length(s.VMIN.idx);
        Av  = sparse(1:ns, s.VMIN.idx, 1, ns, size(mpc.bus,1));
        As  = speye(ns);
        lb  = s.VMIN.sav;
        ub  = Inf(ns,1);

        om.add_lin_constraint('soft_vmin', [Av As], lb, ub, {'Vm', 's_vmin'});
    end
    if ~strcmp(s.VMAX.hl_mod, 'none')
        %%% Vm - s_vmax <= s.VMAX.sav
        ns  = length(s.VMAX.idx);
        Av  = sparse(1:ns, s.VMAX.idx, 1, ns, size(mpc.bus,1));
        As  = speye(ns);
        lb  = -Inf(ns,1);
        ub  = s.VMAX.sav;

        om.add_lin_constraint('soft_vmax', [Av -As], lb, ub, {'Vm', 's_vmax'});
    end
    if ~strcmp(s.QMIN.hl_mod, 'none')
        %%% Qg + s_pmin >= s.QMIN.max
        ns  = length(s.QMIN.idx);
        Av  = sparse(1:ns, s.QMIN.idx, 1, ns, size(mpc.gen,1));
        As  = speye(ns);
        lb  = s.QMIN.sav / mpc.baseMVA;
        ub  = Inf(ns,1);

        om.add_lin_constraint('soft_qmin', [Av As], lb, ub, {'Qg', 's_qmin'});
    end
    if ~strcmp(s.QMAX.hl_mod, 'none')
        %%% Qg - s_pmax <= s.QMAX.sav
        ns  = length(s.QMAX.idx);
        Av  = sparse(1:ns, s.QMAX.idx, 1, ns, size(mpc.gen,1));
        As  = speye(ns);
        lb  = -Inf(ns,1);
        ub  = s.QMAX.sav / mpc.baseMVA;

        om.add_lin_constraint('soft_qmax', [Av -As], lb, ub, {'Qg', 's_qmax'});
    end
    if ~strcmp(s.RATE_A.hl_mod, 'none')
        %% build admittance matrices
        ns  = length(s.RATE_A.idx);
        [Ybus, Yf, Yt] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch);

        fcn = @(x)softlims_fcn(x, mpc, Yf(s.RATE_A.idx, :), Yt(s.RATE_A.idx, :), s.RATE_A.idx, mpopt, s.RATE_A.sav/mpc.baseMVA);
        hess = @(x, lam)softlims_hess(x, lam, mpc, Yf(s.RATE_A.idx, :), Yt(s.RATE_A.idx, :), s.RATE_A.idx, mpopt);
        om.add_nln_constraint({'softSf', 'softSt'}, [ns;ns], 0, fcn, hess, {'Va', 'Vm', 's_rate_a'});
    end
end


%%-----  int2ext  ------------------------------------------------------
function results = userfcn_softlims_int2ext(results, mpopt, args)
%
%   results = userfcn_softlims_int2ext(results, mpopt, args)
%
%   This is the 'int2ext' stage userfcn callback that converts everything
%   back to external indexing and packages up the results. It expects to
%   find a 'softlims' field in the results struct as described for mpc above.
%   It also expects the results to contain solved branch flows and linear
%   constraints named 'softlims' which are used to populate output fields
%   in results.softlims. The optional args are not currently used.

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

% structures used to index into softlims in a loop
lims = softlims_lim2mat();

%% get internal softlims struct and mpopt
isOPF  = isfield(results, 'f') && ~isempty(results.f);

s = results.softlims;
if isOPF
    mpopt = results.om.get_userdata('mpopt');   %% extract and remove mpopt from om
    results.om.userdata = rmfield(results.om.userdata, 'mpopt');
end

%%-----  convert stuff back to external indexing  -----
o = results.order;
nl0 = size(o.ext.branch, 1);    %% original number of branches
nl = size(results.branch, 1);   %% number of on-line branches
o.branch.status.on;
results.softlims = results.order.ext.softlims;

%%-----  restore hard limits  -----
for prop = fieldnames(s).'
    mat  = lims.(prop{:});
    if strcmp(s.(prop{:}).hl_mod, 'none')
        continue
    end
    results.(mat)(s.(prop{:}).idx, eval(prop{:})) = s.(prop{:}).sav;
end

%%-----  remove rval, and sav fields  -----
for prop = fieldnames(s).'
    if strcmp(s.(prop{:}).hl_mod, 'none')
        continue
    end
    results.softlims.(prop{:}) = rmfield(results.softlims.(prop{:}), 'sav');
    results.softlims.(prop{:}) = rmfield(results.softlims.(prop{:}), 'rval');
end
%%-----  results post-processing  -----
%% get overloads and overload costs
if isOPF
    for prop = fieldnames(s).'
        mat  = lims.(prop{:});
        if strcmp(s.(prop{:}).hl_mod, 'none')
            continue
        end
        n0 = size(o.ext.(mat), 1);    %% original number
        n  = size(results.(mat), 1);  %% number on-line
        results.softlims.(prop{:}).overload = zeros(n0, 1);
        results.softlims.(prop{:}).ovl_cost = zeros(n0, 1);
        varname = ['s_', lower(prop{:})];  %% variable name
        if ismember(prop{:}, {'VMIN', 'VMAX'})
            if ~strcmp(mpopt.model, 'AC')
                continue
            end
            var = results.var.val.(varname); %stays in p.u.
            var(var < 1e-8) = 0;
        elseif ismember(prop{:}, {'RATE_A', 'PMIN', 'PMAX', })
            var = results.var.val.(varname);
            var(var < 1e-8) = 0;
            var = var * results.baseMVA; % p.u. -> MW
        elseif ismember(prop{:}, {'QMIN', 'QMAX', })
            if ~strcmp(mpopt.model, 'AC')
                continue
            end
            var = results.var.val.(varname);
            var(var < 1e-8) = 0;
            var = var * results.baseMVA; % p.u. -> MVAr
        elseif ismember(prop{:}, {'ANGMIN', 'ANGMAX'})
            var = results.var.val.(varname);
            var(var < 1e-8) = 0;
            var = var * 180/pi; % rad -> deg
        else
            error('userfcn_soflims_formulation: woops! property %s is unknown ', prop{:})
        end
    %     var(var < 1e-8) = 0;
        % NOTE: o.(mat).status.on is a vector nx1 where n is the INTERNAL number of
        % elements. The entries are the EXTERNAL locations (row numbers).
        if strcmp(mat, 'gen')
            results.softlims.(prop{:}).overload(o.(mat).status.on(o.gen.i2e(s.(prop{:}).idx))) = var;
            results.softlims.(prop{:}).ovl_cost(o.(mat).status.on(o.gen.i2e(s.(prop{:}).idx))) = var .* s.(prop{:}).cost(:,1);
        else
            results.softlims.(prop{:}).overload(o.(mat).status.on(s.(prop{:}).idx)) = var;
            results.softlims.(prop{:}).ovl_cost(o.(mat).status.on(s.(prop{:}).idx)) = var .* s.(prop{:}).cost(:,1);
        end
    end

    %% get shadow prices
    if ~strcmp(s.ANGMAX.hl_mod, 'none')
        results.branch(s.ANGMAX.idx, MU_ANGMAX) = results.lin.mu.u.soft_angmax * pi/180;
    end
    if ~strcmp(s.ANGMIN.hl_mod, 'none')
        results.branch(s.ANGMIN.idx, MU_ANGMIN) = results.lin.mu.l.soft_angmin * pi/180;
    end
    if ~strcmp(s.PMAX.hl_mod, 'none')
        results.gen(s.PMAX.idx, MU_PMAX) = results.lin.mu.u.soft_pmax / results.baseMVA;
    end
    if ~strcmp(s.PMIN.hl_mod, 'none')
        results.gen(s.PMIN.idx, MU_PMIN) = results.lin.mu.l.soft_pmin / results.baseMVA;
    end
    if strcmp(mpopt.model, 'DC')
        if ~strcmp(s.RATE_A.hl_mod, 'none')
            results.branch(s.RATE_A.idx, MU_SF) = results.lin.mu.u.softPf / results.baseMVA;
            results.branch(s.RATE_A.idx, MU_ST) = results.lin.mu.u.softPt / results.baseMVA;

            %if results.success    %% double-check value of overloads being returned (if solution was successful)
            %    vv = results.om.get_idx();
            %    check1 = zeros(nl0, 1);
            %    check1(o.branch.status.on(s.RATE_A.idx)) = results.x(vv.i1.s_rate_a:vv.iN.s_rate_a) * results.baseMVA;
            %    check2 = zeros(nl0, 1);
            %    k = find(results.branch(:, RATE_A) & ...
            %             abs(results.branch(:, PF)) > results.branch(:, RATE_A) );
            %    check2(o.branch.status.on(k)) = ...
            %            abs(results.branch(k, PF)) - results.branch(k, RATE_A);
            %    err1 = norm(results.softlims.RATE_A.overload-check1);
            %    err2 = norm(results.softlims.RATE_A.overload-check2);
            %    errtol = 1e-4;
            %    if err1 > errtol || err2 > errtol
            %        [ results.softlims.RATE_A.overload check1 results.softlims.RATE_A.overload-check1 ]
            %        [ results.softlims.RATE_A.overload check2 results.softlims.RATE_A.overload-check2 ]
            %        error('userfcn_softlims_int2ext: problem with consistency of overload values');
            %    end
            %end
        end
    else %AC model
        if ~strcmp(s.VMAX.hl_mod, 'none')
            results.bus(s.VMAX.idx, MU_VMAX) = results.lin.mu.u.soft_vmax;
        end
        if ~strcmp(s.VMIN.hl_mod, 'none')
            results.bus(s.VMIN.idx, MU_VMIN) = results.lin.mu.l.soft_vmin;
        end
        if ~strcmp(s.QMAX.hl_mod, 'none')
            results.gen(s.QMAX.idx, MU_QMAX) = results.lin.mu.u.soft_qmax / results.baseMVA;
        end
        if ~strcmp(s.QMIN.hl_mod, 'none')
            results.gen(s.QMIN.idx, MU_QMIN) = results.lin.mu.l.soft_qmin / results.baseMVA;
        end
        if ~strcmp(s.RATE_A.hl_mod, 'none')
            if upper(mpopt.opf.flow_lim(1)) == 'P'
                results.branch(s.RATE_A.idx, MU_ST) = results.nli.mu.softSf / results.baseMVA;
                results.branch(s.RATE_A.idx, MU_SF) = results.nli.mu.softSt / results.baseMVA;
            else
                        var = results.var.val.s_rate_a * results.baseMVA;
                        var(var < 1e-8) = 0;
                %% conversion factor for squared constraints (2*F)
                cf = 2 * (s.RATE_A.sav + var)/results.baseMVA;
                %cf = 2 * (s.Pfmax + flv / results.baseMVA);
                results.branch(s.RATE_A.idx, MU_ST) = results.nli.mu.softSf .* cf / results.baseMVA;
                results.branch(s.RATE_A.idx, MU_SF) = results.nli.mu.softSt .* cf / results.baseMVA;
            end

            %if results.success    %% double-check value of overloads being returned (if solution was successful)
            %    vv = results.om.get_idx();
            %    check1 = zeros(nl0, 1);
            %    check1(o.branch.status.on(s.RATE_A.idx)) = results.x(vv.i1.s_rate_a:vv.iN.s_rate_a) * results.baseMVA;
            %    err1 = norm(results.softlims.RATE_A.overload-check1);
            %    errtol = 1e-4;
            %    if err1 > errtol
            %        [ results.softlims.RATE_A.overload check1 results.softlims.RATE_A.overload-check1 ]
            %        error('userfcn_softlims_int2ext: problem with consistency of overload values');
            %    end
            %end
        end
    end
end
%%----- swap bus idx vectors -----
% set s.VMIN.idx and s.VMAX.idx to include external bus numbers (stored in
% busnum. Store row indices in rowidx to be used by the print function.
for prop = {'VMAX', 'VMIN'}
    if ~strcmp(s.(prop{:}).hl_mod, 'none')
        results.softlims.(prop{:}).idx = s.(prop{:}).busnum;
        results.softlims.(prop{:}).rowidx = s.(prop{:}).idx;
        results.softlims.(prop{:}) = rmfield(results.softlims.(prop{:}), 'busnum');
    end
end

%%-----  printpf  ------------------------------------------------------
function results = userfcn_softlims_printpf(results, fd, mpopt, args)
%
%   results = userfcn_softlims_printpf(results, fd, mpopt, args)
%
%   This is the 'printpf' stage userfcn callback that pretty-prints the
%   results. It expects a results struct, a file descriptor and a MATPOWER
%   options struct. The optional args are not currently used.

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
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
    s = results.softlims;
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n|     Soft Limits                                                              |');
    fprintf(fd, '\n================================================================================');
    if ~strcmp(s.RATE_A.hl_mod, 'none')
        k = find(s.RATE_A.overload(s.RATE_A.idx) | sum(results.branch(s.RATE_A.idx, MU_SF:MU_ST), 2) > ptol);
        fprintf(fd, '\nFlow Limits:');
        fprintf(fd, '\n----------------------------------------');
        if isempty(k)
            fprintf(fd,'\nNo violations.\n');
        else
            fprintf(fd, '\nBrnch   From   To     Flow      Limit    Overload     mu');
            fprintf(fd, '\n  #     Bus    Bus    (MW)      (MW)       (MW)     ($/MW)');
            fprintf(fd, '\n-----  -----  -----  --------  --------  --------  ---------');
            fprintf(fd, '\n%4d%7d%7d%10.2f%10.2f%10.2f%11.3f', ...
                    [   s.RATE_A.idx(k), results.branch(s.RATE_A.idx(k), [F_BUS, T_BUS]), ...
                        results.branch(s.RATE_A.idx(k), [PF, RATE_A]), ...
                        s.RATE_A.overload(s.RATE_A.idx(k)), ...
                        sum(results.branch(s.RATE_A.idx(k), MU_SF:MU_ST), 2) ...
                    ]');
            fprintf(fd, '\n                                         --------');
            fprintf(fd, '\n                                Total:%10.2f', ...
                    sum(s.RATE_A.overload(s.RATE_A.idx(k))));
            fprintf(fd, '\n');
        end
    end
    if ~strcmp(s.VMAX.hl_mod,'none') && strcmp(mpopt.model, 'AC')
        k = find(s.VMAX.overload(s.VMAX.rowidx) | results.bus(s.VMAX.rowidx, MU_VMAX) > ptol);
        fprintf(fd, '\nMaximum Voltage Magnitude Limits:');
        fprintf(fd, '\n----------------------------------------');
        if isempty(k)
            fprintf(fd,'\nNo violations.\n');
        else
            fprintf(fd, '\nBus    Voltage   Limit   Overload    mu');
            fprintf(fd, '\n  #    Mag(pu)   (pu)     (pu)     ($/pu)');
            fprintf(fd, '\n-----  -------  -------  -------  ---------');
            fprintf(fd, '\n%5d%8.3f%9.3f%9.3f%11.3f',...
                [ s.VMAX.idx(k), results.bus(s.VMAX.rowidx(k),[VM, VMAX]),...
                  s.VMAX.overload(s.VMAX.rowidx(k)), ...
                  results.bus(s.VMAX.rowidx(k), MU_VMAX)...
                ]');
            fprintf(fd, '\n                        --------');
            fprintf(fd, '\n               Total:%10.2f', ...
                    sum(s.VMAX.overload(s.VMAX.rowidx(k))));
            fprintf(fd, '\n');
        end
    end
    if ~strcmp(s.VMIN.hl_mod,'none') && strcmp(mpopt.model, 'AC')
        k = find(s.VMIN.overload(s.VMIN.rowidx) | results.bus(s.VMIN.rowidx, MU_VMIN) > ptol);
        fprintf(fd, '\nMinimum Voltage Magnitude Limits:');
        fprintf(fd, '\n----------------------------------------');
        if isempty(k)
            fprintf(fd,'\nNo violations.\n');
        else
            fprintf(fd, '\n----------------------------------------');
            fprintf(fd, '\n Bus   Voltage   Limit   Overload    mu');
            fprintf(fd, '\n  #    Mag(pu)   (pu)     (pu)     ($/pu)');
            fprintf(fd, '\n-----  -------  -------  -------  ---------');
            fprintf(fd, '\n%5d%8.3f%9.3f%9.3f%11.3f',...
                [ s.VMIN.idx(k), results.bus(s.VMIN.rowidx(k),[VM, VMIN]),...
                  s.VMIN.overload(s.VMIN.rowidx(k)), ...
                  results.bus(s.VMIN.rowidx(k), MU_VMIN)...
                ]');
            fprintf(fd, '\n                        --------');
            fprintf(fd, '\n               Total:%10.2f', ...
                    sum(s.VMIN.overload(s.VMIN.rowidx(k))));
            fprintf(fd, '\n');
        end
    end
    if ~strcmp(s.PMAX.hl_mod,'none')
        k = find(s.PMAX.overload(s.PMAX.idx) | results.gen(s.PMAX.idx, MU_PMAX) > ptol);
        fprintf(fd, '\nMaximum Generator P Limits:');
        fprintf(fd, '\n----------------------------------------');
        if isempty(k)
            fprintf(fd,'\nNo violations.\n');
        else
            fprintf(fd, '\nGen     Bus  Generation  Limit   Overload    mu');
            fprintf(fd, '\n  #      #     P (MW)    (MW)     (MW)     ($/MW)');
            fprintf(fd, '\n-----  -----  --------  -------  -------  ---------');
            fprintf(fd, '\n%5d%7d%9.2f%9.2f%9.3f%11.3f',...
                [ s.PMAX.idx(k), results.gen(s.PMAX.idx(k),[GEN_BUS, PG, PMAX]),...
                  s.PMAX.overload(s.PMAX.idx(k)), ...
                  results.gen(s.PMAX.idx(k), MU_PMAX)...
                ]');
            fprintf(fd, '\n                                --------');
            fprintf(fd, '\n                       Total:%10.2f', ...
                    sum(s.PMAX.overload(s.PMAX.idx(k))));
            fprintf(fd, '\n');
        end
    end
    if ~strcmp(s.PMIN.hl_mod,'none')
        k = find(s.PMIN.overload(s.PMIN.idx) | results.gen(s.PMIN.idx, MU_PMIN) > ptol);
        fprintf(fd, '\nMinimum Generator P Limits:');
        fprintf(fd, '\n----------------------------------------');
        if isempty(k)
            fprintf(fd,'\nNo violations.\n');
        else
            fprintf(fd, '\nGen     Bus  Generation  Limit   Overload    mu');
            fprintf(fd, '\n  #      #     P (MW)    (MW)     (MW)     ($/MW)');
            fprintf(fd, '\n-----  -----  --------  -------  -------  ---------');
            fprintf(fd, '\n%5d%7d%9.2f%9.2f%9.3f%11.3f',...
                [ s.PMIN.idx(k), results.gen(s.PMIN.idx(k),[GEN_BUS, PG, PMIN]),...
                  s.PMIN.overload(s.PMIN.idx(k)), ...
                  results.gen(s.PMIN.idx(k), MU_PMIN)...
                ]');
            fprintf(fd, '\n                                --------');
            fprintf(fd, '\n                       Total:%10.2f', ...
                    sum(s.PMIN.overload(s.PMIN.idx(k))));
            fprintf(fd, '\n');
        end
    end
    if ~strcmp(s.QMAX.hl_mod,'none') && strcmp(mpopt.model, 'AC')
        k = find(s.QMAX.overload(s.QMAX.idx) | results.gen(s.QMAX.idx, MU_QMAX) > ptol);
        fprintf(fd, '\nMaximum Generator Q Limits:');
        fprintf(fd, '\n----------------------------------------');
        if isempty(k)
            fprintf(fd,'\nNo violations.\n');
        else
            fprintf(fd, '\nGen     Bus  Generation  Limit   Overload    mu');
            fprintf(fd, '\n  #      #    Q (MVAr)  Q (MVAr)  (MVAr)   ($/MVAr)');
            fprintf(fd, '\n-----  -----  --------  -------  -------  ---------');
            fprintf(fd, '\n%5d%7d%9.2f%9.2f%9.3f%11.3f',...
                [ s.QMAX.idx(k), results.gen(s.QMAX.idx(k),[GEN_BUS, QG, QMAX]),...
                  s.QMAX.overload(s.QMAX.idx(k)), ...
                  results.gen(s.QMAX.idx(k), MU_QMAX)...
                ]');
            fprintf(fd, '\n                                --------');
            fprintf(fd, '\n                       Total:%10.2f', ...
                    sum(s.QMAX.overload(s.QMAX.idx(k))));
            fprintf(fd, '\n');
        end
    end
    if ~strcmp(s.QMIN.hl_mod,'none') && strcmp(mpopt.model, 'AC')
        k = find(s.QMIN.overload(s.QMIN.idx) | results.gen(s.QMIN.idx, MU_QMIN) > ptol);
        fprintf(fd, '\nMinimum Generator Q Limits:');
        fprintf(fd, '\n----------------------------------------');
        if isempty(k)
            fprintf(fd,'\nNo violations.\n');
        else
            fprintf(fd, '\nGen     Bus  Generation  Limit   Overload    mu');
            fprintf(fd, '\n  #      #    Q (MVAr)  Q (MVAr)  (MVAr)   ($/MVAr)');
            fprintf(fd, '\n-----  -----  --------  -------  -------  ---------');
            fprintf(fd, '\n%5d%7d%9.2f%9.2f%9.3f%11.3f',...
                [ s.QMIN.idx(k), results.gen(s.QMIN.idx(k),[GEN_BUS, QG, QMIN]),...
                  s.QMIN.overload(s.QMIN.idx(k)), ...
                  results.gen(s.QMIN.idx(k), MU_QMIN)...
                ]');
            fprintf(fd, '\n                                --------');
            fprintf(fd, '\n                       Total:%10.2f', ...
                    sum(s.QMIN.overload(s.QMIN.idx(k))));
            fprintf(fd, '\n');
        end
    end
    if ~strcmp(s.ANGMAX.hl_mod,'none')
        k = find(s.ANGMAX.overload(s.ANGMAX.idx) | results.branch(s.ANGMAX.idx, MU_ANGMAX) > ptol);
        delta = calc_branch_angle(results);
        fprintf(fd, '\nMaximum Angle Difference Limits:');
        fprintf(fd, '\n----------------------------------------');
        if isempty(k)
            fprintf(fd,'\nNo violations.\n');
        else
            fprintf(fd, '\nBrnch   From   To     Angle     Limit    Overload     mu');
            fprintf(fd, '\n  #     Bus    Bus    (deg)     (deg)     (deg)     ($/MW)');
            fprintf(fd, '\n-----  -----  -----  --------  --------  --------  ---------');
            fprintf(fd, '\n%4d%7d%7d%10.3f%10.3f%10.3f%11.3f', ...
                [ s.ANGMAX.idx(k), results.branch(s.ANGMAX.idx(k), [F_BUS, T_BUS]), ...
                  delta(s.ANGMAX.idx(k)), results.branch(s.ANGMAX.idx(k), ANGMAX),...
                  s.ANGMAX.overload(s.ANGMAX.idx(k)),...
                  results.branch(s.ANGMAX.idx(k), MU_ANGMAX)...
                ]');
            fprintf(fd, '\n                                         --------');
            fprintf(fd, '\n                                Total:%10.2f', ...
                    sum(s.ANGMAX.overload(s.ANGMAX.idx(k))));
            fprintf(fd, '\n');
        end
    end
    if ~strcmp(s.ANGMIN.hl_mod,'none')
        k = find(s.ANGMIN.overload(s.ANGMIN.idx) | results.branch(s.ANGMIN.idx, MU_ANGMIN) > ptol);
        delta = calc_branch_angle(results);
        fprintf(fd, '\nMinimum Angle Difference Limits:');
        fprintf(fd, '\n----------------------------------------');
        if isempty(k)
            fprintf(fd,'\nNo violations.\n');
        else
            fprintf(fd, '\nBrnch   From   To     Angle     Limit    Overload     mu');
            fprintf(fd, '\n  #     Bus    Bus    (deg)     (deg)     (deg)     ($/MW)');
            fprintf(fd, '\n-----  -----  -----  --------  --------  --------  ---------');
            fprintf(fd, '\n%4d%7d%7d%10.3f%10.3f%10.3f%11.3f', ...
                [ s.ANGMIN.idx(k), results.branch(s.ANGMIN.idx(k), [F_BUS, T_BUS]), ...
                  delta(s.ANGMIN.idx(k)), results.branch(s.ANGMIN.idx(k), ANGMIN),...
                  s.ANGMIN.overload(s.ANGMIN.idx(k)),...
                  results.branch(s.ANGMIN.idx(k), MU_ANGMIN)...
                ]');
            fprintf(fd, '\n                                         --------');
            fprintf(fd, '\n                                Total:%10.2f', ...
                    sum(s.ANGMIN.overload(s.ANGMIN.idx(k))));
            fprintf(fd, '\n');
        end
    end
end


%%-----  savecase  -----------------------------------------------------
function mpc = userfcn_softlims_savecase(mpc, fd, prefix, args)
%
%   mpc = userfcn_softlims_savecase(mpc, fd, prefix, args)
%
%   This is the 'savecase' stage userfcn callback that prints the M-file
%   code to save the 'softlims' field in the case file. It expects a
%   MATPOWER case struct (mpc), a file descriptor and variable prefix
%   (usually 'mpc.'). The optional args are not currently used.

% structures used to index into softlims in a loop
lims = softlims_lim2mat();

% convenience structure for the different fields
fields = struct('hl_mod', struct('desc','Hard limit modification', 'tok','%s'),...
    'hl_val', struct('desc', 'New hard limit value', 'tok', '%g'),...
    'busnum', struct('desc','bus ids for where soft voltage limit is applied','tok', '%d'),...
    'rowidx', struct('desc','bus matrix row indices', 'tok', '%d'),...
    'idx', struct('desc','%s matrix indices','tok','%d'),...
    'cost', struct('desc','violation cost coefficient','tok','%g'),...
    'ub', struct('desc','Slack variable upper bound','tok','%g'),...
    'sav', struct('desc', 'saved hard limits', 'tok', '%g'),...
    'rval',struct('desc', 'replacment value to deactivate hard limit', 'tok', '%g'),...
    'overload', struct('desc', 'overload (value above hard limit)', 'tok', '%0.10g'),...
    'ovl_cost', struct('desc', 'overload cost', 'tok', '%.10g'));

if isfield(mpc, 'softlims')
    s = mpc.softlims;

    fprintf(fd, '\n%%%%-----  Soft Limit Data  -----%%%%\n');
    
    for prop = fieldnames(lims).'
        fprintf(fd,'\n%%%% Property: %s %%%%\n', prop{:});
        for f = fieldnames(fields).'
            f = f{1};
            if isfield(s.(prop{:}),f)
                if strcmp(f, 'idx')
                    if ismember(prop{:}, {'VMIN', 'VMAX'})
                        desc = fields.busnum.desc;
                    else
                        desc = sprintf(fields.idx.desc, lims.(prop{:}));
                    end
                else
                    desc = fields.(f).desc;
                end
                fprintf(fd,'%% %s\n',desc);
                if strcmp(f,'hl_mod')
                    fprintf(fd, '%ssoftlims.%s.hl_mod = ''%s'';\n', prefix,prop{:}, s.(prop{:}).hl_mod);
                else
                    fprintf(fd, '%ssoftlims.%s.%s = [\n', prefix,prop{:},f);
                    fprintf(fd, ['\t',fields.(f).tok,';\n'], s.(prop{:}).(f));
                    fprintf(fd, '];\n\n');
                end
            end
        end
    end
end

%%-----  softlims_lim2mat  --------------------------------------------
function lim2mat = softlims_lim2mat()
lim2mat = struct(...
    'VMAX', 'bus', ...
    'VMIN','bus', ...
    'ANGMAX', 'branch', ...
    'ANGMIN', 'branch', ...
    'RATE_A', 'branch', ...
    'PMAX', 'gen', ...
    'PMIN', 'gen', ...
    'QMAX', 'gen', ...
    'QMIN', 'gen' ...
);

%%-----  softlims_defaults  --------------------------------------------
function s = softlims_defaults(mpc, mpopt)
%
%   s = softlims_defaults(mpc, mpopt)
%
%   Returns a softlims field for mpc in which any missing inputs have
%   been filled in with defaults, so that each limit type has the all
%   4 input fields (idx, cost, hl_mod, hl_cost) and idx has been expanded
%   to a vector, unless hl_mod == 'none', in which case the others need
%   not be present (and are, in any case, ignored).

%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% initialization
lims = softlims_lim2mat();
s = mpc.softlims;
if isempty(mpopt)
    warning('softlims_defaults: Assuming ''mpopt.opf.softlims.default'' = 1, since mpopt was not provided.');
    use_default = 1;
else
    use_default = mpopt.opf.softlims.default;
end 

%% get maximum generator marginal cost
% since mpc.gen and not mpc.order.ext.gen is used this is the maximum of
% ONLINE generators.
max_gen_cost = max(margcost(mpc.gencost, mpc.gen(:, PMAX)));

%% set defaults for each element
for prop = fieldnames(lims).'
    mat  = lims.(prop{:});
    if isfield(s, prop{:}) && ~isempty(s.(prop{:}))
        s.(prop{:}) = softlims_element_defaults(s.(prop{:}), prop{:}, mpc.order.ext.(mat), max_gen_cost);
    else
        if use_default
            %% pass an empty struct to assign defaults
            s.(prop{:}) = softlims_element_defaults(struct(), prop{:}, mpc.order.ext.(mat), max_gen_cost);
        else
            s.(prop{:}).hl_mod = 'none';
        end
    end
end


%%-----  softlims_init  --------------------------------------------
function s = softlims_init(mpc, mpopt)
%
%   s = softlims_init(mpc, mpopt)
%
%   Returns a softlims field for mpc in which any missing inputs have
%   been filled in with defaults, so that each limit type has the all
%   4 input fields (idx, cost, hl_mod, hl_cost) and idx has been expanded
%   to a vector, unless hl_mod == 'none', in which case the others need
%   not be present (and are, in any case, ignored).

%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% initialization
lims = softlims_lim2mat();
s = mpc.softlims;
if isempty(mpopt)
    warning('softlims_init: Assuming ''mpopt.opf.softlims.default'' = 1, since mpopt was not provided.');
    use_default = 1;
else
    use_default = mpopt.opf.softlims.default;
end 

%% get maximum generator marginal cost
% since mpc.gen and not mpc.order.ext.gen is used this is the maximum of
% ONLINE generators.
max_gen_cost = max(margcost(mpc.gencost, mpc.gen(:, PMAX)));

%% set defaults for each element
for prop = fieldnames(lims).'
    mat  = lims.(prop{:});
    if isfield(s, prop{:}) && ~isempty(s.(prop{:}))
        s.(prop{:}) = softlims_element_init(s.(prop{:}), prop{:}, mpc.order.ext.(mat), max_gen_cost);
    else
        if use_default
            %% pass an empty struct to assign defaults
            s.(prop{:}) = softlims_element_init(struct(), prop{:}, mpc.order.ext.(mat), max_gen_cost);
        else
            s.(prop{:}).hl_mod = 'none';
        end
    end
end


%%-----  softlims_element_defaults  --------------------------------------------
function s = softlims_element_defaults(s, prop, mat, max_gen_cost)
%
%   s = softlims_element_defaults(s, prop, mat, max_gen_cost)
%
%   For each property, we fill in default values as needed for:
%       idx: index of affected buses, branches, or generators (bus numbers
%           rather than indexes for bus)
%       cost: linear cost on corresponding violation variables
%       hl_mod: string value specifying what to do with original hard limits
%               for this property
%           'none'    : leave hard limits unchanged, do not apply soft limits
%           'remove'  : eliminate hard limits completely
%           'replace' : replace hard limit with value specified in 'hl_val'
%           'scale'   : multiply hard limit by value specified in 'hl_val'
%           'shift'   : add value specified in 'hl_val' to hard limit limit
%       hl_val: value used to set/modify hard limit as specified by 'hl_mod'

%% ignore if hl_mod none specified
if isfield(s, 'hl_mod') && strcmp(s.hl_mod, 'none')
    return
end

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% default cost is default_cost * cost_scale.(prop)
default_cost  = 1000;
cost_scale = struct( ...
    'VMAX', 100, 'VMIN', 100, ...   %% $/pu
    'ANGMAX', 1, 'ANGMIN', 1, ...   %% $/deg
    'RATE_A', 1, ...                %% $/MW
    'PMIN', 1, 'PMAX', 1, ...       %% $/MW
    'QMIN', 1, 'QMAX', 1 );         %% $/MVar

if nargin < 4
    max_gen_cost = 0;
end

if ~ismember(prop, {'VMAX', 'VMIN', 'RATE_A', 'PMAX', 'PMIN', 'QMAX', 'QMIN', 'ANGMAX', 'ANGMIN'})
    error('softlims_element_defaults: Unknown limit type ''%s''', prop)
end

%% idx: indicies of bounds to relax
switch prop
    case {'VMAX', 'VMIN'}
        % Initial list of candidates for voltage limits contains all buses
        % NOTE: idxfull contains EXTERNAL bus numbers
        idxfull = mat(:, BUS_I);
    case {'ANGMAX', 'ANGMIN'}
        % Initial list of candidates for branch angle differences contains all
        % active branches with a meaningful limit (not 0, +360 or -360).
        % NOTE: idxfull contains locations in EXTERNAL branch list
        idxfull = find(mat(:, BR_STATUS) > 0 & mat(:, eval(prop)) & abs(mat(:, eval(prop))) < 360 );
    case 'RATE_A'
        % Initial list of candidates for line rating are all active branches
        % with rating not set to infinity (0).
        % NOTE: idxfull contains locations in EXTERNAL branch list
        idxfull = find(mat(:, BR_STATUS) > 0 & mat(:, RATE_A) > 0);
    case {'PMAX', 'QMAX', 'QMIN'}
        % Initial list of candidates for generation limits (except Pmin) are
        % all active generators whose limits are not infinity
        % NOTE: idxfull contains locations in EXTERNAL generator list
        idxfull = find( (mat(:, GEN_STATUS) > 0) & ~isload(mat) & ~isinf(mat(:, eval(prop))) );
    case 'PMIN'
        % Initial list of candidates for Pmin are all active generators that
        % are NOT dispachable loads
        % NOTE: idxfull contains locations in EXTERNAL generator list
        idxfull = find( (mat(:, GEN_STATUS) > 0) & ~isload(mat) );
end

if isfield(s, 'idx')
%     % idxmask is a boolean vector the size of s.idx. Entry i of
%     % idxmask is 1 if s.idx(i) is in idxfull and 0 otherwise.
    if size(s.idx, 2) > 1
        s.idx = s.idx.';
        if size(s.idx, 2) > 1
            error('softlim_defaults: mpc.softlims.%s.idx must be a vector', prop)
        end
    end
%     idxmask = ismember(s.idx, idxfull);
%     s.idx = s.idx(idxmask); %remove possibly irrelevant entries entered by user
else
    % if no indices specified by user use the full list idxfull.
%     idxmask = true(size(idxfull));
    s.idx = idxfull;
end

if ismember(prop, {'VMAX', 'VMIN'})
    % for consistency between all the different limits, we want s.idx for
    % buses to contain locations rather than external bus numbers.
    % External bus numbers are stored in s.busnum and s.idx is rewritten to
    % include the locations (rows) of those buses.
    s.busnum = s.idx;
%     s.idx = find(ismember(mat(:,BUS_I), s.idx));
    s_idx = find(ismember(mat(:,BUS_I), s.idx));
else
    s_idx = s.idx;
end

% if there are no constraints to relax, set hl_mod to none and exit
if isempty(s_idx)
    s.hl_mod = 'none';
    return
end

% %% save original values in s.sav
% s.sav = mat(s.idx, eval(prop));     %% s.idx is row index into relevant matrix

%% cost: cost of violating softlim
if isfield(s, 'cost') && ~isempty(s.cost)
    if isscalar(s.cost)
        s.cost = s.cost * ones(size(s_idx));
%     else    % vector: apply idxmask filter in case elements removed in idx stage
%         try
%             s.cost = s.cost(idxmask);
%         catch ME
%             warning('softlims_element_defaults: something went wrong when handling the cost for property %s. Perhaps the size of the ''cost'' vector didn''t match the size of the ''idx'' vector?', prop)
%             rethrow(ME)
%         end
    end
else
    % default cost is the maximum between a predefined value and
    % 2 times the maximum generation marginal cost. A scaling is applied
    % to account for the fact that we are essentially equating 1MW with
    % 0.01 p.u. voltage or 1 deg angle difference.
    
    ctmp = cost_scale.(prop) * max(default_cost, 2*max_gen_cost);
    s.cost = ctmp * ones(size(s_idx));
end

%% type of limit and upper bound of slack variable
ubsign = struct('VMAX', 1 , 'VMIN', -1, 'RATE_A', 1, ...
                'PMAX', 1, 'PMIN', -1, 'QMAX', 1, 'QMIN', -1, ...
                'ANGMAX', 1, 'ANGMIN', -1);
shift_defaults = struct('VMAX', 0.25 , 'VMIN', 0.25, 'RATE_A', 10, ...
                'PMAX', 10, 'PMIN', 10, 'QMAX', 10, 'QMIN', 10, ...
                'ANGMAX', 10, 'ANGMIN', 10);
if isfield(s, 'hl_mod')
    switch s.hl_mod
        case 'remove'
            % slack upper bound and new limit are both infinite
            s.hl_val = ubsign.(prop)*Inf;
%             s.ub = Inf(size(s.idx));
        case 'scale'
            % new hard limit is hl_val * original limit
            % for softlims in a positive direction the upper limit is:
            %           s.ub = (hl_val) * original limit - original limit
            %                = (hl_val - 1) * original limit
            % for softlims in a negative direction the upper limit is
            %           s.ub = original limit - hl_val*original limit
            %                = (1 - hl_val)*original limit
%             switch vectorcheck(s, idxmask)
%                 case 0
%                     error('softlims_element_defaults: provided hl_val vector does not conform to idx vector. When specifying hl_val as a vector idx must also be explicitly given') 
%                 case 1 % vector of values
%                     s.ub = ubsign.(prop)*(s.hl_val(idxmask) - 1).*mat(s.idx, eval(prop));
%                 case 2 % scalar value
%                     s.ub = ubsign.(prop)*(s.hl_val - 1)*mat(s.idx, eval(prop));
%                 case 3 % no hl_val specified, use default of 2 or 1/2
                if ~isfield(s, 'hl_val')
                    % use 2 if ubsign and original constraint have same sign
                    % otherwise use 1/2
                    orig_lim = mat(s_idx, eval(prop));
                    s.hl_val = 2 * ones(size(s_idx));   % scale up by 2
                    k = find(ubsign.(prop) * orig_lim < 0); % unless opp. sign
                    s.hl_val(k) = 1/2;                  % then scale down by 2
%                     s.ub = abs((s.hl_val - 1) .* orig_lim);
                end
%             end
        case 'shift'
            % new hard limit is original limit + ubsign*hl_val
            % for positive direction limits:
            %    s.ub = original limit + hl_val - original limit = hl_val
            % for negative direction limits:
            %    s.ub = original limit - (original limit - hl_val) = hl_val
            % s.ub is therefore simply hl_val
%             switch vectorcheck(s, idxmask)
%                 case 0
%                     error('softlims_element_defaults: provided hl_val vector does not conform to idx vector. When specifying hl_val as a vector idx must also be explicitly given')
%                 case 1 % vector of values
%                     s.ub = s.hl_val(idxmask);
%                 case 2 % scalar value
%                     s.ub = s.hl_val*ones(size(s.idx));
%                 case 3 % no hl_val specified use defaults in shift_defaults
                if ~isfield(s, 'hl_val')
                    s.hl_val = shift_defaults.(prop);
                end
%                     s.ub = s.hl_val*ones(size(s.idx));
%             end
        case 'replace'
%             % new hard limit is hl_val
%             % for positive direction limits:
%             %    s.ub = hl_val - original limit
%             % for negative direction limits:
%             %    s.ub = original limit - hl_val
%             % Therefore s.ub = ubsign*(hl_val - original limit)
%             switch vectorcheck(s, idxmask)
%                 case 0
%                     error('softlims_element_defaults: provided hl_val vector does not conform to idx vector. When specifying hl_val as a vector idx must also be explicitly given')
%                 case 1 % vector of values
%                     s.ub = ubsign.(prop)*(s.hl_val(idxmask) - mat(s.idx, eval(prop)));
%                 case 2 % scalar value
%                     s.ub = ubsign.(prop)*(s.hl_val - mat(s.idx, eval(prop)));
%                 case 3 % non specified
%                     error('softlims_element_defaults: for hard limit ''replace'' modification, replacement value hl_val must be specified')
%             end
        otherwise
            error('softlims_element_defaults: unknown hard limit modification %s', s.hl_mod)
    end
else
    % defaults
    switch prop
        case 'PMIN'
            % Normal generators (PMIN > 0) can only be relaxed to 0
            % If PMIN < 0 allow PMIN to go to -Inf
            s.hl_mod = 'replace';
            s.hl_val = zeros(size(s_idx));
            s.hl_val(mat(s_idx, PMIN) < 0) = -Inf;
%             s.ub     = mat(s.idx, PMIN) - s.hl_val;
        case 'VMIN'
            % VMIN can only be relaxed to zero, not further.
            s.hl_mod = 'replace';
            s.hl_val = 0;
%             s.ub     = mat(s.idx, VMIN);
        case {'QMIN', 'ANGMIN'}
            s.hl_mod = 'remove';
            s.hl_val = -Inf;
%             s.ub     = Inf(size(s.idx));
        otherwise
            s.hl_mod = 'remove';
            s.hl_val = Inf;
%             s.ub     = Inf(size(s.idx));
    end
end

% % check that all ub are non negative
% if any(s.ub < 0)
%     error('softlims_element_defaults: some soft limit for %s has a negative upper bound. There is most likely a problem in the specification of hl_val.', prop)
% end
% 
% %% rval: replacement value to remove constraint in initial OPF formulation
% if ismember(prop, {'ANGMAX', 'ANGMIN', 'RATE_A'})
%     s.rval = 0;
% elseif ismember(prop, {'VMAX', 'PMAX', 'QMAX'})
%     s.rval = Inf;
% elseif ismember(prop, {'QMIN', 'PMIN','VMIN'})
%     s.rval = -Inf;
% else
%     error('softlims_element_defaults: woops! property %s does not have an ''rval'' assigned', prop)
% end

%%-----  softlims_element_init  --------------------------------------------
function s = softlims_element_init(s, prop, mat, max_gen_cost)
%
%   s = softlims_element_init(s, prop, mat, max_gen_cost)
%
% for each property we want
%   idx: index of affected buses, branches, or generators
%   cost: linear cost to be added for the slack variable
%   hl_mod: 'remove'  : unbounded limit,
%           'replace' : replace hard limt,
%           'scale'   : hard limit is multiple of original limit
%           'shift'   : shift hard limit limit
%           'none'    : Don't apply softlimit to this property
%   hl_val: pairs with hl_mod to determine how the potentially new hard limit is set
%   ub: upper bound of slack variable. Calculated internaly
%   sav: original limits
%   rval: value to place in mpc structure to eliminate the constraint

%% ignore if hl_mod none specified
if isfield(s, 'hl_mod') && strcmp(s.hl_mod, 'none')
    return
end

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% default cost is default_cost * cost_scale.(prop)
default_cost  = 1000;
cost_scale = struct( ...
    'VMAX', 100, 'VMIN', 100, ...   %% $/pu
    'ANGMAX', 1, 'ANGMIN', 1, ...   %% $/deg
    'RATE_A', 1, ...                %% $/MW
    'PMIN', 1, 'PMAX', 1, ...       %% $/MW
    'QMIN', 1, 'QMAX', 1 );         %% $/MVar

if nargin < 4
    max_gen_cost = 0;
end

if ~ismember(prop, {'VMAX', 'VMIN', 'RATE_A', 'PMAX', 'PMIN', 'QMAX', 'QMIN', 'ANGMAX', 'ANGMIN'})
    error('softlims_element_init: Unknown limit type ''%s''', prop)
end

%% idx: indicies of bounds to relax
switch prop
    case {'VMAX', 'VMIN'}
        % Initial list of candidates for voltage limits contains all buses
        % NOTE: idxfull contains EXTERNAL bus numbers
        idxfull = mat(:, BUS_I);
    case {'ANGMAX', 'ANGMIN'}
        % Initial list of candidates for branch angle differences contains all
        % active branches with a meaningful limit (not 0, +360 or -360).
        % NOTE: idxfull contains locations in EXTERNAL branch list
        idxfull = find(mat(:, BR_STATUS) > 0 & mat(:, eval(prop)) & abs(mat(:, eval(prop))) < 360 );
    case 'RATE_A'
        % Initial list of candidates for line rating are all active branches
        % with rating not set to infinity (0).
        % NOTE: idxfull contains locations in EXTERNAL branch list
        idxfull = find(mat(:, BR_STATUS) > 0 & mat(:, RATE_A) > 0);
    case {'PMAX', 'QMAX', 'QMIN'}
        % Initial list of candidates for generation limits (except Pmin) are
        % all active generators whose limits are not infinity
        % NOTE: idxfull contains locations in EXTERNAL generator list
        idxfull = find( (mat(:, GEN_STATUS) > 0) & ~isload(mat) & ~isinf(mat(:, eval(prop))) );
    case 'PMIN'
        % Initial list of candidates for Pmin are all active generators that
        % are NOT dispachable loads
        % NOTE: idxfull contains locations in EXTERNAL generator list
        idxfull = find( (mat(:, GEN_STATUS) > 0) & ~isload(mat) );
end

% if isfield(s, 'idx')
    % idxmask is a boolean vector the size of s.idx. Entry i of
    % idxmask is 1 if s.idx(i) is in idxfull and 0 otherwise.
%     if size(s.idx, 2) > 1
%         s.idx = s.idx.';
%         if size(s.idx, 2) > 1
%             error('softlim_defaults: mpc.softlims.%s.idx must be a vector', prop)
%         end
%     end
    idxmask = ismember(s.idx, idxfull);
    s.idx = s.idx(idxmask); %remove possibly irrelevant entries entered by user
% else
%     % if no indices specified by user use the full list idxfull.
%     idxmask = true(size(idxfull));
%     s.idx = idxfull;
% end

if ismember(prop, {'VMAX', 'VMIN'})
    % for consistency between all the different limits, we want s.idx for
    % buses to contain locations rather than external bus numbers.
    % External bus numbers are stored in s.busnum and s.idx is rewritten to
    % include the locations (rows) of those buses.
    s.busnum = s.idx;
    s.idx =  find(ismember(mat(:,BUS_I), s.idx));
end

% if there are no constraints to relax, set hl_mod to none and exit
if isempty(s.idx)
    s.hl_mod = 'none';
    return
end

%% save original values in s.sav
s.sav = mat(s.idx, eval(prop));     %% s.idx is row index into relevant matrix

%% cost: cost of violating softlim
% if isfield(s, 'cost') && ~isempty(s.cost)
%     if isscalar(s.cost)
%         s.cost = s.cost * ones(size(s.idx));
%     else    % vector: apply idxmask filter in case elements removed in idx stage
        try
            s.cost = s.cost(idxmask);
% s_cost = s.cost
        catch ME
            warning('softlims_element_init: something went wrong when handling the cost for property %s. Perhaps the size of the ''cost'' vector didn''t match the size of the ''idx'' vector?', prop)
            rethrow(ME)
        end
%     end
% else
%     % default cost is the maximum between a predefined value and
%     % 2 times the maximum generation marginal cost. A scaling is applied
%     % to account for the fact that we are essentially equating 1MW with
%     % 0.01 p.u. voltage or 1 deg angle difference.
%     
%     ctmp = cost_scale.(prop) * max(default_cost, 2*max_gen_cost);
%     s.cost = ctmp * ones(size(s.idx));
% end

%% type of limit and upper bound of slack variable
ubsign = struct('VMAX', 1 , 'VMIN', -1, 'RATE_A', 1, ...
                'PMAX', 1, 'PMIN', -1, 'QMAX', 1, 'QMIN', -1, ...
                'ANGMAX', 1, 'ANGMIN', -1);
shift_defaults = struct('VMAX', 0.25 , 'VMIN', 0.25, 'RATE_A', 10, ...
                'PMAX', 10, 'PMIN', 10, 'QMAX', 10, 'QMIN', 10, ...
                'ANGMAX', 10, 'ANGMIN', 10);
if isfield(s, 'hl_mod')
    switch s.hl_mod
        case 'remove'
            % slack upper bound and new limit are both infinite
%             s.hl_val = ubsign.(prop)*Inf;
            s.ub = Inf(size(s.idx));
        case 'scale'
            % new hard limit is hl_val * original limit
            % for softlims in a positive direction the upper limit is:
            %           s.ub = (hl_val) * original limit - original limit
            %                = (hl_val - 1) * original limit
            % for softlims in a negative direction the upper limit is
            %           s.ub = original limit - hl_val*original limit
            %                = (1 - hl_val)*original limit
            switch vectorcheck(s, idxmask)
                case 0
                    error('softlims_element_init: provided hl_val vector does not conform to idx vector. When specifying hl_val as a vector idx must also be explicitly given') 
                case 1 % vector of values
                    s.ub = ubsign.(prop)*(s.hl_val(idxmask) - 1).*mat(s.idx, eval(prop));
                case 2 % scalar value
                    s.ub = ubsign.(prop)*(s.hl_val - 1)*mat(s.idx, eval(prop));
                case 3 % no hl_val specified, use default of 2 or 1/2
                    % use 2 if ubsign and original constraint have same sign
                    % otherwise use 1/2
%                     orig_lim = mat(s.idx, eval(prop));
%                     s.hl_val = 2 * ones(size(s.idx));   % scale up by 2
%                     k = find(ubsign.(prop) * orig_lim < 0); % unless opp. sign
%                     s.hl_val(k) = 1/2;                  % then scale down by 2
                    s.ub = abs((s.hl_val - 1) .* orig_lim);
            end
        case 'shift'
            % new hard limit is original limit + ubsign*hl_val
            % for positive direction limits:
            %    s.ub = original limit + hl_val - original limit = hl_val
            % for negative direction limits:
            %    s.ub = original limit - (original limit - hl_val) = hl_val
            % s.ub is therefore simply hl_val
            switch vectorcheck(s, idxmask)
                case 0
                    error('softlims_element_init: provided hl_val vector does not conform to idx vector. When specifying hl_val as a vector idx must also be explicitly given')
                case 1 % vector of values
                    s.ub = s.hl_val(idxmask);
                case 2 % scalar value
                    s.ub = s.hl_val*ones(size(s.idx));
                case 3 % no hl_val specified use defaults in shift_defaults
%                     s.hl_val = shift_defaults.(prop);
                    s.ub = s.hl_val*ones(size(s.idx));
            end
        case 'replace'
            % new hard limit is hl_val
            % for positive direction limits:
            %    s.ub = hl_val - original limit
            % for negative direction limits:
            %    s.ub = original limit - hl_val
            % Therefore s.ub = ubsign*(hl_val - original limit)
            switch vectorcheck(s, idxmask)
                case 0
                    error('softlims_element_init: provided hl_val vector does not conform to idx vector. When specifying hl_val as a vector idx must also be explicitly given')
                case 1 % vector of values
                    s.ub = ubsign.(prop)*(s.hl_val(idxmask) - mat(s.idx, eval(prop)));
                case 2 % scalar value
                    s.ub = ubsign.(prop)*(s.hl_val - mat(s.idx, eval(prop)));
                case 3 % non specified
                    error('softlims_element_init: for hard limit ''replace'' modification, replacement value hl_val must be specified')
            end
%         otherwise
%             error('softlims_element_init: unknown hard limit modification %s', s.hl_mod)
    end
else
    % defaults
    switch prop
        case 'PMIN'
%             % Normal generators (PMIN > 0) can only be relaxed to 0
%             % If PMIN < 0 allow PMIN to go to -Inf
%             s.hl_mod = 'replace';
%             s.hl_val = zeros(size(s.idx));
%             s.hl_val(mat(s.idx, PMIN) < 0) = -Inf;
            s.ub     = mat(s.idx, PMIN) - s.hl_val;
        case 'VMIN'
%             % VMIN can only be relaxed to zero, not further.
%             s.hl_mod = 'replace';
%             s.hl_val = 0;
            s.ub     = mat(s.idx, VMIN);
        case {'QMIN', 'ANGMIN'}
%             s.hl_mod = 'remove';
%             s.hl_val = -Inf;
            s.ub     = Inf(size(s.idx));
        otherwise
%             s.hl_mod = 'remove';
%             s.hl_val = Inf;
            s.ub     = Inf(size(s.idx));
    end
end

% check that all ub are non negative
if any(s.ub < 0)
    error('softlims_element_init: some soft limit for %s has a negative upper bound. There is most likely a problem in the specification of hl_val.', prop)
end

%% rval: replacement value to remove constraint in initial OPF formulation
if ismember(prop, {'ANGMAX', 'ANGMIN', 'RATE_A'})
    s.rval = 0;
elseif ismember(prop, {'VMAX', 'PMAX', 'QMAX'})
    s.rval = Inf;
elseif ismember(prop, {'QMIN', 'PMIN','VMIN'})
    s.rval = -Inf;
else
    error('softlims_element_init: woops! property %s does not have an ''rval'' assigned', prop)
end

%%----- vector check ---------------------------------------------------
function out = vectorcheck(s, idx)
% Utility function for softlims_element_defaults. Checks whether the hl_val
% field in the softlims element struct is a scalar or a vector whose size
% conforms with the idx vector
% OUTPUT:
%           out      0  hl_val is not scalar and does NOT conform to idx
%                    1  hl_val is not scalar and DOES conform to idx
%                    2  hl_val is a scalar
%                    3  hl_val is not given (use default)

if ~isfield(s, 'hl_val')
    out = 3;
elseif isscalar(s.hl_val)
    out = 2;
else
    out = all(size(s.hl_val) == size(idx));
end

%%-----  softlims_fcn  -------------------------------------------------
function [h, dh] = softlims_fcn(x, mpc, Yf, Yt, il, mpopt, Fmax)
%
%   Evaluates AC branch flow soft limit constraints and Jacobian.

%% options
lim_type = upper(mpopt.opf.flow_lim(1));

%% form constraints
%% compute flow (opf.flow_lim = 'P') or square of flow ('S', 'I', '2')
%% note that the Fmax used by opf_branch_flow_fcn() from branch matrix is zero
if nargout == 1
    h = opf_branch_flow_fcn(x(1:2), mpc, Yf, Yt, il, mpopt);
else
    [h, dh] = opf_branch_flow_fcn(x(1:2), mpc, Yf, Yt, il, mpopt);
end
flv = x{3};     %% flow limit violation variable
if lim_type == 'P'
    %%   Ff(Va,Vm) - flv <=  Fmax ===> Ff(Va,Vm) - flv - Fmax <= 0
    %%   Ft(Va,Vm) - flv <=  Fmax ===> Ff(Va,Vm) - flv - Fmax <= 0
    tmp2 = Fmax + flv;
else    %% lim_type == 'S', 'I', '2'
    %%   |Ff(Va,Vm)| - flv <=  Fmax ===> |Ff(Va,Vm)|.^2 - (flv + Fmax).^2 <= 0
    %%   |Ft(Va,Vm)| - flv <=  Fmax ===> |Ff(Va,Vm)|.^2 - (flv + Fmax).^2 <= 0
    tmp1 = Fmax + flv;
    tmp2 = tmp1.^2;
end
h = h - [tmp2; tmp2];
if nargout == 2
    ns = length(il);        %% number of soft limits
    if lim_type == 'P'
        tmp3 = -speye(ns, ns);
    else    %% lim_type == 'S', 'I', '2'
        tmp3 = spdiags(-2*tmp1, 0, ns, ns);
    end
    dh = [ dh [tmp3; tmp3] ];
end

%%-----  softlims_hess  ------------------------------------------------
function d2H = softlims_hess(x, lambda, mpc, Yf, Yt, il, mpopt)
%
%   Evaluates AC branch flow soft limit constraint Hessian.

%% options
lim_type = upper(mpopt.opf.flow_lim(1));

%% form Hessian
d2H = opf_branch_flow_hess(x(1:2), lambda, mpc, Yf, Yt, il, mpopt);
nh = size(d2H, 1);
ns = length(il);        %% number of soft limits

if lim_type == 'P'
    tmp = sparse(ns, ns);
else    %% lim_type == 'S', 'I', '2'
    tmp = spdiags(-2*lambda, 0, ns, ns);
end
d2H = [     d2H         sparse(nh, ns);
        sparse(ns, nh)      tmp         ];
