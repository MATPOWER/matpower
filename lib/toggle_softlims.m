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

%% check for cartesian voltage coordinates
if ~isempty(mpopt) && isfield(mpopt, 'opf') && mpopt.opf.v_cartesian
    error('userfcn_softlims_ext2int: TOGGLE_SOFTLIMS is not implemented for OPF with cartesian voltages. Please set the MPOPT.opf.v_cartesian option to 0 for use with TOGGLE_SOFTLIMS.');
end

%% structures used to index into softlims in a loop
lims = softlims_lim2mat();
mat2lims = struct(...
    'bus', {{'VMAX', 'VMIN'}}, ...
    'branch', {{'ANGMAX', 'ANGMIN', 'RATE_A'}}, ...
    'gen', {{'PMAX', 'PMIN', 'QMAX', 'QMIN'}} ...
);

%% set up softlims defaults
if ~isfield(mpc, 'softlims')
    mpc.softlims = struct();
end
mpc.softlims = softlims_defaults(mpc, mpopt);

%% initialize some things
softlims = softlims_init(mpc, mpopt);
o = mpc.order;

%% save softlims struct with external indexing
mpc.order.ext.softlims = softlims;

%%-----  convert stuff to internal indexing  -----
for m = {'bus', 'branch', 'gen'}
    mat = m{:};
    n0  = size(o.ext.(mat), 1);     %% original number
    n   = size(mpc.(mat), 1);       %% on-line number
    e2i = zeros(n0, 1);             %% ext->int index mapping
    if strcmp(mat, 'gen')
        %% for gens, account for permutation of gen matrix
        e2i(o.gen.status.on(o.gen.i2e)) = (1:n)';
    else
        e2i(o.(mat).status.on) = (1:n)';
    end
    for p = mat2lims.(mat)
        prop = p{:};
        s = softlims.(prop);
        if ~strcmp( s.hl_mod, 'none')
            s.idx = e2i(s.idx);
            k = find(s.idx == 0);    %% find idxs corresponding to off-line elements
            s.idx(k)     = [];       %% delete them
            s.cost(k, :) = [];
            s.sav(k)     = [];
            s.ub(k)      = [];
            if isfield(s, 'hl_val') && ~isscalar(s.hl_val)
                s.hl_val(k)  = [];
            end
            softlims.(prop) = s;
        end
    end
end

%%-----  remove hard limits on elements with soft limits  -----
for p = fieldnames(lims).'
    prop = p{:};
    s = softlims.(prop);
    if ~strcmp(s.hl_mod, 'none')
        mat = lims.(prop);      %% mpc sub matrix
        mpc.(mat)(s.idx, eval(prop)) = s.rval;
    end
end

mpc.softlims = softlims;
mpc.order.int.softlims = softlims;


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
nb = size(bus, 1);
ng = size(mpc.gen, 1);

%% temporarily save mpopt in om for use by int2ext callback
om.userdata.mpopt = mpopt;

%%-----  add variables, costs, and constraints  -----
for p = fieldnames(mpc.softlims).'
    prop = p{:};
    s = mpc.softlims.(prop);
    if strcmp(s.hl_mod, 'none')
        continue;
    end
    varname = ['s_', lower(prop)];
    cstname = ['cs_', lower(prop)];
    ns = length(s.idx);  %% number of softlims

    %% variables and costs
    switch prop
        case {'VMIN', 'VMAX'}
            ub = s.ub;                      %% bound already in pu
            Cw = s.cost(:, 1);              %% cost already in $/pu
        case {'RATE_A', 'PMIN', 'PMAX', 'QMIN', 'QMAX'}
            ub = s.ub / baseMVA;            %% bound MVA -> pu
            Cw = s.cost(:, 1) * baseMVA;    %% cost in $/MVA -> $/pu
        case {'ANGMIN', 'ANGMAX'}
            ub = s.ub * pi/180;             %% bound deg -> rad
            Cw = s.cost(:, 1) * 180/pi;     %% cost in $/deg -> $/rad
        otherwise
            error('userfcn_soflims_formulation: property %s is unknown ', prop)
    end
    om.add_var(varname, ns, zeros(ns, 1), zeros(ns, 1), ub);
    om.add_quad_cost(cstname, [], Cw, 0, {varname});

    %% linear constraints that are identical for DC and AC formulations
    switch prop
        case 'ANGMIN'
            %% theta_f - theta_t + s_angmin >= s.sav
            ns = length(s.idx);
            Av = sparse([1:ns,1:ns].', ...
                    [branch(s.idx, F_BUS); branch(s.idx, T_BUS)], ...
                    [ones(ns,1);-ones(ns,1)], ns, nb);
            As = speye(ns);
            lb = s.sav * pi/180;
            ub = Inf(ns,1);

            om.add_lin_constraint('soft_angmin', [Av As], lb, ub, ...
                {'Va', 's_angmin'});
        case 'ANGMAX'
            %% theta_f - theta_t - s_angmax <= s.sav
            ns = length(s.idx);
            Av = sparse([1:ns,1:ns].', ...
                    [branch(s.idx, F_BUS); branch(s.idx, T_BUS)], ...
                    [ones(ns,1);-ones(ns,1)], ns, nb);
            As = speye(ns);
            lb = -Inf(ns,1);
            ub = s.sav * pi/180;

            om.add_lin_constraint('soft_angmax', [Av -As], lb, ub, ...
                {'Va', 's_angmax'});
        case 'PMIN'
            %% Pg + s_pmin >= s.sav
            ns = length(s.idx);
            Av = sparse(1:ns, s.idx, 1, ns, ng);
            As = speye(ns);
            lb = s.sav / baseMVA;
            ub = Inf(ns,1);

            om.add_lin_constraint('soft_pmin', [Av As], lb, ub, ...
                {'Pg', 's_pmin'});
        case 'PMAX'
            %% Pg - s_pmax <= s.sav
            ns = length(s.idx);
            Av = sparse(1:ns, s.idx, 1, ns, ng);
            As = speye(ns);
            lb = -Inf(ns,1);
            ub = s.sav / baseMVA;

            om.add_lin_constraint('soft_pmax', [Av -As], lb, ub, ...
                {'Pg', 's_pmax'});
    end
end

%% RATE_A (different for DC and AC) and other AC-only constraints
isDC = strcmp(mpopt.model, 'DC');
props = {'RATE_A'};
if ~isDC
    props = {'VMIN', 'VMAX', 'QMIN', 'QMAX', props{:}};
end
for p = props
    prop = p{:};
    s = mpc.softlims.(prop);
    if strcmp(s.hl_mod, 'none') %% skip
        continue;
    end
    ns = length(s.idx);         %% number of softlims

    switch prop
        case 'RATE_A'
            if isDC
                %% fetch Bf matrix for DC model
                Bf = om.get_userdata('Bf');
                Pfinj = om.get_userdata('Pfinj');

                %% form constraints
                %%    Bf * Va - s_rate_a <= -Pfinj + Pfmax
                %%   -Bf * Va - s_rate_a <=  Pfinj + Pfmax
                I = speye(ns, ns);
                Asf = [ Bf(s.idx, :) -I];
                Ast = [-Bf(s.idx, :) -I];
                lsf = -Inf(ns, 1);
                lst = lsf;
                usf = -Pfinj(s.idx) + s.sav/baseMVA;
                ust =  Pfinj(s.idx) + s.sav/baseMVA;
                vs = {'Va', 's_rate_a'};

                om.add_lin_constraint('softPf', Asf, lsf, usf, vs);    %% ns
                om.add_lin_constraint('softPt', Ast, lst, ust, vs);    %% ns
            else
                [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

                fcn = @(x)softlims_fcn(x, mpc, Yf(s.idx, :), Yt(s.idx, :), ...
                                       s.idx, mpopt, s.sav/baseMVA);
                hess = @(x, lam)softlims_hess(x, lam, mpc, Yf(s.idx, :), ...
                                              Yt(s.idx, :), s.idx, mpopt);
                om.add_nln_constraint({'softSf', 'softSt'}, [ns;ns], 0, ...
                                        fcn, hess, {'Va', 'Vm', 's_rate_a'});
            end
        case 'VMIN'
            %% Vm + s_vmin >= s.sav
            Av  = sparse(1:ns, s.idx, 1, ns, nb);
            As  = speye(ns);
            lb  = s.sav;
            ub  = Inf(ns, 1);

            om.add_lin_constraint('soft_vmin', [Av As], lb, ub, ...
                                    {'Vm', 's_vmin'});
        case 'VMAX'
            %% Vm - s_vmax <= s.sav
            Av  = sparse(1:ns, s.idx, 1, ns, nb);
            As  = speye(ns);
            lb  = -Inf(ns, 1);
            ub  = s.sav;

            om.add_lin_constraint('soft_vmax', [Av -As], lb, ub, ...
                                    {'Vm', 's_vmax'});
        case 'QMIN'
            %% Qg + s_pmin >= s.max
            Av  = sparse(1:ns, s.idx, 1, ns, ng);
            As  = speye(ns);
            lb  = s.sav / baseMVA;
            ub  = Inf(ns, 1);

            om.add_lin_constraint('soft_qmin', [Av As], lb, ub, ...
                                    {'Qg', 's_qmin'});
        case 'QMAX'
            %% Qg - s_pmax <= s.sav
            Av  = sparse(1:ns, s.idx, 1, ns, ng);
            As  = speye(ns);
            lb  = -Inf(ns, 1);
            ub  = s.sav / baseMVA;

            om.add_lin_constraint('soft_qmax', [Av -As], lb, ub, ...
                                    {'Qg', 's_qmax'});
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
%   It also expects the results to contain values in bus, gen and branch
%   and constraints for the specified softlims, with the following possible
%   names, which are used to populate output fields in results.softlims:
%       soft_Pf, softPt (DC only)
%       soft_Sf, softSt (AC only)
%       soft_angmin (both)
%       soft_angmax (both)
%       soft_pmin (both)
%       soft_pmax (both)
%       soft_vmin (AC only)
%       soft_vmax (AC only)
%       soft_qmin (AC only)
%       soft_qmax (AC only)
%   The optional args are not currently used.

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% structures used to index into softlims in a loop
lims = softlims_lim2mat();

%% get internal softlims struct and mpopt
softlims = results.softlims;
isOPF = isfield(results, 'f') && ~isempty(results.f);
if isOPF
    mpopt = results.om.get_userdata('mpopt');   %% extract and remove mpopt from om
    results.om.userdata = rmfield(results.om.userdata, 'mpopt');
    isDC = strcmp(mpopt.model, 'DC');
end

%%-----  convert stuff back to external indexing  -----
o = results.order;
results.softlims = o.ext.softlims;

%%-----  restore hard limits  -----
for p = fieldnames(softlims).'
    prop = p{:};
    s = softlims.(prop);
    if strcmp(s.hl_mod, 'none')
        continue;
    end
    mat = lims.(prop);
    results.(mat)(s.idx, eval(prop)) = s.sav;
end

%%-----  remove rval and sav fields  -----
for p = fieldnames(softlims).'
    prop = p{:};
    if strcmp(softlims.(prop).hl_mod, 'none')
        continue;
    end
    results.softlims.(prop) = rmfield(results.softlims.(prop), 'sav');
    results.softlims.(prop) = rmfield(results.softlims.(prop), 'rval');
end

%%-----  results post-processing  -----
%% get overloads and overload costs
tol = 1e-8;
if isOPF
    for p = fieldnames(softlims).'
        prop = p{:};
        s = softlims.(prop);
        if strcmp(s.hl_mod, 'none')
            continue;
        end

        mat  = lims.(prop);
        n0 = size(o.ext.(mat), 1);      %% original number
        n  = size(results.(mat), 1);    %% number on-line
        results.softlims.(prop).overload = zeros(n0, 1);
        results.softlims.(prop).ovl_cost = zeros(n0, 1);
        varname = ['s_', lower(prop)];  %% variable name
        if isfield(results.var.val, varname)
            var = results.var.val.(varname);    %% violation in p.u.
            var(var < tol) = 0;                 %% zero out violations < tol
        end
        switch prop
            case {'VMIN', 'VMAX'}               %% p.u. (no change)
                if isDC
                    continue;
                end
            case {'RATE_A', 'PMIN', 'PMAX'}     %% p.u. -> MVA|MW
                var = var * results.baseMVA;
            case {'QMIN', 'QMAX'}               %% p.u. -> MVAr
                if isDC
                    continue;
                end
                var = var * results.baseMVA;
            case {'ANGMIN', 'ANGMAX'}           %% rad -> deg
                var = var * 180/pi;
            otherwise
                error('userfcn_soflims_formulation: property %s is unknown ', prop)
        end

        %% NOTE: o.(mat).status.on is a vector nx1 where n is the INTERNAL number of
        %% elements. The entries are the EXTERNAL locations (row numbers).
        if strcmp(mat, 'gen')
            results.softlims.(prop).overload(o.gen.status.on(o.gen.i2e(s.idx))) = var;
            results.softlims.(prop).ovl_cost(o.gen.status.on(o.gen.i2e(s.idx))) = var .* s.cost(:, 1);
        else
            results.softlims.(prop).overload(o.(mat).status.on(s.idx)) = var;
            results.softlims.(prop).ovl_cost(o.(mat).status.on(s.idx)) = var .* s.cost(:, 1);
        end
    end

    %% get shadow prices
    if ~strcmp(softlims.ANGMAX.hl_mod, 'none')
        results.branch(softlims.ANGMAX.idx, MU_ANGMAX) = results.lin.mu.u.soft_angmax * pi/180;
    end
    if ~strcmp(softlims.ANGMIN.hl_mod, 'none')
        results.branch(softlims.ANGMIN.idx, MU_ANGMIN) = results.lin.mu.l.soft_angmin * pi/180;
    end
    if ~strcmp(softlims.PMAX.hl_mod, 'none')
        results.gen(softlims.PMAX.idx, MU_PMAX) = results.lin.mu.u.soft_pmax / results.baseMVA;
    end
    if ~strcmp(softlims.PMIN.hl_mod, 'none')
        results.gen(softlims.PMIN.idx, MU_PMIN) = results.lin.mu.l.soft_pmin / results.baseMVA;
    end
    if isDC
        if ~strcmp(softlims.RATE_A.hl_mod, 'none')
            results.branch(softlims.RATE_A.idx, MU_SF) = results.lin.mu.u.softPf / results.baseMVA;
            results.branch(softlims.RATE_A.idx, MU_ST) = results.lin.mu.u.softPt / results.baseMVA;
        end
    else        %% AC model
        if ~strcmp(softlims.VMAX.hl_mod, 'none')
            results.bus(softlims.VMAX.idx, MU_VMAX) = results.lin.mu.u.soft_vmax;
        end
        if ~strcmp(softlims.VMIN.hl_mod, 'none')
            results.bus(softlims.VMIN.idx, MU_VMIN) = results.lin.mu.l.soft_vmin;
        end
        if ~strcmp(softlims.QMAX.hl_mod, 'none')
            results.gen(softlims.QMAX.idx, MU_QMAX) = results.lin.mu.u.soft_qmax / results.baseMVA;
        end
        if ~strcmp(softlims.QMIN.hl_mod, 'none')
            results.gen(softlims.QMIN.idx, MU_QMIN) = results.lin.mu.l.soft_qmin / results.baseMVA;
        end
        if ~strcmp(softlims.RATE_A.hl_mod, 'none')
            if upper(mpopt.opf.flow_lim(1)) == 'P'
                results.branch(softlims.RATE_A.idx, MU_ST) = results.nli.mu.softSf / results.baseMVA;
                results.branch(softlims.RATE_A.idx, MU_SF) = results.nli.mu.softSt / results.baseMVA;
            else
                s = softlims.RATE_A;
                var = results.softlims.RATE_A.overload(o.branch.status.on(s.idx));
                %% conversion factor for squared constraints (2*F)
                cf = 2 * (s.sav + var) / results.baseMVA;
                results.branch(s.idx, MU_ST) = results.nli.mu.softSf .* cf / results.baseMVA;
                results.branch(s.idx, MU_SF) = results.nli.mu.softSt .* cf / results.baseMVA;
            end
        end
    end
end

%%-----  swap bus idx vectors  -----
%% restore softlims.VMIN.idx and softlims.VMAX.idx to external bus numbers
%% (stored in busnum). Store row indices in rowidx for use by printpf function.
for p = {'VMAX', 'VMIN'}
    prop = p{:};
    s = softlims.(prop);
    if ~strcmp(s.hl_mod, 'none')
        results.softlims.(prop).idx = s.busnum;
        results.softlims.(prop).rowidx = s.idx;
        results.softlims.(prop) = rmfield(results.softlims.(prop), 'busnum');
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
OUT_V_LIM       = OUT_ALL == 1 || (OUT_ALL == -1 && ~SUPPRESS && mpopt.out.lim.v);
OUT_LINE_LIM    = OUT_ALL == 1 || (OUT_ALL == -1 && ~SUPPRESS && mpopt.out.lim.line);
OUT_PG_LIM      = OUT_ALL == 1 || (OUT_ALL == -1 && ~SUPPRESS && mpopt.out.lim.pg);
OUT_QG_LIM      = OUT_ALL == 1 || (OUT_ALL == -1 && ~SUPPRESS && mpopt.out.lim.qg);

if isOPF && (results.success || OUT_FORCE)
    isAC = strcmp(mpopt.model, 'AC');
    sl = results.softlims;
    s = sl.RATE_A;
    if OUT_LINE_LIM && ~strcmp(s.hl_mod, 'none')
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n|     Soft Branch Flow Limits                                                  |');
        fprintf(fd, '\n================================================================================');
        k = find(s.overload(s.idx) | sum(results.branch(s.idx, MU_SF:MU_ST), 2) > ptol);
        if isempty(k)
            fprintf(fd,'\nNo violations.\n');
        else
            fprintf(fd, '\nBrnch   From   To     Flow      Limit    Overload     mu');
            fprintf(fd, '\n  #     Bus    Bus    (MW)      (MW)       (MW)     ($/MW)');
            fprintf(fd, '\n-----  -----  -----  --------  --------  --------  ---------');
            fprintf(fd, '\n%4d%7d%7d%10.2f%10.2f%10.2f%11.3f', ...
                    [   s.idx(k), results.branch(s.idx(k), [F_BUS, T_BUS]), ...
                        results.branch(s.idx(k), [PF, RATE_A]), ...
                        s.overload(s.idx(k)), ...
                        sum(results.branch(s.idx(k), MU_SF:MU_ST), 2) ...
                    ]');
            fprintf(fd, '\n                                         --------');
            fprintf(fd, '\n                                Total:%10.2f', ...
                    sum(s.overload(s.idx(k))));
            fprintf(fd, '\n');
        end
    end
    s = sl.VMAX;
    if isAC && OUT_V_LIM && ~strcmp(s.hl_mod,'none')
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n|     Soft Voltage Upper Bounds                                                |');
        fprintf(fd, '\n================================================================================');
        k = find(s.overload(s.rowidx) | results.bus(s.rowidx, MU_VMAX) > ptol);
        if isempty(k)
            fprintf(fd,'\nNo violations.\n');
        else
            fprintf(fd, '\nBus    Voltage   Limit   Overload    mu');
            fprintf(fd, '\n  #    Mag(pu)   (pu)     (pu)     ($/pu)');
            fprintf(fd, '\n-----  -------  -------  -------  ---------');
            fprintf(fd, '\n%5d%8.3f%9.3f%9.3f%11.3f',...
                [ s.idx(k), results.bus(s.rowidx(k),[VM, VMAX]),...
                  s.overload(s.rowidx(k)), ...
                  results.bus(s.rowidx(k), MU_VMAX)...
                ]');
            fprintf(fd, '\n                        --------');
            fprintf(fd, '\n               Total:%10.2f', ...
                    sum(s.overload(s.rowidx(k))));
            fprintf(fd, '\n');
        end
    end
    s = sl.VMIN;
    if isAC && OUT_V_LIM && ~strcmp(s.hl_mod,'none')
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n|     Soft Voltage Lower Bounds                                                |');
        fprintf(fd, '\n================================================================================');
        k = find(s.overload(s.rowidx) | results.bus(s.rowidx, MU_VMIN) > ptol);
        if isempty(k)
            fprintf(fd,'\nNo violations.\n');
        else
            fprintf(fd, '\n----------------------------------------');
            fprintf(fd, '\n Bus   Voltage   Limit   Overload    mu');
            fprintf(fd, '\n  #    Mag(pu)   (pu)     (pu)     ($/pu)');
            fprintf(fd, '\n-----  -------  -------  -------  ---------');
            fprintf(fd, '\n%5d%8.3f%9.3f%9.3f%11.3f',...
                [ s.idx(k), results.bus(s.rowidx(k),[VM, VMIN]),...
                  s.overload(s.rowidx(k)), ...
                  results.bus(s.rowidx(k), MU_VMIN)...
                ]');
            fprintf(fd, '\n                        --------');
            fprintf(fd, '\n               Total:%10.2f', ...
                    sum(s.overload(s.rowidx(k))));
            fprintf(fd, '\n');
        end
    end
    s = sl.PMAX;
    if OUT_PG_LIM && ~strcmp(s.hl_mod,'none')
        k = find(s.overload(s.idx) | results.gen(s.idx, MU_PMAX) > ptol);
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n|     Soft Generator Active Power Upper Bounds                                 |');
        fprintf(fd, '\n================================================================================');
        if isempty(k)
            fprintf(fd,'\nNo violations.\n');
        else
            fprintf(fd, '\nGen     Bus  Generation  Limit   Overload    mu');
            fprintf(fd, '\n  #      #     P (MW)    (MW)     (MW)     ($/MW)');
            fprintf(fd, '\n-----  -----  --------  -------  -------  ---------');
            fprintf(fd, '\n%5d%7d%9.2f%9.2f%9.3f%11.3f',...
                [ s.idx(k), results.gen(s.idx(k),[GEN_BUS, PG, PMAX]),...
                  s.overload(s.idx(k)), ...
                  results.gen(s.idx(k), MU_PMAX)...
                ]');
            fprintf(fd, '\n                                --------');
            fprintf(fd, '\n                       Total:%10.2f', ...
                    sum(s.overload(s.idx(k))));
            fprintf(fd, '\n');
        end
    end
    s = sl.PMIN;
    if OUT_PG_LIM && ~strcmp(s.hl_mod,'none')
        k = find(s.overload(s.idx) | results.gen(s.idx, MU_PMIN) > ptol);
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n|     Soft Generator Active Power Lower Bounds                                 |');
        fprintf(fd, '\n================================================================================');
        if isempty(k)
            fprintf(fd,'\nNo violations.\n');
        else
            fprintf(fd, '\nGen     Bus  Generation  Limit   Overload    mu');
            fprintf(fd, '\n  #      #     P (MW)    (MW)     (MW)     ($/MW)');
            fprintf(fd, '\n-----  -----  --------  -------  -------  ---------');
            fprintf(fd, '\n%5d%7d%9.2f%9.2f%9.3f%11.3f',...
                [ s.idx(k), results.gen(s.idx(k),[GEN_BUS, PG, PMIN]),...
                  s.overload(s.idx(k)), ...
                  results.gen(s.idx(k), MU_PMIN)...
                ]');
            fprintf(fd, '\n                                --------');
            fprintf(fd, '\n                       Total:%10.2f', ...
                    sum(s.overload(s.idx(k))));
            fprintf(fd, '\n');
        end
    end
    s = sl.QMAX;
    if isAC && OUT_QG_LIM && ~strcmp(s.hl_mod,'none')
        k = find(s.overload(s.idx) | results.gen(s.idx, MU_QMAX) > ptol);
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n|     Soft Generator Reactive Power Upper Bounds                               |');
        fprintf(fd, '\n================================================================================');
        if isempty(k)
            fprintf(fd,'\nNo violations.\n');
        else
            fprintf(fd, '\nGen     Bus  Generation  Limit   Overload    mu');
            fprintf(fd, '\n  #      #    Q (MVAr)  Q (MVAr)  (MVAr)   ($/MVAr)');
            fprintf(fd, '\n-----  -----  --------  -------  -------  ---------');
            fprintf(fd, '\n%5d%7d%9.2f%9.2f%9.3f%11.3f',...
                [ s.idx(k), results.gen(s.idx(k),[GEN_BUS, QG, QMAX]),...
                  s.overload(s.idx(k)), ...
                  results.gen(s.idx(k), MU_QMAX)...
                ]');
            fprintf(fd, '\n                                --------');
            fprintf(fd, '\n                       Total:%10.2f', ...
                    sum(s.overload(s.idx(k))));
            fprintf(fd, '\n');
        end
    end
    s = sl.QMIN;
    if isAC && OUT_QG_LIM && ~strcmp(s.hl_mod,'none')
        k = find(s.overload(s.idx) | results.gen(s.idx, MU_QMIN) > ptol);
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n|     Soft Generator Reactive Power Lower Bounds                               |');
        fprintf(fd, '\n================================================================================');
        if isempty(k)
            fprintf(fd,'\nNo violations.\n');
        else
            fprintf(fd, '\nGen     Bus  Generation  Limit   Overload    mu');
            fprintf(fd, '\n  #      #    Q (MVAr)  Q (MVAr)  (MVAr)   ($/MVAr)');
            fprintf(fd, '\n-----  -----  --------  -------  -------  ---------');
            fprintf(fd, '\n%5d%7d%9.2f%9.2f%9.3f%11.3f',...
                [ s.idx(k), results.gen(s.idx(k),[GEN_BUS, QG, QMIN]),...
                  s.overload(s.idx(k)), ...
                  results.gen(s.idx(k), MU_QMIN)...
                ]');
            fprintf(fd, '\n                                --------');
            fprintf(fd, '\n                       Total:%10.2f', ...
                    sum(s.overload(s.idx(k))));
            fprintf(fd, '\n');
        end
    end
    s = sl.ANGMAX;
    if OUT_V_LIM && ~strcmp(s.hl_mod,'none')
        k = find(s.overload(s.idx) | results.branch(s.idx, MU_ANGMAX) > ptol);
        delta = calc_branch_angle(results);
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n|     Soft Maximum Angle Difference Limits                                     |');
        fprintf(fd, '\n================================================================================');
        if isempty(k)
            fprintf(fd,'\nNo violations.\n');
        else
            fprintf(fd, '\nBrnch   From   To     Angle     Limit    Overload     mu');
            fprintf(fd, '\n  #     Bus    Bus    (deg)     (deg)     (deg)     ($/MW)');
            fprintf(fd, '\n-----  -----  -----  --------  --------  --------  ---------');
            fprintf(fd, '\n%4d%7d%7d%10.3f%10.3f%10.3f%11.3f', ...
                [ s.idx(k), results.branch(s.idx(k), [F_BUS, T_BUS]), ...
                  delta(s.idx(k)), results.branch(s.idx(k), ANGMAX),...
                  s.overload(s.idx(k)),...
                  results.branch(s.idx(k), MU_ANGMAX)...
                ]');
            fprintf(fd, '\n                                         --------');
            fprintf(fd, '\n                                Total:%10.2f', ...
                    sum(s.overload(s.idx(k))));
            fprintf(fd, '\n');
        end
    end
    s = sl.ANGMIN;
    if OUT_V_LIM && ~strcmp(s.hl_mod,'none')
        k = find(s.overload(s.idx) | results.branch(s.idx, MU_ANGMIN) > ptol);
        delta = calc_branch_angle(results);
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n|     Soft Minimum Angle Difference Limits                                     |');
        fprintf(fd, '\n================================================================================');
        if isempty(k)
            fprintf(fd,'\nNo violations.\n');
        else
            fprintf(fd, '\nBrnch   From   To     Angle     Limit    Overload     mu');
            fprintf(fd, '\n  #     Bus    Bus    (deg)     (deg)     (deg)     ($/MW)');
            fprintf(fd, '\n-----  -----  -----  --------  --------  --------  ---------');
            fprintf(fd, '\n%4d%7d%7d%10.3f%10.3f%10.3f%11.3f', ...
                [ s.idx(k), results.branch(s.idx(k), [F_BUS, T_BUS]), ...
                  delta(s.idx(k)), results.branch(s.idx(k), ANGMIN),...
                  s.overload(s.idx(k)),...
                  results.branch(s.idx(k), MU_ANGMIN)...
                ]');
            fprintf(fd, '\n                                         --------');
            fprintf(fd, '\n                                Total:%10.2f', ...
                    sum(s.overload(s.idx(k))));
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
function softlims = softlims_defaults(mpc, mpopt)
%
%   softlims = softlims_defaults(mpc, mpopt)
%
%   Returns a softlims field for mpc in which any missing inputs have
%   been filled in with defaults, so that each limit type has the all
%   4 input fields (idx, cost, hl_mod, hl_cost) and idx has been expanded
%   to a vector, unless hl_mod == 'none', in which case the others need
%   not be present (and are, in any case, ignored).

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% initialization
lims = softlims_lim2mat();
softlims = mpc.softlims;
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
for p = fieldnames(lims).'
    prop = p{:};
    mat  = mpc.order.ext.(lims.(prop));
    specified = isfield(softlims, prop) && ~isempty(softlims.(prop));
    if ~specified && ~use_default
        softlims.(prop).hl_mod = 'none';
    else
        if specified
            s = softlims.(prop);

            %% ignore if hl_mod = 'none'
            if isfield(s, 'hl_mod') && strcmp(s.hl_mod, 'none')
                continue;
            end
        else
            s = struct();
        end
%         softlims.(prop) = softlims_element_defaults(s, prop, mat, max_gen_cost);

        %% default base cost is the maximum between a predefined value and
        %% 2 times the maximum generation marginal cost. A scaling is applied
        %% to account for the fact that we are essentially equating 1MW with
        %% 0.01 p.u. voltage or 1 deg angle difference.
        default_base_cost = max(1000, 2 * max_gen_cost);
        cost_scale = struct( ...
            'VMAX', 100, 'VMIN', 100, ...   %% $/pu
            'ANGMAX', 1, 'ANGMIN', 1, ...   %% $/deg
            'RATE_A', 1, ...                %% $/MW
            'PMIN', 1, 'PMAX', 1, ...       %% $/MW
            'QMIN', 1, 'QMAX', 1 );         %% $/MVar
        default_cost = cost_scale.(prop) * default_base_cost;

        %% idx: indices of bounds to relax
        %% idxfull is the full list of candidate values for idx for given constraint
        %% type, all are row indices into EXTERNAL matrix, except for bus matrix,
        %% external bus numbers are used instead of row indices
        switch prop
            case {'VMAX', 'VMIN'}
                %% all buses
                idxfull = mat(:, BUS_I);
            case {'ANGMAX', 'ANGMIN'}
                %% all active branches with meaninful limit (not 0, +360 or -360)
                idxfull = find(mat(:, BR_STATUS) > 0 & mat(:, eval(prop)) & abs(mat(:, eval(prop))) < 360 );
            case 'RATE_A'
                %% all active branches with flow limit (RATE_A not 0)
                idxfull = find(mat(:, BR_STATUS) > 0 & mat(:, RATE_A) > 0);
            case {'PMAX', 'QMAX', 'QMIN'}
                %% all active generators (excluding dispatchable loads) w/finite limits
                idxfull = find( (mat(:, GEN_STATUS) > 0) & ~isload(mat) & ~isinf(mat(:, eval(prop))) );
            case 'PMIN'
                %% all active generators (excluding dispatchable loads)
                idxfull = find( (mat(:, GEN_STATUS) > 0) & ~isload(mat) );
        end

        %% check that idx is a scalar or vector
        if isfield(s, 'idx')
            if size(s.idx, 2) > 1
                s.idx = s.idx.';        %% make row vector into column vector
                if size(s.idx, 2) > 1
                    error('softlim_defaults: mpc.softlims.%s.idx must be a scalar or vector', prop);
                end
            end
        else    %% set default value for s.idx to full set of constraints
            s.idx = idxfull;
        end

        %% get row indexes (as opposed to bus numbers) for bus constraints
        if ismember(prop, {'VMAX', 'VMIN'})
            s_idx = find(ismember(mat(:,BUS_I), s.idx));
        else
            s_idx = s.idx;
        end

        %% if there are no constraints to relax, set hl_mod to 'none' and skip
        if isempty(s_idx)
            s.hl_mod = 'none';
            softlims.(prop) = s;
            continue;
        end

        %% cost: cost of violating softlim
        if isfield(s, 'cost') && ~isempty(s.cost)
            if isscalar(s.cost)     %% expand to vector if specified as scalar
                s.cost = s.cost * ones(size(s_idx));
            end
        else    %% not specified, use default cost
            s.cost = default_cost * ones(size(s_idx));
        end

        %% type of limit
        ubsign = struct('VMAX', 1 , 'VMIN', -1, 'RATE_A', 1, ...
                        'PMAX', 1, 'PMIN', -1, 'QMAX', 1, 'QMIN', -1, ...
                        'ANGMAX', 1, 'ANGMIN', -1);
        shift_defaults = struct('VMAX', 0.25 , 'VMIN', 0.25, 'RATE_A', 10, ...
                        'PMAX', 10, 'PMIN', 10, 'QMAX', 10, 'QMIN', 10, ...
                        'ANGMAX', 10, 'ANGMIN', 10);
        if isfield(s, 'hl_mod') %% use specified soft limits
            switch s.hl_mod
                case 'remove'   %% new hard limit is infinite (hl_val unused)
                case 'scale'    %% new hard limit is original limit * hl_val
                    if ~isfield(s, 'hl_val')
                        % use 2 if ubsign and original constraint have same sign
                        % otherwise use 1/2
                        orig_lim = mat(s_idx, eval(prop));
                        s.hl_val = 2 * ones(size(s_idx));   % scale up by 2
                        k = find(ubsign.(prop) * orig_lim < 0); % unless opp. sign
                        s.hl_val(k) = 1/2;                  % then scale down by 2
                    end
                case 'shift'    %% new hard limit is original limit + ubsign * hl_val
                    if ~isfield(s, 'hl_val')
                        s.hl_val = shift_defaults.(prop);
                    end
                case 'replace'  %% new hard limit is hl_val (no default value)
                otherwise
                    error('softlims_element_defaults: unknown hard limit modification %s', s.hl_mod);
            end
        else                    %% use default soft limits
            switch prop
                case 'PMIN'
                    %% for normal gens (PMIN > 0) PMIN is relaxed to 0, not removed
                    %% for gens with negative PMIN, it is relaxed to -Inf
                    s.hl_mod = 'replace';
                    s.hl_val = zeros(size(s_idx));
                    s.hl_val(mat(s_idx, PMIN) < 0) = -Inf;
                case 'VMIN'
                    %% VMIN is relaxed to zero, not removed
                    s.hl_mod = 'replace';
                    s.hl_val = 0;
                case {'QMIN', 'ANGMIN'}
                    s.hl_mod = 'remove';
                    s.hl_val = -Inf;
                otherwise
                    s.hl_mod = 'remove';
                    s.hl_val = Inf;
            end
        end

        softlims.(prop) = s;
    end
end


%%-----  softlims_init  --------------------------------------------
function softlims = softlims_init(mpc, mpopt)
%
%   softlims = softlims_init(mpc, mpopt)
%
%   Returns a softlims field for mpc in which any missing inputs have
%   been filled in with defaults, so that each limit type has the all
%   4 input fields (idx, cost, hl_mod, hl_cost) and idx has been expanded
%   to a vector, unless hl_mod == 'none', in which case the others need
%   not be present (and are, in any case, ignored).

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% initialization
lims = softlims_lim2mat();
softlims = mpc.softlims;

%% set defaults for each element
for p = fieldnames(lims).'
    prop = p{:};
    s = softlims.(prop);

    %% ignore if hl_mod = 'none'
    if isfield(s, 'hl_mod') && strcmp(s.hl_mod, 'none')
        continue;
    end

    mat  = mpc.order.ext.(lims.(prop));

%     softlims.(prop) = softlims_element_init(s, prop, mat);

    %% idx: indices of bounds to relax
    %% idxfull is the full list of candidate values for idx for given constraint
    %% type, all are row indices into EXTERNAL matrix, except for bus matrix,
    %% external bus numbers are used instead of row indices
    switch prop
        case {'VMAX', 'VMIN'}
            %% all buses
            idxfull = mat(:, BUS_I);
        case {'ANGMAX', 'ANGMIN'}
            %% all active branches with meaninful limit (not 0, +360 or -360)
            idxfull = find(mat(:, BR_STATUS) > 0 & mat(:, eval(prop)) & abs(mat(:, eval(prop))) < 360 );
        case 'RATE_A'
            %% all active branches with flow limit (RATE_A not 0)
            idxfull = find(mat(:, BR_STATUS) > 0 & mat(:, RATE_A) > 0);
        case {'PMAX', 'QMAX', 'QMIN'}
            %% all active generators (excluding dispatchable loads) w/finite limits
            idxfull = find( (mat(:, GEN_STATUS) > 0) & ~isload(mat) & ~isinf(mat(:, eval(prop))) );
        case 'PMIN'
            %% all active generators (excluding dispatchable loads)
            idxfull = find( (mat(:, GEN_STATUS) > 0) & ~isload(mat) );
    end

    %% idxmask is a boolean vector the size of s.idx, where
    %% idxmask(i) is 1 if s.idx(i) is in idxfull and 0 otherwise
    idxmask = ismember(s.idx, idxfull);
    s.idx = s.idx(idxmask); %% remove possibly irrelevant entries entered by user

    %% convert bus numbers to row indices for bus constraints
    if ismember(prop, {'VMAX', 'VMIN'})
        % External bus numbers are stored in s.busnum and s.idx is rewritten to
        % include the row indices of those buses.
        s.busnum = s.idx;
        s.idx =  find(ismember(mat(:,BUS_I), s.idx));
    end

    %% if there are no constraints to relax, skip
    if isempty(s.idx)
        softlims.(prop) = s;
        continue;
    end

    %% save original bounds in s.sav
    s.sav = mat(s.idx, eval(prop));

    %% cost: cost of violating softlim
    %% filter cost to include only valid constraints
    try
        s.cost = s.cost(idxmask);
    catch ME
        warning('softlims_init: something went wrong when handling the cost for property %s. Perhaps the size of the ''cost'' vector didn''t match the size of the ''idx'' vector?', prop);
        rethrow(ME);
    end

    %% upper bound on constraint violation variable (always a vector)
    %% ub = abs( new limit - original limit )
    if isfield(s, 'hl_mod')
        switch s.hl_mod
            case 'remove'   %% upper bound is infinite
                s.ub = Inf(size(s.idx));
            case 'scale'
                %% new limit = original limit * hl_val
                %% ub = abs( original limit * (hl_val - 1) )
                switch vectorcheck(s, idxmask)
                    case 0
                        error('softlims_init: provided hl_val vector does not conform to idx vector. When specifying hl_val as a vector idx must also be explicitly given');
                    case 1          %% vector
                        s.ub = abs( (s.hl_val(idxmask) - 1) .* s.sav );
                    case 2          %% scalar
                        s.ub = abs( (s.hl_val - 1)  * s.sav );
                    case {1, 3}     %% default, which is a vector
                        s.ub = abs( (s.hl_val - 1) .* s.sav );
                end
            case 'shift'
                %% new limit = original limit + ubsign * hl_val
                %% ub = hl_val
                switch vectorcheck(s, idxmask)
                    case 0
                        error('softlims_init: provided hl_val vector does not conform to idx vector. When specifying hl_val as a vector idx must also be explicitly given');
                    case 1          %% vector
                        s.ub = s.hl_val(idxmask);
                    case {2, 3}     %% scalar (including default which is a scalar)
                        s.ub = s.hl_val * ones(size(s.idx));
                end

                %% check that all ub are non negative
                if any(s.ub < 0)
                    error('softlims_init: some soft limit for %s has a negative upper bound. There is most likely a problem in the specification of hl_val.', prop);
                end
            case 'replace'
                %% new limit = hl_val
                %% ub = abs( hl_val - original limit )
                switch vectorcheck(s, idxmask)
                    case 0
                        error('softlims_init: provided hl_val vector does not conform to idx vector. When specifying hl_val as a vector idx must also be explicitly given');
                    case 1          %% vector
                        s.ub = abs( s.hl_val(idxmask) - s.sav );
                    case 2          %% scalar
                        s.ub = abs( s.hl_val - s.sav );
                    case 3          %% no default for 'replace'
                        error('softlims_init: for hard limit ''replace'' modification, replacement value hl_val must be specified');
                end
        end
    end

    %% rval: replacement value used to eliminate original hard constraint
    switch prop
        case {'ANGMAX', 'ANGMIN', 'RATE_A'}
            s.rval = 0;
        case {'VMAX', 'PMAX', 'QMAX'}
            s.rval = Inf;
        case {'QMIN', 'PMIN','VMIN'}
            s.rval = -Inf;
        otherwise
            error('softlims_init: property %s does not have an ''rval'' assigned', prop);
    end
    softlims.(prop) = s;
end


% %%-----  softlims_element_defaults  --------------------------------------------
% function s = softlims_element_defaults(s, prop, mat, max_gen_cost)
% %
% %   s = softlims_element_defaults(s, prop, mat, max_gen_cost)
% %
% %   For each property, we fill in default values as needed for:
% %       idx: index of affected buses, branches, or generators (bus numbers
% %           rather than indexes for bus)
% %       cost: linear cost on corresponding violation variables
% %       hl_mod: string value specifying what to do with original hard limits
% %               for this property
% %           'none'    : leave hard limits unchanged, do not apply soft limits
% %           'remove'  : eliminate hard limits completely
% %           'replace' : replace hard limit with value specified in 'hl_val'
% %           'scale'   : multiply hard limit by value specified in 'hl_val'
% %           'shift'   : add value specified in 'hl_val' to hard limit limit
% %       hl_val: value used to set/modify hard limit as specified by 'hl_mod'
% 
% %% ignore if hl_mod none specified
% if isfield(s, 'hl_mod') && strcmp(s.hl_mod, 'none')
%     return
% end
% 
% %% define named indices into data matrices
% [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
%     VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
% [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
%     MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
%     QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
% [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
%     TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
%     ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
% 
% if nargin < 4
%     max_gen_cost = 0;
% end
% 
% %% default base cost is the maximum between a predefined value and
% %% 2 times the maximum generation marginal cost. A scaling is applied
% %% to account for the fact that we are essentially equating 1MW with
% %% 0.01 p.u. voltage or 1 deg angle difference.
% default_base_cost = max(1000, 2 * max_gen_cost);
% cost_scale = struct( ...
%     'VMAX', 100, 'VMIN', 100, ...   %% $/pu
%     'ANGMAX', 1, 'ANGMIN', 1, ...   %% $/deg
%     'RATE_A', 1, ...                %% $/MW
%     'PMIN', 1, 'PMAX', 1, ...       %% $/MW
%     'QMIN', 1, 'QMAX', 1 );         %% $/MVar
% default_cost = cost_scale.(prop) * default_base_cost;
% 
% if ~ismember(prop, {'VMAX', 'VMIN', 'RATE_A', 'PMAX', 'PMIN', 'QMAX', 'QMIN', 'ANGMAX', 'ANGMIN'})
%     error('softlims_element_defaults: Unknown limit type ''%s''', prop);
% end
% 
% %% idx: indices of bounds to relax
% %% idxfull is the full list of candidate values for idx for given constraint
% %% type, all are row indices into EXTERNAL matrix, except for bus matrix,
% %% external bus numbers are used instead of row indices
% switch prop
%     case {'VMAX', 'VMIN'}
%         %% all buses
%         idxfull = mat(:, BUS_I);
%     case {'ANGMAX', 'ANGMIN'}
%         %% all active branches with meaninful limit (not 0, +360 or -360)
%         idxfull = find(mat(:, BR_STATUS) > 0 & mat(:, eval(prop)) & abs(mat(:, eval(prop))) < 360 );
%     case 'RATE_A'
%         %% all active branches with flow limit (RATE_A not 0)
%         idxfull = find(mat(:, BR_STATUS) > 0 & mat(:, RATE_A) > 0);
%     case {'PMAX', 'QMAX', 'QMIN'}
%         %% all active generators (excluding dispatchable loads) w/finite limits
%         idxfull = find( (mat(:, GEN_STATUS) > 0) & ~isload(mat) & ~isinf(mat(:, eval(prop))) );
%     case 'PMIN'
%         %% all active generators (excluding dispatchable loads)
%         idxfull = find( (mat(:, GEN_STATUS) > 0) & ~isload(mat) );
% end
% 
% %% check that idx is a scalar or vector
% if isfield(s, 'idx')
%     if size(s.idx, 2) > 1
%         s.idx = s.idx.';        %% make row vector into column vector
%         if size(s.idx, 2) > 1
%             error('softlim_defaults: mpc.softlims.%s.idx must be a scalar or vector', prop);
%         end
%     end
% else    %% set default value for s.idx to full set of constraints
%     s.idx = idxfull;
% end
% 
% %% get row indexes (as opposed to bus numbers) for bus constraints
% if ismember(prop, {'VMAX', 'VMIN'})
%     s_idx = find(ismember(mat(:,BUS_I), s.idx));
% else
%     s_idx = s.idx;
% end
% 
% %% if there are no constraints to relax, set hl_mod to none and exit
% if isempty(s_idx)
%     s.hl_mod = 'none';
%     return;
% end
% 
% %% cost: cost of violating softlim
% if isfield(s, 'cost') && ~isempty(s.cost)
%     if isscalar(s.cost)     %% expand to vector if specified as scalar
%         s.cost = s.cost * ones(size(s_idx));
%     end
% else    %% not specified, use default cost
%     s.cost = default_cost * ones(size(s_idx));
% end
% 
% %% type of limit
% ubsign = struct('VMAX', 1 , 'VMIN', -1, 'RATE_A', 1, ...
%                 'PMAX', 1, 'PMIN', -1, 'QMAX', 1, 'QMIN', -1, ...
%                 'ANGMAX', 1, 'ANGMIN', -1);
% shift_defaults = struct('VMAX', 0.25 , 'VMIN', 0.25, 'RATE_A', 10, ...
%                 'PMAX', 10, 'PMIN', 10, 'QMAX', 10, 'QMIN', 10, ...
%                 'ANGMAX', 10, 'ANGMIN', 10);
% if isfield(s, 'hl_mod') %% use specified soft limits
%     switch s.hl_mod
%         case 'remove'   %% new hard limit is infinite (hl_val unused)
%         case 'scale'    %% new hard limit is original limit * hl_val
%             if ~isfield(s, 'hl_val')
%                 % use 2 if ubsign and original constraint have same sign
%                 % otherwise use 1/2
%                 orig_lim = mat(s_idx, eval(prop));
%                 s.hl_val = 2 * ones(size(s_idx));   % scale up by 2
%                 k = find(ubsign.(prop) * orig_lim < 0); % unless opp. sign
%                 s.hl_val(k) = 1/2;                  % then scale down by 2
%             end
%         case 'shift'    %% new hard limit is original limit + ubsign * hl_val
%             if ~isfield(s, 'hl_val')
%                 s.hl_val = shift_defaults.(prop);
%             end
%         case 'replace'  %% new hard limit is hl_val (no default value)
%         otherwise
%             error('softlims_element_defaults: unknown hard limit modification %s', s.hl_mod);
%     end
% else                    %% use default soft limits
%     switch prop
%         case 'PMIN'
%             %% for normal gens (PMIN > 0) PMIN is relaxed to 0, not removed
%             %% for gens with negative PMIN, it is relaxed to -Inf
%             s.hl_mod = 'replace';
%             s.hl_val = zeros(size(s_idx));
%             s.hl_val(mat(s_idx, PMIN) < 0) = -Inf;
%         case 'VMIN'
%             %% VMIN is relaxed to zero, not removed
%             s.hl_mod = 'replace';
%             s.hl_val = 0;
%         case {'QMIN', 'ANGMIN'}
%             s.hl_mod = 'remove';
%             s.hl_val = -Inf;
%         otherwise
%             s.hl_mod = 'remove';
%             s.hl_val = Inf;
%     end
% end


% %%-----  softlims_element_init  --------------------------------------------
% function s = softlims_element_init(s, prop, mat)
% %
% %   s = softlims_element_init(s, prop, mat)
% %
% %   For each property, do some additional pre-processing, including:
% %       - Determine upper bound on softlim violation variables (ub)
% %       - Save values of original limits (sav)
% %       - Convert bus numbers (busnum) with bus row indices (idx)
% %       - Determine replacement value (rval) to use in mpc data matrices to
% %         eliminate the constraint
% 
% %% ignore if hl_mod none specified
% if isfield(s, 'hl_mod') && strcmp(s.hl_mod, 'none')
%     return
% end
% 
% %% define named indices into data matrices
% [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
%     VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
% [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
%     MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
%     QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
% [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
%     TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
%     ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
% 
% if ~ismember(prop, {'VMAX', 'VMIN', 'RATE_A', 'PMAX', 'PMIN', 'QMAX', 'QMIN', 'ANGMAX', 'ANGMIN'})
%     error('softlims_element_init: Unknown limit type ''%s''', prop);
% end
% 
% %% idx: indices of bounds to relax
% %% idxfull is the full list of candidate values for idx for given constraint
% %% type, all are row indices into EXTERNAL matrix, except for bus matrix,
% %% external bus numbers are used instead of row indices
% switch prop
%     case {'VMAX', 'VMIN'}
%         %% all buses
%         idxfull = mat(:, BUS_I);
%     case {'ANGMAX', 'ANGMIN'}
%         %% all active branches with meaninful limit (not 0, +360 or -360)
%         idxfull = find(mat(:, BR_STATUS) > 0 & mat(:, eval(prop)) & abs(mat(:, eval(prop))) < 360 );
%     case 'RATE_A'
%         %% all active branches with flow limit (RATE_A not 0)
%         idxfull = find(mat(:, BR_STATUS) > 0 & mat(:, RATE_A) > 0);
%     case {'PMAX', 'QMAX', 'QMIN'}
%         %% all active generators (excluding dispatchable loads) w/finite limits
%         idxfull = find( (mat(:, GEN_STATUS) > 0) & ~isload(mat) & ~isinf(mat(:, eval(prop))) );
%     case 'PMIN'
%         %% all active generators (excluding dispatchable loads)
%         idxfull = find( (mat(:, GEN_STATUS) > 0) & ~isload(mat) );
% end
% 
% %% idxmask is a boolean vector the size of s.idx, where
% %% idxmask(i) is 1 if s.idx(i) is in idxfull and 0 otherwise
% idxmask = ismember(s.idx, idxfull);
% % if ~all(size(s.idx) == size(s.idx(idxmask))) || any(s.idx ~= s.idx(idxmask))
% %     s_idx = s.idx
% %     s_idx_masked = s.idx(idxmask)
% % end
% s.idx = s.idx(idxmask); %remove possibly irrelevant entries entered by user
% 
% %% convert bus numbers to row indices for bus constraints
% if ismember(prop, {'VMAX', 'VMIN'})
%     % External bus numbers are stored in s.busnum and s.idx is rewritten to
%     % include the row indices of those buses.
%     s.busnum = s.idx;
%     s.idx =  find(ismember(mat(:,BUS_I), s.idx));
% end
% 
% %% if there are no constraints to relax, set hl_mod to none and exit
% if isempty(s.idx)
%     return;
% end
% 
% %% save original bounds in s.sav
% s.sav = mat(s.idx, eval(prop));
% 
% %% cost: cost of violating softlim
% %% filter cost to include only valid constraints
% try
%     s.cost = s.cost(idxmask);
% catch ME
%     warning('softlims_element_init: something went wrong when handling the cost for property %s. Perhaps the size of the ''cost'' vector didn''t match the size of the ''idx'' vector?', prop);
%     rethrow(ME);
% end
% 
% %% upper bound on constraint violation variable
% %% ub = abs( new limit - original limit )
% if isfield(s, 'hl_mod')
%     switch s.hl_mod
%         case 'remove'   %% upper bound is infinite
%             s.ub = Inf(size(s.idx));
%         case 'scale'
%             %% new limit = original limit * hl_val
%             %% ub = abs( original limit * (hl_val - 1) )
%             switch vectorcheck(s, idxmask)
%                 case 0
%                     error('softlims_element_init: provided hl_val vector does not conform to idx vector. When specifying hl_val as a vector idx must also be explicitly given');
%                 case 1          %% vector
%                     s.ub = abs( (s.hl_val(idxmask) - 1) .* s.sav );
%                 case 2          %% scalar
%                     s.ub = abs( (s.hl_val - 1)  * s.sav );
%                 case {1, 3}     %% default, which is a vector
%                     s.ub = abs( (s.hl_val - 1) .* s.sav );
%             end
%         case 'shift'
%             %% new limit = original limit + ubsign * hl_val
%             %% ub = hl_val
%             switch vectorcheck(s, idxmask)
%                 case 0
%                     error('softlims_element_init: provided hl_val vector does not conform to idx vector. When specifying hl_val as a vector idx must also be explicitly given');
%                 case 1          %% vector
%                     s.ub = s.hl_val(idxmask);
%                 case {2, 3}     %% scalar (including default which is a scalar)
%                     s.ub = s.hl_val * ones(size(s.idx));
%             end
% 
%             %% check that all ub are non negative
%             if any(s.ub < 0)
%                 error('softlims_element_init: some soft limit for %s has a negative upper bound. There is most likely a problem in the specification of hl_val.', prop);
%             end
%         case 'replace'
%             %% new limit = hl_val
%             %% ub = abs( hl_val - original limit )
%             switch vectorcheck(s, idxmask)
%                 case 0
%                     error('softlims_element_init: provided hl_val vector does not conform to idx vector. When specifying hl_val as a vector idx must also be explicitly given');
%                 case 1          %% vector
%                     s.ub = abs( s.hl_val(idxmask) - s.sav );
%                 case 2          %% scalar
%                     s.ub = abs( s.hl_val - s.sav );
%                 case 3          %% no default for 'replace'
%                     error('softlims_element_init: for hard limit ''replace'' modification, replacement value hl_val must be specified')
%             end
%     end
% end
% 
% %% rval: replacement value used to eliminate original hard constraint
% switch prop
%     case {'ANGMAX', 'ANGMIN', 'RATE_A'}
%         s.rval = 0;
%     case {'VMAX', 'PMAX', 'QMAX'}
%         s.rval = Inf;
%     case {'QMIN', 'PMIN','VMIN'}
%         s.rval = -Inf;
%     otherwise
%         error('softlims_element_init: property %s does not have an ''rval'' assigned', prop);
% end

%%----- vector check ---------------------------------------------------
function out = vectorcheck(s, idx)
% Utility function for softlims_init(). Checks whether the hl_val
% field in the softlims element struct is a scalar or a vector whose size
% conforms with the idx vector (prior to filtering with idxmask)
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
