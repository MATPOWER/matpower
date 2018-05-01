function mpc = toggle_softlims(mpc, on_off, varargin)
%TOGGLE_SOFTLIMS Relax DC optimal power flow branch limits.
%   MPC = TOGGLE_SOFTLIMS(MPC, 'on')
%   MPC = TOGGLE_SOFTLIMS(MPC, 'off')
%   MPC = TOGGLE_SOFTLIMS(MPC, 'on', MPOPT)
%   T_F = TOGGLE_SOFTLIMS(MPC, 'status')
%
%   Enables, disables or checks the status of a set of OPF userfcn
%   callbacks to implement relaxed inequality constraints for an OPF model.
%
%   These callbacks expect to find a 'softlims' field in the input MPC,
%   where MPC.softlims is a struct with fields corresponding to the
%   possible limits, namely:
%       VMIN, VMAX, RATE_A, PMIN, PMAX, QMIN, QMAX, ANGMAX, ANGMIN,
%   Each of these is a structure in its own right with the following
%   fields:
%       idx     index of affected buses, branches, or generators. When
%               specifying buses, these should be bus numbers. For all
%               others these are indexes into the respective matrix. The
%               default are all elements that are not unbounded.
%
%       cost    linear cost to be added for the slack variable. Defaults
%               are:
%                   $100,000 $/pu   for VMAX and VMIN
%                      $1000 $/MW   for RATE_A, PMAX, and PMIN
%                      $1000 $/MVAr for QMAX, QMIN
%                      $1000 $/deg  for ANGMAX, ANGMIN
%
%       type    type of slack variable options are:
%                   'unbnd': unbounded limit
%                   'cnst' : constant upper bound for slack variable
%                   'frac' : multiplier of current limit
%                   'none' : No softlimit for this property
%
%       ub      slack variable upper bound. If a VECTOR is passed it is
%               assumed these are already the desired upperbounds (should
%               be done in conjunction with type 'cnst' or 'frac'). If a
%               SCALAR is passed in conjuction with 'cnst' the resulting ub
%               is: ub - abs(hard limit). The sign is flipped for VMIN and
%               PMIN, since they need to be positive quantities.
%               If a SCALAR is passed in conjunction with 'frac', the
%               resulting ub is: ub*abs(hard limit).
%
%       sav     original limits (this is handled in the defaults function).
%
%       rval    value to place in mpc structure to effectively eliminate
%               the constraint. Also handled in the default function.
%
%   Default values are assigned to any property that is not specified. The
%   default structure is:
%   softlims
%           .VMIN
%               .type = 'frac'
%               .ub   = 0.5
%           .VMAX
%               .type = 'frac'
%               .ub   = 0.5
%           .RATE_A
%               .type = 'frac'
%               .ub   = 0.5
%           .PMIN
%               .type = 'frac'
%               .ub   = 1
%           .PMAX
%               .type = 'unbnd'
%           .QMIN
%               .type = 'unbnd'
%           .QMAX
%               .type = 'unbnd'
%           .ANGMIN
%               .type = 'cnst'
%               .ub   = 360
%           .ANGMAX
%               .type = 'cnst'
%               .ub   = 360
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
%   MPC = TOGGLE_SOFTLIMS(MPC, 'on', MPOPT) additionally passes the
%   matpower options structure, which is added as an optional arg in the
%   ext2int callback. If mpopt.opf.softlims.default = 1 (default) the
%   default softlimits are applied to unspecified limits in the softlims
%   structure. If mpopt.opf.softlims.default = 0 the unspecified softlims
%   are ignored. In this case, to initialize the default soft limit a blank
%   structure should be created. For example, to use the default settings
%   for VMAX the command, mpc.softlims.VMAX = struct(), should be issued.
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
    if ~isempty(varargin)
        args = varargin{1};
    else
        args = [];
    end
    mpc = add_userfcn(mpc, 'ext2int', @userfcn_softlims_ext2int, args);
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

% structures used to index into softlims in a loop
lims = struct(...
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
mat2lims = struct(...
    'bus', {{'VMAX', 'VMIN'}}, ...
    'branch', {{'ANGMAX', 'ANGMIN', 'RATE_A'}}, ...
    'gen', {{'PMAX', 'PMIN', 'QMAX', 'QMIN'}} ...
);

%% check opf.softlims.default
use_default = 1;
if ~isempty(args)
    try 
        use_default = args.opf.softlims.default;
    catch
    end
end 
%% check for proper softlims inputs
if isfield(mpc, 'softlims')
    % loop over limit types
    for prop = fieldnames(lims).'
        mat  = lims.(prop{:});
        if isfield(mpc.softlims, prop{:})
            mpc.softlims.(prop{:}) = softlims_defaults(mpc.softlims.(prop{:}), prop{:}, mpc.order.ext.(mat));
        else
            if use_default
                % passing an empty struct, results in assignment of defaults
                mpc.softlims.(prop{:}) = softlims_defaults(struct(), prop{:}, mpc.order.ext.(mat));
            else
                mpc.softlims.(prop{:}).type = 'none';
            end
        end
    end
else
    % looping over limit types and assign defaults
    % passing an empty struct, results in assignment of defaults
    for prop = fieldnames(lims).'
        if use_default
            mat  = lims.(prop{:});
            mpc.softlims.(prop{:}) = softlims_defaults(struct(), prop{:}, mpc.order.ext.(mat));
        else
            mpc.softlims.(prop{:}).type = 'none';
        end
    end
end

%% initialize some things
s = mpc.softlims;
o = mpc.order;
%nl0 = size(o.ext.branch, 1);     %% original number of branches
%nl  = size(mpc.branch, 1);       %% number of on-line branches

%% save softlims struct for external indexing
mpc.order.ext.softlims = s;

%%-----  convert stuff to internal indexing  -----
for mat = {'bus', 'branch', 'gen'}
    mat = mat{1};
    n0  = size(o.ext.(mat), 1);  %% original number
    n   = size(mpc.(mat), 1);    %% on-line number
    e2i = zeros(n0, 1);
    e2i(o.(mat).status.on) = (1:n)';  %% ext->int index mapping
    for prop = mat2lims.(mat)
        if ~strcmp( s.(prop{:}).type, 'none')
            s.(prop{:}).idx = e2i(s.(prop{:}).idx);
            k = find(s.(prop{:}).idx == 0); %% find idxes corresponding to off-line elements
            s.(prop{:}).idx(k)     = [];    %% delete them
            s.(prop{:}).cost(k, :) = [];
            s.(prop{:}).sav(k)     = [];
            if ~isscalar(s.(prop{:}).ub)
                s.(prop{:}).ub(k)  = [];
            end
            if isempty(s.(prop{:}).idx)
                s.(prop{:}).type = 'none';
            end
        end
    end
end

% permute generators, since they are reordered inside of ext2int
for prop = mat2lims.gen
    if ~strcmp( s.(prop{:}).type, 'none')
        s.(prop{:}).idx   = s.(prop{:}).idx(o.gen.e2i);
    end
end


%%%%-------- remove hard limits on elements with soft limits
for prop = fieldnames(lims).'
    if ~strcmp(s.(prop{:}).type, 'none')
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
    if strcmp(s.(prop{:}).type, 'none')
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

if strcmp(mpopt.model, 'DC') && ~strcmp(s.RATE_A.type, 'none')
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
    if ~strcmp(s.VMIN.type, 'none')
        %%% Vm + s_vmin >= s.VMIN.sav
        ns  = length(s.VMIN.idx);
        Av  = sparse(1:ns, s.VMIN.idx, 1, ns, size(mpc.bus,1));
        As  = speye(ns);
        lb  = s.VMIN.sav;
        ub  = Inf(ns,1);

        om.add_lin_constraint('soft_vmin', [Av As], lb, ub, {'Vm', 's_vmin'});
    end
    if ~strcmp(s.VMAX.type, 'none')
        %%% Vm - s_vmax <= s.VMAX.sav
        ns  = length(s.VMAX.idx);
        Av  = sparse(1:ns, s.VMAX.idx, 1, ns, size(mpc.bus,1));
        As  = speye(ns);
        lb  = -Inf(ns,1);
        ub  = s.VMAX.sav;

        om.add_lin_constraint('soft_vmax', [Av -As], lb, ub, {'Vm', 's_vmax'});
    end
    if ~strcmp(s.QMIN.type, 'none')
        %%% Qg + s_pmin >= s.QMIN.max
        ns  = length(s.QMIN.idx);
        Av  = sparse(1:ns, s.QMIN.idx, 1, ns, size(mpc.gen,1));
        As  = speye(ns);
        lb  = s.QMIN.sav / mpc.baseMVA;
        ub  = Inf(ns,1);

        om.add_lin_constraint('soft_qmin', [Av As], lb, ub, {'Qg', 's_qmin'});
    end
    if ~strcmp(s.QMAX.type, 'none')
        %%% Qg - s_pmax <= s.QMAX.sav
        ns  = length(s.QMAX.idx);
        Av  = sparse(1:ns, s.QMAX.idx, 1, ns, size(mpc.gen,1));
        As  = speye(ns);
        lb  = -Inf(ns,1);
        ub  = s.QMAX.sav / mpc.baseMVA;

        om.add_lin_constraint('soft_qmax', [Av -As], lb, ub, {'Qg', 's_qmax'});
    end
    if ~strcmp(s.RATE_A.type, 'none')
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
lims = struct(...
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

%% get internal softlims struct and mpopt
isOPF  = isfield(results, 'f') && ~isempty(results.f);
if ~isOPF
    return
end
s = results.softlims;
mpopt = results.om.get_userdata('mpopt');   %% extract and remove mpopt from om
results.om.userdata = rmfield(results.om.userdata, 'mpopt');

%%-----  convert stuff back to external indexing  -----
o = results.order;
nl0 = size(o.ext.branch, 1);    %% original number of branches
nl = size(results.branch, 1);   %% number of on-line branches
o.branch.status.on;
results.softlims = results.order.ext.softlims;

%%-----  restore hard limits  -----
for prop = fieldnames(s).'
    mat  = lims.(prop{:});
    if strcmp(s.(prop{:}).type, 'none')
        continue
    end
    results.(mat)(s.(prop{:}).idx, eval(prop{:})) = s.(prop{:}).sav;
end

%%-----  results post-processing  -----
%% get overloads and overload costs
for prop = fieldnames(s).'
    mat  = lims.(prop{:});
    if strcmp(s.(prop{:}).type, 'none')
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
    results.softlims.(prop{:}).overload(o.(mat).status.on(s.(prop{:}).idx)) = var;
    results.softlims.(prop{:}).ovl_cost(o.(mat).status.on(s.(prop{:}).idx)) = var .* s.(prop{:}).cost(:,1);
end

%% get shadow prices
if ~strcmp(s.ANGMAX.type, 'none')
    results.branch(s.ANGMAX.idx, MU_ANGMAX) = results.lin.mu.u.soft_angmax * pi/180;
end
if ~strcmp(s.ANGMIN.type, 'none')
    results.branch(s.ANGMIN.idx, MU_ANGMIN) = results.lin.mu.l.soft_angmin * pi/180;
end
if ~strcmp(s.PMAX.type, 'none')
    results.gen(s.PMAX.idx, MU_PMAX) = results.lin.mu.u.soft_pmax / results.baseMVA;
end
if ~strcmp(s.PMIN.type, 'none')
    results.gen(s.PMIN.idx, MU_PMIN) = results.lin.mu.l.soft_pmin / results.baseMVA;
end
if strcmp(mpopt.model, 'DC')
    if ~strcmp(s.RATE_A.type, 'none')
        results.branch(s.RATE_A.idx, MU_SF) = results.lin.mu.u.softPf / results.baseMVA;
        results.branch(s.RATE_A.idx, MU_ST) = results.lin.mu.u.softPt / results.baseMVA;

        if 1    %% double-check value of overloads being returned
            vv = results.om.get_idx();
            check1 = zeros(nl0, 1);
            check1(o.branch.status.on(s.RATE_A.idx)) = results.x(vv.i1.s_rate_a:vv.iN.s_rate_a) * results.baseMVA;
            check2 = zeros(nl0, 1);
            k = find(results.branch(:, RATE_A) & ...
                     abs(results.branch(:, PF)) > results.branch(:, RATE_A) );
            check2(o.branch.status.on(k)) = ...
                    abs(results.branch(k, PF)) - results.branch(k, RATE_A);
            err1 = norm(results.softlims.RATE_A.overload-check1);
            err2 = norm(results.softlims.RATE_A.overload-check2);
            errtol = 1e-4;
            if err1 > errtol || err2 > errtol
                [ results.softlims.RATE_A.overload check1 results.softlims.RATE_A.overload-check1 ]
                [ results.softlims.RATE_A.overload check2 results.softlims.RATE_A.overload-check2 ]
                error('userfcn_softlims_int2ext: problem with consistency of overload values');
            end
        end
    end
else %AC model
    if ~strcmp(s.VMAX.type, 'none')
        results.bus(s.VMAX.idx, MU_VMAX) = results.lin.mu.u.soft_vmax;
    end
    if ~strcmp(s.VMIN.type, 'none')
        results.bus(s.VMIN.idx, MU_VMIN) = results.lin.mu.l.soft_vmin;
    end
    if ~strcmp(s.QMAX.type, 'none')
        results.gen(s.QMAX.idx, MU_QMAX) = results.lin.mu.u.soft_qmax / results.baseMVA;
    end
    if ~strcmp(s.QMIN.type, 'none')
        results.gen(s.QMIN.idx, MU_QMIN) = results.lin.mu.l.soft_qmin / results.baseMVA;
    end
    if ~strcmp(s.RATE_A.type, 'none')
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

        if 1    %% double-check value of overloads being returned
            vv = results.om.get_idx();
            check1 = zeros(nl0, 1);
            check1(o.branch.status.on(s.RATE_A.idx)) = results.x(vv.i1.s_rate_a:vv.iN.s_rate_a) * results.baseMVA;
            err1 = norm(results.softlims.RATE_A.overload-check1);
            errtol = 1e-4;
            if err1 > errtol
                [ results.softlims.RATE_A.overload check1 results.softlims.RATE_A.overload-check1 ]
                error('userfcn_softlims_int2ext: problem with consistency of overload values');
            end
        end
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
    if ~strcmp(s.RATE_A.type, 'none')
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
    if ~strcmp(s.VMAX.type,'none') && strcmp(mpopt.model, 'AC')
        k = find(s.VMAX.overload(s.VMAX.idx) | results.bus(s.VMAX.idx, MU_VMAX) > ptol);
        fprintf(fd, '\nMaximum Voltage Magnitude Limits:');
        fprintf(fd, '\n----------------------------------------');
        if isempty(k)
            fprintf(fd,'\nNo violations.\n');
        else
            fprintf(fd, '\nBus    Voltage   Limit   Overload    mu');
            fprintf(fd, '\n  #    Mag(pu)   (pu)     (pu)     ($/pu)');
            fprintf(fd, '\n-----  -------  -------  -------  ---------');
            fprintf(fd, '\n%5d%8.3f%9.3f%9.3f%11.3f',...
                [ s.VMAX.busidx(k), results.bus(s.VMAX.idx(k),[VM, VMAX]),...
                  s.VMAX.overload(s.VMAX.idx(k)), ...
                  results.bus(s.VMAX.idx(k), MU_VMAX)...
                ]');
            fprintf(fd, '\n                        --------');
            fprintf(fd, '\n               Total:%10.2f', ...
                    sum(s.VMAX.overload(s.VMAX.idx(k))));
            fprintf(fd, '\n');
        end
    end
    if ~strcmp(s.VMIN.type,'none') && strcmp(mpopt.model, 'AC')
        k = find(s.VMIN.overload(s.VMIN.idx) | results.bus(s.VMIN.idx, MU_VMIN) > ptol);
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
                [ s.VMIN.busidx(k), results.bus(s.VMIN.idx(k),[VM, VMIN]),...
                  s.VMIN.overload(s.VMIN.idx(k)), ...
                  results.bus(s.VMIN.idx(k), MU_VMIN)...
                ]');
            fprintf(fd, '\n                        --------');
            fprintf(fd, '\n               Total:%10.2f', ...
                    sum(s.VMIN.overload(s.VMIN.idx(k))));
            fprintf(fd, '\n');
        end
    end
    if ~strcmp(s.PMAX.type,'none')
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
    if ~strcmp(s.PMIN.type,'none')
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
    if ~strcmp(s.QMAX.type,'none') && strcmp(mpopt.model, 'AC')
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
    if ~strcmp(s.QMIN.type,'none') && strcmp(mpopt.model, 'AC')
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
    if ~strcmp(s.ANGMAX.type,'none')
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
    if ~strcmp(s.ANGMIN.type,'none')
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
lims = struct(...
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

% convenience structure for the different fields
fields = struct('type', struct('desc','Soft limit type', 'tok','%s'),...
    'busidx', struct('desc','bus ids for where soft voltage limit is applied','tok', '%d'),...
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
                    desc = sprintf(fields.idx.desc, lims.(prop{:}));
                else
                    desc = fields.(f).desc;
                end
                fprintf(fd,'%% %s\n',desc);
                if strcmp(f,'type')
                    fprintf(fd, '%ssoftlims.%s.type = ''%s'';\n', prefix,prop{:}, s.(prop{:}).type);
                else
                    fprintf(fd, '%ssoftlims.%s.%s = [\n', prefix,prop{:},f);
                    fprintf(fd, ['\t',fields.(f).tok,';\n'], s.(prop{:}).(f));
                    fprintf(fd, '];\n\n');
                end
            end
        end
    end
end

%%-----  softlims_defaults  --------------------------------------------
function s = softlims_defaults(s, prop, mat)
% for each property we want
%   idx: index of affected buses, branches, or generators
%	cost: linear cost to be added for the slack variable
%	type: 'unbnd': unbounded limit,
%         'cnst' : constant upper bound,
%		  'frac' : multiplier of current limit
%		  'none' : Don't apply softlimit to this property
%	ub: upper bound of slack variable, if 'cnst' it can be a scalar or a vector.
%   sav: original limits
%	rval: value to place in mpc structure to effectively eliminate the constraint
%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

default_cost = struct('VMAX', 100000, 'VMIN', 100000,...% $/pu
                      'ANGMAX', 1000, 'ANGMIN', 1000,...% $/deg
                      'RATE_A', 1000, ... %$/MW
                      'PMIN', 1000, 'PMAX', 1000, ... %$/MW
                      'QMIN', 1000, 'QMAX', 1000); %$/MW

%% ignore if type none specified
if isfield(s, 'type')
    if strcmp(s.type, 'none')
        return
    end
end

%% idx: indicies of bounds to relax
if ismember(prop, {'VMAX', 'VMIN'})
    % Initial list of candidates for voltage limits contains all buses
    % NOTE: idxfull contains EXTERNAL bus numbers
    
    idxfull = mat(:, BUS_I);
    
elseif ismember(prop, {'ANGMAX', 'ANGMIN'})
    % Initial list of candidates for branch angle differences contains all
    % active branches with a meaningful limit (not 0, +360 or -360).
    % NOTE: idxfull contains locations in EXTERNAL branch list
    
    idxfull = find(mat(:, BR_STATUS) > 0 & abs(mat(:, eval(prop))) < 360 & mat(:, eval(prop)) );
    
elseif strcmp(prop, 'RATE_A')
    % Initial list of candidates for line rating are all active branches
    % with rating not set to infinit (0).
    % NOTE: idxfull contains locations in EXTERNAL branch list
    
    idxfull = find(mat(:, BR_STATUS) > 0 & mat(:, RATE_A) > 0);
    
elseif ismember(prop, {'PMAX', 'QMAX', 'QMIN'})
    % Initial list of candidates for generation limits (except Pmin) are
    % all active generators
    % NOTE: idxfull contains locations in EXTERNAL generator list
    
    idxfull = find(mat(:, GEN_STATUS) > 0);
    
elseif strcmp(prop, 'PMIN')
    % Initial list of candidates for Pmin are all active generators with
    % non-zero Pmin.
    % NOTE: idxfull contains locations in EXTERNAL generator list
    
    idxfull = find(mat(:, GEN_STATUS) > 0 & mat(:, PMIN) > 0);
    
end

if isfield(s, 'idx')
    %idxmask returns a boolean vector the size of s.idx. Entry i of
    %idxmask is 1 if s.idx(i) is in idxfull and 0 otherwise.
    idxmask = ismember(s.idx, idxfull);
    s.idx = s.idx(idxmask); %remove possibly irrelevant entries entered by user
else
    % if no indices specified by user use the full list idxfull.
    idxmask = true(size(idxfull));
    s.idx = idxfull;
end

if ismember(prop, {'VMAX', 'VMIN'})
    % for consistency between all the different limits, we want s.idx for
    % buses to contain locations rather than external bus numbers.
    % External bus numbers are stored in s.busidx and s.idx is rewritten to
    % include the locations (rows) of those buses.
    s.busidx = s.idx;
    s.idx =  find(ismember(mat(:,BUS_I), s.idx));
end

% if there are no constraints to relax, set type to none and exit
if isempty(s.idx)
    s.type = 'none';
    return
end



%% sav: saves original values
% if ismember(prop, {'VMAX', 'VMIN'})
%     % save original bus numbers.
%     s.sav = mat(ismember(mat(:,BUS_I), s.idx), eval(prop));
if ismember(prop, {'VMAX', 'VMIN', 'RATE_A', 'PMAX', 'PMIN', 'QMAX', 'QMIN', 'ANGMAX', 'ANGMIN'})
    % Again note that here s.idx contains locations in the relevant matrix
    s.sav = mat(s.idx, eval(prop)) ;
else
    error('softlims_defaults: woops! property %s is not known.', prop)
end
%% cost: cost of violating softlim
if isfield(s, 'cost')
    if ~isscalar(s.cost)
        % if vector is specified aply the idxmask filter in case some
        % elements were removed in the idx stage.
        try
            s.cost = s.cost(idxmask);
        catch ME
            warning('softlims_defaults: something went wrong when handling the cost for property %s. Perhaps the size of the ''cost'' vector didn''t match the size of the ''idx'' vector?', prop)
            rethrow(ME)
        end
    else
        s.cost = s.cost * ones(size(s.idx));
    end
else
    s.cost = default_cost.(prop) * ones(size(s.idx));
end


%% type of limit and upper bound of slack variable
if isfield(s, 'type')
    if strcmp(s.type, 'unbnd')
        s.ub = Inf(size(s.idx));
    else
        if isfield(s, 'ub')
            if all(size(s.ub) == size(idxmask))
                % vector of upper bounds specified
                s.ub = s.ub(idxmask);
            elseif isscalar(s.ub)
                switch s.type
                    case 'cnst'
                        s.ub = s.ub - abs(mat(s.idx, eval(prop)));
                        if ismember(prop,{'VMIN','PMIN'})
                            % since VMIN and PMIN are expected to be
                            % positive qantities the sign needs to be
                            % flipped.
                            s.ub = -s.ub;
                        end
                    case 'frac'
                        s.ub = s.ub*abs(mat(s.idx, eval(prop)));
                    otherwise
                        error('softlims_defaults: unknown upper bound type %s', s.type)
                end
            else
                error('soflims_defaults: when specifying an upper bound ''ub'' must be the same shape as ''idx'' or a scalar.')
            end
        else
            error('sofltims_defaults: when specifying a non-unbounded type, field ''ub'' must be specified.')
        end
    end
else
    % defaults
    if ismember(prop, {'VMAX', 'VMIN'})
        s.type = 'frac';
        s.ub     = 0.5*mat(s.idx, eval(prop));
    elseif ismember(prop, {'ANGMAX', 'ANGMIN'})
        s.type = 'cnst';
        s.ub   = 360 - abs(mat(s.idx, eval(prop)));
    elseif strcmp(prop, 'RATE_A')
        s.type = 'frac';
        s.ub   = 0.5*mat(s.idx, eval(prop));
    elseif ismember(prop, {'PMAX', 'QMAX', 'QMIN'})
        s.type = 'unbnd';
        s.ub   = Inf(size(s.idx));
    elseif strcmp(prop, 'PMIN')
        s.type = 'frac';
        s.ub   = mat(s.idx, PMIN);
    end
end

% check that all ub are non negative
if any(s.ub < 0)
    error('softlims_defaults: negative upper bound found for softlim %s', prop)
end

%% rval: replacement value to remove constraint in initial OPF formulation
if ismember(prop, {'VMIN', 'PMIN', 'ANGMAX', 'ANGMIN', 'RATE_A'})
    s.rval = 0;
elseif ismember(prop, {'VMAX', 'PMAX', 'QMAX'})
    s.rval = Inf;
elseif strcmp(prop, 'QMIN')
    s.rval = -Inf;
else
    error('softlims_defaults: woops! property %s does not have an ''rval'' assigned', prop)
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
