function mpc = toggle_reserves(mpc, on_off)
%TOGGLE_RESERVES Enable, disable or check status of fixed reserve requirements.
%   MPC = TOGGLE_RESERVES(MPC, 'on')
%   MPC = TOGGLE_RESERVES(MPC, 'off')
%   T_F = TOGGLE_RESERVES(MPC, 'status')
%
%   Enables, disables or checks the status of a set of OPF userfcn
%   callbacks to implement co-optimization of reserves with fixed zonal
%   reserve requirements.
%
%   These callbacks expect to find a 'reserves' field in the input MPC,
%   where MPC.reserves is a struct with the following fields:
%       zones   nrz x ng, zone(i, j) = 1, if gen j belongs to zone i
%                                      0, otherwise
%       req     nrz x 1, zonal reserve requirement in MW
%       cost    (ng or ngr) x 1, cost of reserves in $/MW
%       qty     (ng or ngr) x 1, max quantity of reserves in MW (optional)
%   where nrz is the number of reserve zones and ngr is the number of
%   generators belonging to at least one reserve zone and ng is the total
%   number of generators.
%
%   The 'int2ext' callback also packages up results and stores them in
%   the following output fields of results.reserves:
%       R       - ng x 1, reserves provided by each gen in MW
%       Rmin    - ng x 1, lower limit on reserves provided by each gen, (MW)
%       Rmax    - ng x 1, upper limit on reserves provided by each gen, (MW)
%       mu.l    - ng x 1, shadow price on reserve lower limit, ($/MW)
%       mu.u    - ng x 1, shadow price on reserve upper limit, ($/MW)
%       mu.Pmax - ng x 1, shadow price on Pg + R <= Pmax constraint, ($/MW)
%       prc     - ng x 1, reserve price for each gen equal to maximum of the
%                         shadow prices on the zonal requirement constraint
%                         for each zone the generator belongs to
%
%   See also RUNOPF_W_RES, ADD_USERFCN, REMOVE_USERFCN, RUN_USERFCN,
%   T_CASE30_USERFCNS.

%   MATPOWER
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if strcmp(upper(on_off), 'ON')
    %% check for proper reserve inputs
    if ~isfield(mpc, 'reserves') || ~isstruct(mpc.reserves) || ...
            ~isfield(mpc.reserves, 'zones') || ...
            ~isfield(mpc.reserves, 'req') || ...
            ~isfield(mpc.reserves, 'cost')
        error('toggle_reserves: case must contain a ''reserves'' field, a struct defining ''zones'', ''req'' and ''cost''');
    end
    
    %% add callback functions
    %% note: assumes all necessary data included in 1st arg (mpc, om, results)
    %%       so, no additional explicit args are needed
    mpc = add_userfcn(mpc, 'ext2int', @userfcn_reserves_ext2int);
    mpc = add_userfcn(mpc, 'formulation', @userfcn_reserves_formulation);
    mpc = add_userfcn(mpc, 'int2ext', @userfcn_reserves_int2ext);
    mpc = add_userfcn(mpc, 'printpf', @userfcn_reserves_printpf);
    mpc = add_userfcn(mpc, 'savecase', @userfcn_reserves_savecase);
    mpc.userfcn.status.reserves = 1;
elseif strcmp(upper(on_off), 'OFF')
    mpc = remove_userfcn(mpc, 'savecase', @userfcn_reserves_savecase);
    mpc = remove_userfcn(mpc, 'printpf', @userfcn_reserves_printpf);
    mpc = remove_userfcn(mpc, 'int2ext', @userfcn_reserves_int2ext);
    mpc = remove_userfcn(mpc, 'formulation', @userfcn_reserves_formulation);
    mpc = remove_userfcn(mpc, 'ext2int', @userfcn_reserves_ext2int);
    mpc.userfcn.status.reserves = 0;
elseif strcmp(upper(on_off), 'STATUS')
    if isfield(mpc, 'userfcn') && isfield(mpc.userfcn, 'status') && ...
            isfield(mpc.userfcn.status, 'reserves')
        mpc = mpc.userfcn.status.reserves;
    else
        mpc = 0;
    end
else
    error('toggle_reserves: 2nd argument must be ''on'', ''off'' or ''status''');
end


%%-----  ext2int  ------------------------------------------------------
function mpc = userfcn_reserves_ext2int(mpc, mpopt, args)
%
%   mpc = userfcn_reserves_ext2int(mpc, mpopt, args)
%
%   This is the 'ext2int' stage userfcn callback that prepares the input
%   data for the formulation stage. It expects to find a 'reserves' field
%   in mpc as described above. The optional args are not currently used.

%% initialize some things
r = mpc.reserves;
o = mpc.order;
ng0 = size(o.ext.gen, 1);   %% number of original gens (+ disp loads)
nrz = size(r.req, 1);       %% number of reserve zones
if nrz > 1
    mpc.reserves.rgens = any(r.zones);  %% mask of gens available to provide reserves
else
    mpc.reserves.rgens = r.zones;
end
igr = find(mpc.reserves.rgens); %% indices of gens available to provide reserves
ngr = length(igr);              %% number of gens available to provide reserves

%% check data for consistent dimensions
if size(r.zones, 1) ~= nrz
    error('userfcn_reserves_ext2int: the number of rows in mpc.reserves.req (%d) and mpc.reserves.zones (%d) must match', nrz, size(r.zones, 1));
end
if size(r.cost, 1) ~= ng0 && size(r.cost, 1) ~= ngr
    error('userfcn_reserves_ext2int: the number of rows in mpc.reserves.cost (%d) must equal the total number of generators (%d) or the number of generators able to provide reserves (%d)', size(r.cost, 1), ng0, ngr);
end
if isfield(r, 'qty') && size(r.qty, 1) ~= size(r.cost, 1)
    error('userfcn_reserves_ext2int: mpc.reserves.cost (%d x 1) and mpc.reserves.qty (%d x 1) must be the same dimension', size(r.cost, 1), size(r.qty, 1));
end

%% convert both cost and qty from ngr x 1 to full ng x 1 vectors if necessary
if size(r.cost, 1) < ng0
    mpc.reserves.original.cost = r.cost;    %% save original
    cost = zeros(ng0, 1);
    cost(igr) = r.cost;
    mpc.reserves.cost = cost;
    if isfield(r, 'qty')
        mpc.reserves.original.qty = r.qty;  %% save original
        qty = zeros(ng0, 1);
        qty(igr) = r.qty;
        mpc.reserves.qty = qty;
    end
end

%%-----  convert stuff to internal indexing  -----
%% convert all reserve parameters (zones, costs, qty, rgens)
if isfield(r, 'qty')
    mpc = e2i_field(mpc, {'reserves', 'qty'}, 'gen');
end
mpc = e2i_field(mpc, {'reserves', 'cost'}, 'gen');
mpc = e2i_field(mpc, {'reserves', 'zones'}, 'gen', 2);
mpc = e2i_field(mpc, {'reserves', 'rgens'}, 'gen', 2);

%% save indices of gens available to provide reserves
mpc.order.ext.reserves.igr = igr;               %% external indexing
mpc.reserves.igr = find(mpc.reserves.rgens);    %% internal indexing


%%-----  formulation  --------------------------------------------------
function om = userfcn_reserves_formulation(om, mpopt, args)
%
%   om = userfcn_reserves_formulation(om, mpopt, args)
%
%   This is the 'formulation' stage userfcn callback that defines the
%   user costs and constraints for fixed reserves. It expects to find
%   a 'reserves' field in the mpc stored in om, as described above.
%   By the time it is passed to this callback, mpc.reserves should
%   have two additional fields:
%       igr     1 x ngr, indices of generators available for reserves
%       rgens   1 x ng, 1 if gen avaiable for reserves, 0 otherwise
%   It is also assumed that if cost or qty were ngr x 1, they have been
%   expanded to ng x 1 and that everything has been converted to
%   internal indexing, i.e. all gens are on-line (by the 'ext2int'
%   callback). The optional args are not currently used.

%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% initialize some things
mpc = om.get_mpc();
r = mpc.reserves;
igr = r.igr;                %% indices of gens available to provide reserves
ngr = length(igr);          %% number of gens available to provide reserves
ng  = size(mpc.gen, 1);     %% number of on-line gens (+ disp loads)

%% variable bounds
Rmin = zeros(ngr, 1);               %% bound below by 0
Rmax = Inf(ngr, 1);                 %% bound above by ...
k = find(mpc.gen(igr, RAMP_10));
Rmax(k) = mpc.gen(igr(k), RAMP_10); %% ... ramp rate and ...
if isfield(r, 'qty')
    k = find(r.qty(igr) < Rmax);
    Rmax(k) = r.qty(igr(k));        %% ... stated max reserve qty
end
Rmax = Rmax / mpc.baseMVA;

%% constraints
I = speye(ngr);                     %% identity matrix
Ar = [sparse(1:ngr, igr, 1, ngr, ng) I];
ur = mpc.gen(igr, PMAX) / mpc.baseMVA;
lreq = r.req / mpc.baseMVA;

%% cost
Cw = r.cost(igr) * mpc.baseMVA;     %% per unit cost coefficients

%% add them to the model
om.add_var('R', ngr, [], Rmin, Rmax);
om.add_lin_constraint('Pg_plus_R', Ar, [], ur, {'Pg', 'R'});
om.add_lin_constraint('Rreq', r.zones(:, igr), lreq, [], {'R'});
om.add_quad_cost('Rcost', [], Cw, 0, {'R'});


%%-----  int2ext  ------------------------------------------------------
function results = userfcn_reserves_int2ext(results, mpopt, args)
%
%   results = userfcn_reserves_int2ext(results, mpopt, args)
%
%   This is the 'int2ext' stage userfcn callback that converts everything
%   back to external indexing and packages up the results. It expects to
%   find a 'reserves' field in the results struct as described for mpc
%   above, including the two additional fields 'igr' and 'rgens'. It also
%   expects the results to contain a variable 'R' and linear constraints
%   'Pg_plus_R' and 'Rreq' which are used to populate output fields in
%   results.reserves. The optional args are not currently used.

%% initialize some things
r = results.reserves;

%% grab some info in internal indexing order
igr = r.igr;                %% indices of gens available to provide reserves
ng  = size(results.gen, 1); %% number of on-line gens (+ disp loads)

%%-----  convert stuff back to external indexing  -----
%% convert all reserve parameters (zones, costs, qty, rgens)
if isfield(r, 'qty')
    results = i2e_field(results, {'reserves', 'qty'}, 'gen');
end
results = i2e_field(results, {'reserves', 'cost'}, 'gen');
results = i2e_field(results, {'reserves', 'zones'}, 'gen', 2);
results = i2e_field(results, {'reserves', 'rgens'}, 'gen', 2);
results.order.int.reserves.igr = results.reserves.igr;  %% save internal version
results.reserves.igr = results.order.ext.reserves.igr;  %% use external version
r = results.reserves;       %% update
o = results.order;          %% update

%% grab same info in external indexing order
igr0 = r.igr;               %% indices of gens available to provide reserves
ng0  = size(o.ext.gen, 1);  %% number of gens (+ disp loads)

%%-----  results post-processing  -----
%% get the results (per gen reserves, multipliers) with internal gen indexing
%% and convert from p.u. to per MW units
[R0, Rl, Ru] = results.om.params_var('R');
R       = zeros(ng, 1);
Rmin    = zeros(ng, 1);
Rmax    = zeros(ng, 1);
mu_l    = zeros(ng, 1);
mu_u    = zeros(ng, 1);
mu_Pmax = zeros(ng, 1);
R(igr)       = results.var.val.R * results.baseMVA;
Rmin(igr)    = Rl * results.baseMVA;
Rmax(igr)    = Ru * results.baseMVA;
mu_l(igr)    = results.var.mu.l.R / results.baseMVA;
mu_u(igr)    = results.var.mu.u.R / results.baseMVA;
mu_Pmax(igr) = results.lin.mu.u.Pg_plus_R / results.baseMVA;

%% store in results in results struct
z = zeros(ng0, 1);
results.reserves.R       = i2e_data(results, R, z, 'gen');
results.reserves.Rmin    = i2e_data(results, Rmin, z, 'gen');
results.reserves.Rmax    = i2e_data(results, Rmax, z, 'gen');
results.reserves.mu.l    = i2e_data(results, mu_l, z, 'gen');
results.reserves.mu.u    = i2e_data(results, mu_u, z, 'gen');
results.reserves.mu.Pmax = i2e_data(results, mu_Pmax, z, 'gen');
results.reserves.prc     = z;
for k = igr0
    iz = find(r.zones(:, k));
    results.reserves.prc(k) = sum(results.lin.mu.l.Rreq(iz)) / results.baseMVA;
end
results.reserves.totalcost = sum(results.qdc.Rcost);

%% replace ng x 1 cost, qty with ngr x 1 originals
if isfield(r, 'original')
    if isfield(r, 'qty')
        results.reserves.qty = r.original.qty;
    end
    results.reserves.cost = r.original.cost;
    results.reserves = rmfield(results.reserves, 'original');
end


%%-----  printpf  ------------------------------------------------------
function results = userfcn_reserves_printpf(results, fd, mpopt, args)
%
%   results = userfcn_reserves_printpf(results, fd, mpopt, args)
%
%   This is the 'printpf' stage userfcn callback that pretty-prints the
%   results. It expects a results struct, a file descriptor and a MATPOWER
%   options struct. The optional args are not currently used.

%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%%-----  print results  -----
r = results.reserves;
nrz = size(r.req, 1);
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
OUT_RES         = OUT_ALL == 1 || (OUT_ALL == -1 && ~SUPPRESS && (mpopt.out.bus || mpopt.out.gen));

if isOPF && OUT_RES && (results.success || OUT_FORCE)
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n|     Reserves                                                                 |');
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n Gen   Bus   Status  Reserves   Price');
    fprintf(fd, '\n  #     #              (MW)     ($/MW)     Included in Zones ...');
    fprintf(fd, '\n----  -----  ------  --------  --------   ------------------------');
    for k = r.igr
        iz = find(r.zones(:, k));
        fprintf(fd, '\n%3d %6d     %2d ', k, results.gen(k, GEN_BUS), results.gen(k, GEN_STATUS));
        if results.gen(k, GEN_STATUS) > 0 && abs(results.reserves.R(k)) > 1e-6
            fprintf(fd, '%10.2f', results.reserves.R(k));
        else
            fprintf(fd, '       -  ');
        end
        fprintf(fd, '%10.2f     ', results.reserves.prc(k));
        for i = 1:length(iz)
            if i ~= 1
                fprintf(fd, ', ');
            end
            fprintf(fd, '%d', iz(i));
        end
    end
    fprintf(fd, '\n                     --------');
    fprintf(fd, '\n            Total:%10.2f              Total Cost: $%.2f', ...
        sum(results.reserves.R(r.igr)), results.reserves.totalcost);
    fprintf(fd, '\n');
    
    fprintf(fd, '\nZone  Reserves   Price  ');
    fprintf(fd, '\n  #     (MW)     ($/MW) ');
    fprintf(fd, '\n----  --------  --------');
    for k = 1:nrz
        iz = find(r.zones(k, :));     %% gens in zone k
        fprintf(fd, '\n%3d%10.2f%10.2f', k, sum(results.reserves.R(iz)), ...
                    results.lin.mu.l.Rreq(k) / results.baseMVA);
    end
    fprintf(fd, '\n');
    
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n|     Reserve Limits                                                           |');
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n Gen   Bus   Status  Rmin mu     Rmin    Reserves    Rmax    Rmax mu   Pmax mu ');
    fprintf(fd, '\n  #     #             ($/MW)     (MW)      (MW)      (MW)     ($/MW)    ($/MW) ');
    fprintf(fd, '\n----  -----  ------  --------  --------  --------  --------  --------  --------');
    for k = r.igr
        fprintf(fd, '\n%3d %6d     %2d ', k, results.gen(k, GEN_BUS), results.gen(k, GEN_STATUS));
        if results.gen(k, GEN_STATUS) > 0 && results.reserves.mu.l(k) > 1e-6
            fprintf(fd, '%10.2f', results.reserves.mu.l(k));
        else
            fprintf(fd, '       -  ');
        end
        fprintf(fd, '%10.2f', results.reserves.Rmin(k));
        if results.gen(k, GEN_STATUS) > 0 && abs(results.reserves.R(k)) > 1e-6
            fprintf(fd, '%10.2f', results.reserves.R(k));
        else
            fprintf(fd, '       -  ');
        end
        fprintf(fd, '%10.2f', results.reserves.Rmax(k));
        if results.gen(k, GEN_STATUS) > 0 && results.reserves.mu.u(k) > 1e-6
            fprintf(fd, '%10.2f', results.reserves.mu.u(k));
        else
            fprintf(fd, '       -  ');
        end
        if results.gen(k, GEN_STATUS) > 0 && results.reserves.mu.Pmax(k) > 1e-6
            fprintf(fd, '%10.2f', results.reserves.mu.Pmax(k));
        else
            fprintf(fd, '       -  ');
        end
    end
    fprintf(fd, '\n                                         --------');
    fprintf(fd, '\n                                Total:%10.2f', sum(results.reserves.R(r.igr)));
    fprintf(fd, '\n');
end


%%-----  savecase  -----------------------------------------------------
function mpc = userfcn_reserves_savecase(mpc, fd, prefix, args)
%
%   mpc = userfcn_reserves_savecase(mpc, fd, prefix, args)
%
%   This is the 'savecase' stage userfcn callback that prints the M-file
%   code to save the 'reserves' field in the case file. It expects a
%   MATPOWER case struct (mpc), a file descriptor and variable prefix
%   (usually 'mpc.'). The optional args are not currently used.

r = mpc.reserves;

fprintf(fd, '\n%%%%-----  Reserve Data  -----%%%%\n');
fprintf(fd, '%%%% reserve zones, element i, j is 1 if gen j is in zone i, 0 otherwise\n');
fprintf(fd, '%sreserves.zones = [\n', prefix);
template = '';
for i = 1:size(r.zones, 2)
    template = [template, '\t%d'];
end
template = [template, ';\n'];
fprintf(fd, template, r.zones.');
fprintf(fd, '];\n');

fprintf(fd, '\n%%%% reserve requirements for each zone in MW\n');
fprintf(fd, '%sreserves.req = [\t%g', prefix, r.req(1));
if length(r.req) > 1
    fprintf(fd, ';\t%g', r.req(2:end));
end
fprintf(fd, '\t];\n');

fprintf(fd, '\n%%%% reserve costs in $/MW for each gen that belongs to at least 1 zone\n');
fprintf(fd, '%%%% (same order as gens, but skipping any gen that does not belong to any zone)\n');
fprintf(fd, '%sreserves.cost = [\t%g', prefix, r.cost(1));
if length(r.cost) > 1
    fprintf(fd, ';\t%g', r.cost(2:end));
end
fprintf(fd, '\t];\n');

if isfield(r, 'qty')
    fprintf(fd, '\n%%%% OPTIONAL max reserve quantities for each gen that belongs to at least 1 zone\n');
    fprintf(fd, '%%%% (same order as gens, but skipping any gen that does not belong to any zone)\n');
    fprintf(fd, '%sreserves.qty = [\t%g', prefix, r.qty(1));
    if length(r.qty) > 1
        fprintf(fd, ';\t%g', r.qty(2:end));
    end
    fprintf(fd, '\t];\n');
end

%% save output fields for solved case
if isfield(r, 'R')
    if exist('serialize', 'file') == 2
        fprintf(fd, '\n%%%% solved values\n');
        fprintf(fd, '%sreserves.R = %s\n', prefix, serialize(r.R));
        fprintf(fd, '%sreserves.Rmin = %s\n', prefix, serialize(r.Rmin));
        fprintf(fd, '%sreserves.Rmax = %s\n', prefix, serialize(r.Rmax));
        fprintf(fd, '%sreserves.mu.l = %s\n', prefix, serialize(r.mu.l));
        fprintf(fd, '%sreserves.mu.u = %s\n', prefix, serialize(r.mu.u));
        fprintf(fd, '%sreserves.prc = %s\n', prefix, serialize(r.prc));
        fprintf(fd, '%sreserves.totalcost = %s\n', prefix, serialize(r.totalcost));
    else
        url = 'http://www.mathworks.com/matlabcentral/fileexchange/12063';
        warning('MATPOWER:serialize', ...
            'userfcn_reserves_savecase: Cannot save the ''reserves'' output fields without the ''serialize'' function, which is available as a free download from:\n<%s>\n\n', url);
    end
end
