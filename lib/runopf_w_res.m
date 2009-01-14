function [varargout] = runopf_w_res(casename, mpopt, fname, solvedcase)
%RUNOPF_W_RES  Runs an optimal power flow with fixed zonal reserves.
%
%   results = runopf_w_res(casename, mpopt, fname, solvedcase)
%   [results, success] = runopf_w_res(casename, mpopt, fname, solvedcase)
%
%   The case file or struct must define a userfcn field which is a struct
%   with fields as follows:
%       userfcn.name = 'userfcn_reserves'
%       userfcn.args = struct( ... )
%   The 'name' field defines the name of a function that will be called by
%   opf() to construct additional vars, constraints, costs for the OPF. The
%   value must be 'userfcn_reserves'.
%   The 'arg' field is an argument that will be passed to the user function.
%   See 'help userfcn_reserves' for details on the value of args expected
%   by userfcn_reserves(). See 'case30_reserves.m' for an example case file
%   with fixed reserves.
%
%   Note: runopf() can solve the same case correctly, the only difference
%   is that runopf_w_res() checks to make sure the case references
%   userfcn_reserves() and it does some post processing to package and print
%   the reserve results more nicely.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2008 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% default arguments
if nargin < 4
    solvedcase = '';                %% don't save solved case
    if nargin < 3
        fname = '';                 %% don't print results to a file
        if nargin < 2
            mpopt = mpoption;       %% use default options
            if nargin < 1
                casename = 'case9'; %% default data file is 'case9.m'
            end
        end
    end
end

%%-----  load and check for reserve inputs  -----
mpc = loadcase(casename);
if ~isfield(mpc, 'userfcn') || ~isfield(mpc.userfcn, 'name') || ...
        ~isstr(mpc.userfcn.name) || ...
        ~strcmp(mpc.userfcn.name, 'userfcn_reserves')
    error('runopf_w_res: case must contain a userfcn field with userfcn.name = ''userfcn_reserves''');
end
if ~isfield(mpc, 'userfcn') || ~isfield(mpc.userfcn, 'args')
    error('runopf_w_res: case must contain a userfcn field with a userfcn.args struct defining the reserve zones, requirements and costs');
end

%%-----  run it  -----
results = runopf(mpc, mpopt, fname, solvedcase);

%%-----  package up results  -----
args = mpc.userfcn.args;
nr = size(args.req, 1);
ng = size(results.gen, 1);

%% find gens with reserves
if nr > 1
    rgens = any(args.zones);
else
    rgens = args.zones;
end
igr = find(rgens);

%% adjust reserve zones, costs, qty, etc for off-line gens and gen reordering
rgens2 = rgens(results.status.gen.on(results.reorder.gen));
igr2 = find(rgens2);
ng2 = length(results.status.gen.on);

%% per gen reserves, multipliers with internal gen indexing
[R0, Rl, Ru] = getv(results.om, 'R');
Rg      = zeros(ng2, 1);
Rgmin   = zeros(ng2, 1);
Rgmax   = zeros(ng2, 1);
mu_l    = zeros(ng2, 1);
mu_u    = zeros(ng2, 1);
mu_Pmax = zeros(ng2, 1);
Rg(igr2)      = results.var.val.R * mpc.baseMVA;
Rgmin(igr2)   = Rl * mpc.baseMVA;
Rgmax(igr2)   = Ru * mpc.baseMVA;
mu_l(igr2)    = results.var.mu.l.R;
mu_u(igr2)    = results.var.mu.u.R;
mu_Pmax(igr2) = results.lin.mu.u.Pg_plus_R;

%% store in results
results.reserves.R       = zeros(ng, 1);
results.reserves.mu.l    = zeros(ng, 1);
results.reserves.mu.u    = zeros(ng, 1);
results.reserves.mu.Pmax = zeros(ng, 1);
results.reserves.R(results.status.gen.on)       = Rg(results.reorder.invgen);
results.reserves.mu.l(results.status.gen.on)    = mu_l(results.reorder.invgen);
results.reserves.mu.u(results.status.gen.on)    = mu_u(results.reorder.invgen);
results.reserves.mu.Pmax(results.status.gen.on) = mu_Pmax(results.reorder.invgen);
results.reserves.prc = zeros(ng, 1);
for k = igr
    iz = find(args.zones(:, k));
    results.reserves.prc(k) = max(results.lin.mu.l.Rreq(iz));
end

Rmin = zeros(ng, 1);
Rmax = zeros(ng, 1);
Rmin(results.status.gen.on) = Rgmin(results.reorder.invgen);
Rmax(results.status.gen.on) = Rgmax(results.reorder.invgen);

%%-----  print results  -----
if mpopt(32) ~= 0
    fd = 1;
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n|     Reserves                                                                 |');
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n Gen   Bus   Status  Reserves   Price');
    fprintf(fd, '\n  #     #              (MW)     ($/MW)     Included in Zones ...');
    fprintf(fd, '\n----  -----  ------  --------  --------   ------------------------');
    for k = igr
        iz = find(args.zones(:, k));
        fprintf(fd, '\n%3d %6d     %2d ', k, results.gen(k, GEN_BUS), results.gen(k, GEN_STATUS));
        if results.gen(k, GEN_STATUS) > 0 & abs(results.reserves.R(k)) > 1e-6
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
    fprintf(fd, '\n            Total:%10.2f', sum(results.reserves.R(igr)));
    fprintf(fd, '\n');
    
    fprintf(fd, '\nZone  Reserves   Price  ');
    fprintf(fd, '\n  #     (MW)     ($/MW) ');
    fprintf(fd, '\n----  --------  --------');
    for k = 1:nr
        iz = find(args.zones(k, :));     %% gens in zone k
        fprintf(fd, '\n%3d%10.2f%10.2f', k, sum(results.reserves.R(iz)), results.lin.mu.l.Rreq(k));
    end
    fprintf(fd, '\n');
    
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n|     Reserve Limits                                                           |');
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n Gen   Bus   Status  Rmin mu     Rmin    Reserves    Rmax    Rmax mu   Pmax mu ');
    fprintf(fd, '\n  #     #             ($/MW)     (MW)      (MW)      (MW)     ($/MW)    ($/MW) ');
    fprintf(fd, '\n----  -----  ------  --------  --------  --------  --------  --------  --------');
    for k = igr
        fprintf(fd, '\n%3d %6d     %2d ', k, results.gen(k, GEN_BUS), results.gen(k, GEN_STATUS));
        if results.gen(k, GEN_STATUS) > 0 & results.reserves.mu.l(k) > 1e-6
            fprintf(fd, '%10.2f', results.reserves.mu.l(k));
        else
            fprintf(fd, '       -  ');
        end
        fprintf(fd, '%10.2f', Rmin(k));
        if results.gen(k, GEN_STATUS) > 0 & abs(results.reserves.R(k)) > 1e-6
            fprintf(fd, '%10.2f', results.reserves.R(k));
        else
            fprintf(fd, '       -  ');
        end
        fprintf(fd, '%10.2f', Rmax(k));
        if results.gen(k, GEN_STATUS) > 0 & results.reserves.mu.u(k) > 1e-6
            fprintf(fd, '%10.2f', results.reserves.mu.u(k));
        else
            fprintf(fd, '       -  ');
        end
        if results.gen(k, GEN_STATUS) > 0 & results.reserves.mu.Pmax(k) > 1e-6
            fprintf(fd, '%10.2f', results.reserves.mu.Pmax(k));
        else
            fprintf(fd, '       -  ');
        end
    end
    fprintf(fd, '\n                                         --------');
    fprintf(fd, '\n                                Total:%10.2f', sum(results.reserves.R(igr)));
    fprintf(fd, '\n');
end

if nargout > 0
    varargout{1} = results;
    if nargout > 1
        varargout{2} = results.success;
    end
end

return;
