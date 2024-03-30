function task_rv = t_run_mp(quiet)
% t_run_mp - Tests for run_mp and simple creation and solve of models.

%   MATPOWER
%   Copyright (c) 2021-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 1
    quiet = 0;
end

% define_constants;
if quiet
    verbose = 0;
else
    verbose = 1;
end

t_begin(25, quiet);

casefile = 'case9';
casefilet = 'case9target';
opt = struct('verbose', 0);
mpopt = mpoption(opt);
mpopt = mpoption(mpopt, 'out.all', 0);

tasks = {'PF', 'CPF', 'OPF'};
task_classes = {@mp.task_pf, @mp.task_cpf, @mp.task_opf};
d = {casefile, {casefile, casefilet}, casefile};
dmc_class = @mp.dm_converter_mpc2;
dm_classes = {
    @mp.data_model,
    @mp.data_model_cpf,
    @mp.data_model_opf
};
nm_class = @mp.net_model_acp;
mm_classes = {
    @mp.math_model_pf_acps,
    @mp.math_model_cpf_acps,
    @mp.math_model_opf_acps
};

t = 'build data model converter';
dmc = dmc_class();
dmc.build();
t_ok(isa(dmc, 'mp.dm_converter'), t);

for k = 1:length(tasks)
    t = sprintf('%s : ', tasks{k});

    %% build data model
    dm_class = dm_classes{k};
    dm = dm_class();
    dm.build('case9', dmc);
    t_ok(isa(dm, 'mp.data_model'), [t 'build data model']);

    %% build network model
    nm = nm_class();
    nm.build(dm);
    t_ok(isa(nm, 'mp.net_model'), [t 'build network model']);

    if strcmp(tasks{k}, 'CPF')
        dm.userdata.target = dm_class();
        dm.userdata.target.build('case9target', dmc);
        nm.userdata.target = nm_class();
        nm.userdata.target.build(dm.userdata.target);
    end

    %% build math model
    mm_class = mm_classes{k};
    mm = mm_class();
    mm.build(nm, dm, mpopt);
    t_ok(isa(mm, 'mp.math_model'), [t 'build math model']);

    %% solve math model
    % opt = mm.solve_opts(nm, dm, mpopt);
    mm.solve(opt);
    t_is(mm.soln.eflag, 1, 12, [t 'solve math model']);

    %% network model solution
    nm = mm.network_model_x_soln(nm);
    t_ok(isfield(nm.soln, 'x'), [t 'network model x soln']);
    nm.port_inj_soln();
    t_ok(isfield(nm.soln, 'gs_'), [t 'network model port inj soln']);

    %% data model update
    dm = mm.data_model_update(nm, dm, mpopt);
    t_ok(all(dm.elements.branch.tab.pl_fr ~= 0), [t 'data model update']);

    %% run_mp
    task = run_mp(task_classes{k}, d{k}, mpopt);
    t_is(task.success, 1, 12, [t 'run_mp : success']);
end

if nargout  %% set output arg
    task_rv = task;
end

t_end;
