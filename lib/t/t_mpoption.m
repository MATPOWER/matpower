function t_mpoption(quiet)
%T_MPOPTION  Tests for MPOPTION.

%   MATPOWER
%   Copyright (c) 2013-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

v = 19;

t_begin(152, quiet);

%% default options struct
t = 'mpoption() : ';
mpopt = mpoption();
t_ok(isstruct(mpopt), [t 'isstruct mpopt']);
t_ok(isfield(mpopt, 'v'), [t ' isfield mpopt.v']);
t_is(mpopt.v, v, 12, [t '         mpopt.v == ' sprintf('%d', v)]);

t_ok(isfield(mpopt, 'model'), [t ' isfield mpopt.model']);
t_ok(strcmp(mpopt.model, 'AC'), [t '         mpopt.model = ''AC''']);

t_ok(isfield(mpopt, 'pf'), [t ' isfield mpopt.pf']);
t_ok(isstruct(mpopt.pf), [t 'isstruct mpopt.pf']);
t_ok(isfield(mpopt.pf, 'alg'), [t ' isfield mpopt.pf.alg']);
t_ok(strcmp(mpopt.pf.alg, 'NR'), [t '         mpopt.pf.alg = ''NR''']);
t_ok(isfield(mpopt.pf, 'tol'), [t ' isfield mpopt.pf.tol']);
t_is(mpopt.pf.tol, 1e-8, 12, [t '         mpopt.pf.tol == 1e-8']);
t_ok(isfield(mpopt.pf, 'nr'), [t ' isfield mpopt.pf.nr']);
t_ok(isstruct(mpopt.pf.nr), [t 'isstruct mpopt.pf.nr']);
t_ok(isfield(mpopt.pf.nr, 'max_it'), [t ' isfield mpopt.pf.nr.max_it']);
t_is(mpopt.pf.nr.max_it, 10, 12, [t '         mpopt.pf.nr.max_it == 10']);
t_ok(isfield(mpopt.pf, 'fd'), [t ' isfield mpopt.pf.fd']);
t_ok(isstruct(mpopt.pf.fd), [t 'isstruct mpopt.pf.fd']);
t_ok(isfield(mpopt.pf.fd, 'max_it'), [t ' isfield mpopt.pf.fd.max_it']);
t_is(mpopt.pf.fd.max_it, 30, 12, [t '         mpopt.pf.fd.max_it == 30']);
t_ok(isfield(mpopt.pf, 'gs'), [t ' isfield mpopt.pf.gs']);
t_ok(isstruct(mpopt.pf.gs), [t 'isstruct mpopt.pf.gs']);
t_ok(isfield(mpopt.pf.gs, 'max_it'), [t ' isfield mpopt.pf.gs.max_it']);
t_is(mpopt.pf.gs.max_it, 1000, 12, [t '         mpopt.pf.gs.max_it == 1000']);
t_ok(isfield(mpopt.pf, 'enforce_q_lims'), [t ' isfield mpopt.pf.enforce_q_lims']);
t_is(mpopt.pf.enforce_q_lims, 0, 12, [t '         mpopt.pf.enforce_q_lims == 0']);

t_ok(isfield(mpopt, 'opf'), [t ' isfield mpopt.opf']);
t_ok(isstruct(mpopt.opf), [t 'isstruct mpopt.opf']);
t_ok(isfield(mpopt.opf, 'ac'), [t ' isfield mpopt.opf.ac']);
t_ok(isstruct(mpopt.opf.ac), [t 'isstruct mpopt.opf.ac']);
t_ok(isfield(mpopt.opf.ac, 'solver'), [t ' isfield mpopt.opf.ac.solver']);
t_ok(strcmp(mpopt.opf.ac.solver, 'DEFAULT'), [t '         mpopt.opf.ac.solver = ''DEFAULT''']);
t_ok(isfield(mpopt.opf, 'dc'), [t ' isfield mpopt.opf.dc']);
t_ok(isstruct(mpopt.opf.dc), [t 'isstruct mpopt.opf.dc']);
t_ok(isfield(mpopt.opf.dc, 'solver'), [t ' isfield mpopt.opf.dc.solver']);
t_ok(strcmp(mpopt.opf.dc.solver, 'DEFAULT'), [t '         mpopt.opf.dc.solver = ''DEFAULT''']);
t_ok(isfield(mpopt.opf, 'violation'), [t ' isfield mpopt.opf.violation']);
t_is(mpopt.opf.violation, 5e-6, 12, [t '         mpopt.opf.violation == 5e-6']);
t_ok(isfield(mpopt.opf, 'flow_lim'), [t ' isfield mpopt.opf.flow_lim']);
t_ok(strcmp(upper(mpopt.opf.flow_lim), 'S'), [t '         mpopt.opf.flow_lim = ''S''']);
t_ok(isfield(mpopt.opf, 'ignore_angle_lim'), [t ' isfield mpopt.opf.ignore_angle_lim']);
t_is(mpopt.opf.ignore_angle_lim, 0, 12, [t '         mpopt.opf.ignore_angle_lim == 0']);
t_ok(isfield(mpopt.opf, 'return_raw_der'), [t ' isfield mpopt.opf.return_raw_der']);
t_is(mpopt.opf.return_raw_der, 0, 12, [t '         mpopt.opf.return_raw_der == 0']);

t_ok(isfield(mpopt, 'out'), [t ' isfield mpopt.out']);
t_ok(isstruct(mpopt.out), [t 'isstruct mpopt.out']);
t_ok(isfield(mpopt.out, 'all'), [t ' isfield mpopt.out.all']);
t_is(mpopt.out.all, -1, 12, [t '         mpopt.out.all == -1']);
t_ok(isfield(mpopt.out, 'sys_sum'), [t ' isfield mpopt.out.sys_sum']);
t_is(mpopt.out.sys_sum, 1, 12, [t '         mpopt.out.sys_sum == 1']);
t_ok(isfield(mpopt.out, 'area_sum'), [t ' isfield mpopt.out.area_sum']);
t_is(mpopt.out.area_sum, 0, 12, [t '         mpopt.out.area_sum == 0']);
t_ok(isfield(mpopt.out, 'bus'), [t ' isfield mpopt.out.bus']);
t_is(mpopt.out.bus, 1, 12, [t '         mpopt.out.bus == 1']);
t_ok(isfield(mpopt.out, 'branch'), [t ' isfield mpopt.out.branch']);
t_is(mpopt.out.branch, 1, 12, [t '         mpopt.out.branch == 1']);
t_ok(isfield(mpopt.out, 'gen'), [t ' isfield mpopt.out.gen']);
t_is(mpopt.out.gen, 0, 12, [t '         mpopt.out.gen == 0']);
t_ok(isfield(mpopt.out, 'lim'), [t ' isfield mpopt.out.lim']);
t_ok(isstruct(mpopt.out.lim), [t 'isstruct mpopt.out.lim']);
t_ok(isfield(mpopt.out.lim, 'all'), [t ' isfield mpopt.out.lim.all']);
t_is(mpopt.out.lim.all, -1, 12, [t '         mpopt.out.lim.all == -1']);
t_ok(isfield(mpopt.out.lim, 'v'), [t ' isfield mpopt.out.lim.v']);
t_is(mpopt.out.lim.v, 1, 12, [t '         mpopt.out.lim.v == 1']);
t_ok(isfield(mpopt.out.lim, 'line'), [t ' isfield mpopt.out.lim.line']);
t_is(mpopt.out.lim.line, 1, 12, [t '         mpopt.out.lim.line == 1']);
t_ok(isfield(mpopt.out.lim, 'pg'), [t ' isfield mpopt.out.lim.pg']);
t_is(mpopt.out.lim.pg, 1, 12, [t '         mpopt.out.lim.pg == 1']);
t_ok(isfield(mpopt.out.lim, 'qg'), [t ' isfield mpopt.out.lim.qg']);
t_is(mpopt.out.lim.qg, 1, 12, [t '         mpopt.out.lim.qg == 1']);
t_ok(isfield(mpopt.out, 'force'), [t ' isfield mpopt.out.force']);
t_is(mpopt.out.force, 0, 12, [t '         mpopt.out.force == 0']);

t_ok(isfield(mpopt, 'verbose'), [t ' isfield mpopt.verbose']);
t_is(mpopt.verbose, 1, 12, [t '         mpopt.verbose == 1']);

t_ok(isfield(mpopt, 'exp'), [t ' isfield mpopt.exp']);
t_ok(isstruct(mpopt.exp), [t '         isstruct(mpopt.exp)']);
t_ok(isfield(mpopt.exp, 'sys_wide_zip_loads'), [t ' isfield mpopt.exp.sys_wide_zip_loads']);
t_ok(isstruct(mpopt.exp.sys_wide_zip_loads), [t '         isstruct(mpopt.exp.sys_wide_zip_loads)']);
t_ok(isfield(mpopt.exp.sys_wide_zip_loads, 'pw'), [t ' isfield mpopt.exp.sys_wide_zip_loads.pw']);
t_ok(isempty(mpopt.exp.sys_wide_zip_loads.pw), [t '         isempty(mpopt.exp.sys_wide_zip_loads.pw)']);
t_ok(isfield(mpopt.exp.sys_wide_zip_loads, 'qw'), [t ' isfield mpopt.exp.sys_wide_zip_loads.qw']);
t_ok(isempty(mpopt.exp.sys_wide_zip_loads.pw), [t '         isempty(mpopt.exp.sys_wide_zip_loads.pw)']);

mpopt0 = mpopt;

t = 'mpoption(mpoption(), []) == mpoption_old()';
mpopt_v = mpoption_old;
t_is(mpoption(mpopt0, []), mpopt_v, 12, t);

t = 'mpoption(ov) : ';
ov = struct('verbose', 3, 'model', 'DC', 'opf', struct('dc', struct('solver', 'MIPS')));
mpopt = mpoption(ov);
t_is(mpopt.verbose, 3, 12, [t 'mpopt.verbose']);
t_ok(strcmp(upper(mpopt.model), 'DC'), [t 'mpopt.model']);
t_ok(strcmp(mpopt.opf.dc.solver, 'MIPS'), [t 'mpopt.opf.dc.solver']);
mpopt.verbose = 1;
mpopt.model = 'AC';
mpopt.opf.dc.solver = 'DEFAULT';
%% The following line appears to work around a bizarre bug in MATLAB 7.0.4 (Mac)
%% that caused the next test, and subsequent 'everything else' tests, to fail
%% mysteriously (but not if a debugger breakpoint was set).
isequal(mpopt, mpopt0);
t_ok(isequal(mpopt, mpopt0), [t 'everything else']);

t = 'mpoption(''t_mpoption_ov'') : ';
mpopt = mpoption('t_mpoption_ov');
t_is(mpopt.verbose, 2, 12, [t 'mpopt.verbose']);
t_ok(strcmp(upper(mpopt.model), 'DC'), [t 'mpopt.model']);
t_ok(strcmp(mpopt.opf.dc.solver, 'CPLEX'), [t 'mpopt.opf.dc.solver']);
mpopt.verbose = 1;
mpopt.model = 'AC';
mpopt.opf.dc.solver = 'DEFAULT';
t_ok(isequal(mpopt, mpopt0), [t 'everything else']);

t = 'mpoption(<new-style pairs>) : ';
mpopt = mpoption('verbose', 3, 'opf.dc.solver', 'MIPS', 'model', 'DC');
t_is(mpopt.verbose, 3, 12, [t 'mpopt.verbose']);
t_ok(strcmp(upper(mpopt.model), 'DC'), [t 'mpopt.model']);
t_ok(strcmp(mpopt.opf.dc.solver, 'MIPS'), [t 'mpopt.opf.dc.solver']);
mpopt.verbose = 1;
mpopt.model = 'AC';
mpopt.opf.dc.solver = 'DEFAULT';
t_ok(isequal(mpopt, mpopt0), [t 'everything else']);

t = 'mpoption(<old-style pairs>) : ';
mpopt = mpoption('VERBOSE', 0, 'PF_DC', 1, 'OPF_ALG_DC', 250);
t_is(mpopt.verbose, 0, 12, [t 'mpopt.verbose']);
t_ok(strcmp(upper(mpopt.model), 'DC'), [t 'mpopt.model']);
t_ok(strcmp(mpopt.opf.dc.solver, 'MIPS'), [t 'mpopt.opf.dc.solver']);
t_is(mpopt.mips.step_control, 1, 12, [t 'mpopt.mips.step_control']);
mpopt.verbose = 1;
mpopt.model = 'AC';
mpopt.opf.dc.solver = 'DEFAULT';
mpopt.mips.step_control = 0;
mpopt = delete_missing_optional_fields(mpopt, mpopt0);
t_ok(isequal(mpopt, mpopt0), [t 'everything else']);

t = 'mpoption(default_mpopt_v) : ';
mpopt = mpoption(mpopt_v);
mpopt = delete_missing_optional_fields(mpopt, mpopt0);
t_ok(isequal(mpopt, mpopt0), t);

t = 'mpoption(mpopt_v) : ';
mpopt_v1 = mpoption_old('VERBOSE', 0, 'PF_DC', 1, 'OPF_ALG_DC', 250);
mpopt = mpoption(mpopt_v1);
t_is(mpopt.verbose, 0, 12, [t 'mpopt.verbose']);
t_ok(strcmp(upper(mpopt.model), 'DC'), [t 'mpopt.model']);
t_ok(strcmp(mpopt.opf.dc.solver, 'MIPS'), [t 'mpopt.opf.dc.solver']);
t_is(mpopt.mips.step_control, 1, 12, [t 'mpopt.mips.step_control']);
mpopt.verbose = 1;
mpopt.model = 'AC';
mpopt.opf.dc.solver = 'DEFAULT';
mpopt.mips.step_control = 0;
mpopt = delete_missing_optional_fields(mpopt, mpopt0);
t_ok(isequal(mpopt, mpopt0), [t 'everything else']);

t = 'mpoption(mpopt_v, ov) : ';
ov = struct('verbose', 3, 'model', 'AC', 'opf', struct('ac', struct('solver', 'MIPS')));
mpopt = mpoption(mpopt_v1(1:116), ov);      %% use only first 116 elements, i.e MATPOWER 4.0 vector
t_is(mpopt.verbose, 3, 12, [t 'mpopt.verbose']);
t_ok(strcmp(upper(mpopt.model), 'AC'), [t 'mpopt.model']);
t_ok(strcmp(mpopt.opf.dc.solver, 'MIPS'), [t 'mpopt.opf.dc.solver']);
t_is(mpopt.mips.step_control, 1, 12, [t 'mpopt.mips.step_control']);
t_ok(strcmp(mpopt.opf.ac.solver, 'MIPS'), [t 'mpopt.opf.ac.solver']);
mpopt.verbose = 1;
mpopt.opf.dc.solver = 'DEFAULT';
mpopt.mips.step_control = 0;
mpopt.opf.ac.solver = 'DEFAULT';
mpopt = delete_missing_optional_fields(mpopt, mpopt0);
t_ok(isequal(mpopt, mpopt0), [t 'everything else']);

t = 'mpoption(mpopt_v, ''t_mpoption_ov'') : ';
mpopt = mpoption(mpopt_v1(1:93), 't_mpoption_ov');      %% use only first 93 elements, i.e MATPOWER 3.2 vector
t_is(mpopt.verbose, 2, 12, [t 'mpopt.verbose']);
t_ok(strcmp(upper(mpopt.model), 'DC'), [t 'mpopt.model']);
t_ok(strcmp(mpopt.opf.dc.solver, 'CPLEX'), [t 'mpopt.opf.dc.solver']);
t_is(mpopt.mips.step_control, 1, 12, [t 'mpopt.mips.step_control']);
mpopt.verbose = 1;
mpopt.model = 'AC';
mpopt.opf.dc.solver = 'DEFAULT';
mpopt.mips.step_control = 0;
mpopt = delete_missing_optional_fields(mpopt, mpopt0);
t_ok(isequal(mpopt, mpopt0), [t 'everything else']);

t = 'mpoption(mpopt_v, <new-style pairs>) : ';
mpopt = mpoption(mpopt_v1, 'verbose', 3, 'model', 'AC', 'opf.ac.solver', 'MIPS');
t_is(mpopt.verbose, 3, 12, [t 'mpopt.verbose']);
t_ok(strcmp(upper(mpopt.model), 'AC'), [t 'mpopt.model']);
t_ok(strcmp(mpopt.opf.dc.solver, 'MIPS'), [t 'mpopt.opf.dc.solver']);
t_is(mpopt.mips.step_control, 1, 12, [t 'mpopt.mips.step_control']);
t_ok(strcmp(mpopt.opf.ac.solver, 'MIPS'), [t 'mpopt.opf.ac.solver']);
mpopt.verbose = 1;
mpopt.opf.dc.solver = 'DEFAULT';
mpopt.mips.step_control = 0;
mpopt.opf.ac.solver = 'DEFAULT';
mpopt = delete_missing_optional_fields(mpopt, mpopt0);
t_ok(isequal(mpopt, mpopt0), [t 'everything else']);

t = 'mpoption(mpopt_v, <old-style pairs>) : ';
mpopt = mpoption(mpopt_v1, 'VERBOSE', 3, 'PF_DC', 0, 'OPF_ALG', 520);
t_is(mpopt.verbose, 3, 12, [t 'mpopt.verbose']);
t_ok(strcmp(upper(mpopt.model), 'AC'), [t 'mpopt.model']);
t_ok(strcmp(mpopt.opf.dc.solver, 'MIPS'), [t 'mpopt.opf.dc.solver']);
t_is(mpopt.mips.step_control, 1, 12, [t 'mpopt.mips.step_control']);
t_ok(strcmp(mpopt.opf.ac.solver, 'FMINCON'), [t 'mpopt.opf.ac.solver']);
mpopt.verbose = 1;
mpopt.opf.dc.solver = 'DEFAULT';
mpopt.mips.step_control = 0;
mpopt.opf.ac.solver = 'DEFAULT';
mpopt = delete_missing_optional_fields(mpopt, mpopt0);
t_ok(isequal(mpopt, mpopt0), [t 'everything else']);

t = 'mpoption(mpopt) : ';
mpopt1 = mpoption('verbose', 0, 'model', 'DC', 'out.gen', 1, 'opf.dc.solver', 'GUROBI');
mpopt = mpoption(mpopt1);
t_ok(isequal(mpopt, mpopt1), [t 'unchanged']);

t = 'mpoption(mpopt, ov) : ';
ov = struct('verbose', 3, 'model', 'AC', 'opf', struct('ac', struct('solver', 'KNITRO')));
mpopt = mpoption(mpopt1, ov);
t_is(mpopt.verbose, 3, 12, [t 'mpopt.verbose']);
t_ok(strcmp(upper(mpopt.model), 'AC'), [t 'mpopt.model']);
t_ok(strcmp(mpopt.opf.dc.solver, 'GUROBI'), [t 'mpopt.opf.dc.solver']);
t_is(mpopt.out.gen, 1, 12, [t 'mpopt.out.gen']);
t_ok(strcmp(mpopt.opf.ac.solver, 'KNITRO'), [t 'mpopt.opf.ac.solver']);
mpopt.verbose = 1;
mpopt.opf.dc.solver = 'DEFAULT';
mpopt.out.gen = 0;
mpopt.opf.ac.solver = 'DEFAULT';
t_ok(isequal(mpopt, mpopt0), [t 'everything else']);

t = 'mpoption(mpopt, ''t_mpoption_ov'') : ';
mpopt = mpoption(mpopt1, 't_mpoption_ov');
t_is(mpopt.verbose, 2, 12, [t 'mpopt.verbose']);
t_ok(strcmp(upper(mpopt.model), 'DC'), [t 'mpopt.model']);
t_ok(strcmp(mpopt.opf.dc.solver, 'CPLEX'), [t 'mpopt.opf.dc.solver']);
t_is(mpopt.out.gen, 1, 12, [t 'mpopt.out.gen']);
mpopt.verbose = 1;
mpopt.model = 'AC';
mpopt.opf.dc.solver = 'DEFAULT';
mpopt.out.gen = 0;
t_ok(isequal(mpopt, mpopt0), [t 'everything else']);

t = 'mpoption(mpopt, <new-style pairs>) : ';
mpopt = mpoption(mpopt1, 'verbose', 3, 'model', 'AC', 'opf.ac.solver', 'KNITRO');
t_is(mpopt.verbose, 3, 12, [t 'mpopt.verbose']);
t_ok(strcmp(upper(mpopt.model), 'AC'), [t 'mpopt.model']);
t_ok(strcmp(mpopt.opf.dc.solver, 'GUROBI'), [t 'mpopt.opf.dc.solver']);
t_is(mpopt.out.gen, 1, 12, [t 'mpopt.out.gen']);
t_ok(strcmp(mpopt.opf.ac.solver, 'KNITRO'), [t 'mpopt.opf.ac.solver']);
mpopt.verbose = 1;
mpopt.opf.dc.solver = 'DEFAULT';
mpopt.out.gen = 0;
mpopt.opf.ac.solver = 'DEFAULT';
t_ok(isequal(mpopt, mpopt0), [t 'everything else']);

t = 'mpoption(mpopt, <old-style pairs>) : ';
mpopt = mpoption(mpopt1, 'VERBOSE', 3, 'PF_DC', 0, 'OPF_ALG', 580);
t_is(mpopt.verbose, 3, 12, [t 'mpopt.verbose']);
t_ok(strcmp(upper(mpopt.model), 'AC'), [t 'mpopt.model']);
t_ok(strcmp(mpopt.opf.dc.solver, 'GUROBI'), [t 'mpopt.opf.dc.solver']);
t_is(mpopt.out.gen, 1, 12, [t 'mpopt.out.gen']);
t_ok(strcmp(mpopt.opf.ac.solver, 'IPOPT'), [t 'mpopt.opf.ac.solver']);
mpopt.verbose = 1;
mpopt.opf.dc.solver = 'DEFAULT';
mpopt.out.gen = 0;
mpopt.opf.ac.solver = 'DEFAULT';
t_ok(isequal(mpopt, mpopt0), [t 'everything else']);

t_end;


function opt = delete_missing_optional_fields(opt, unless)
%% deletes the fields from opt, unless they are found in unless
%% (which is empty by default)
%pkgs = {'cplex', 'sdp_pf', 'sopf', 'yalmip'};
pkgs = {...
    'cplex', 'fmincon', 'gurobi', 'ipopt', 'knitro', 'minopf', ...
    'mosek', 'sopf', 'pdipm', 'tralm', ...
};
if nargin < 2
	unless = struct;
end
for k = 1:length(pkgs)
	if ~isfield(unless, pkgs{k}) && isfield(opt, pkgs{k})
		opt = rmfield(opt, pkgs{k});
	end
end
