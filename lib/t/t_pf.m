function t_pf(quiet)
%T_PF  Tests for power flow solvers.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2004-2010 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

if nargin < 1
    quiet = 0;
end

t_begin(25, quiet);

casefile = 't_case9_pf';
if quiet
    verbose = 0;
else
    verbose = 1;
end
if have_fcn('octave')
	s1 = warning('query', 'Octave:load-file-in-path');
    warning('off', 'Octave:load-file-in-path');
end
mpopt = mpoption('OUT_ALL', 0, 'VERBOSE', verbose);

[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% get solved AC power flow case from MAT-file
load soln9_pf;      %% defines bus_soln, gen_soln, branch_soln

%% run Newton PF
t = 'Newton PF : ';
mpopt = mpoption(mpopt, 'PF_ALG', 1);
[baseMVA, bus, gen, branch, success, et] = runpf(casefile, mpopt);
t_ok(success, [t 'success']);
t_is(bus, bus_soln, 6, [t 'bus']);
t_is(gen, gen_soln, 6, [t 'gen']);
t_is(branch, branch_soln, 6, [t 'branch']);

%% run fast-decoupled PF (XB version)
t = 'Fast Decoupled (XB) PF : ';
mpopt = mpoption(mpopt, 'PF_ALG', 2);
[baseMVA, bus, gen, branch, success, et] = runpf(casefile, mpopt);
t_ok(success, [t 'success']);
t_is(bus, bus_soln, 6, [t 'bus']);
t_is(gen, gen_soln, 6, [t 'gen']);
t_is(branch, branch_soln, 6, [t 'branch']);

%% run fast-decoupled PF (BX version)
t = 'Fast Decoupled (BX) PF : ';
mpopt = mpoption(mpopt, 'PF_ALG', 3);
[baseMVA, bus, gen, branch, success, et] = runpf(casefile, mpopt);
t_ok(success, [t 'success']);
t_is(bus, bus_soln, 6, [t 'bus']);
t_is(gen, gen_soln, 6, [t 'gen']);
t_is(branch, branch_soln, 6, [t 'branch']);

%% run Gauss-Seidel PF
t = 'Gauss-Seidel PF : ';
mpopt = mpoption(mpopt, 'PF_ALG', 4);
[baseMVA, bus, gen, branch, success, et] = runpf(casefile, mpopt);
t_ok(success, [t 'success']);
t_is(bus, bus_soln, 5, [t 'bus']);
t_is(gen, gen_soln, 5, [t 'gen']);
t_is(branch, branch_soln, 5, [t 'branch']);

%% get solved AC power flow case from MAT-file
load soln9_dcpf;        %% defines bus_soln, gen_soln, branch_soln

%% run DC PF
t = 'DC PF : ';
[baseMVA, bus, gen, branch, success, et] = rundcpf(casefile, mpopt);
t_ok(success, [t 'success']);
t_is(bus, bus_soln, 6, [t 'bus']);
t_is(gen, gen_soln, 6, [t 'gen']);
t_is(branch, branch_soln, 6, [t 'branch']);

%% check Qg distribution, when Qmin = Qmax
t = 'check Qg : ';
mpopt = mpoption(mpopt, 'PF_ALG', 1, 'VERBOSE', 0);
mpc = loadcase(casefile);
mpc.gen(1, [QMIN QMAX]) = [20 20];
[baseMVA, bus, gen, branch, success, et] = runpf(mpc, mpopt);
t_is(gen(1, QG), 24.07, 2, [t 'single gen, Qmin = Qmax']);

mpc.gen = [mpc.gen(1, :); mpc.gen];
mpc.gen(1, [QMIN QMAX]) = [10 10];
mpc.gen(2, [QMIN QMAX]) = [0 50];
[baseMVA, bus, gen, branch, success, et] = runpf(mpc, mpopt);
t_is(gen(1:2, QG), [10; 14.07], 2, [t '2 gens, Qmin = Qmax for one']);

mpc.gen(1, [QMIN QMAX]) = [10 10];
mpc.gen(2, [QMIN QMAX]) = [-50 -50];
[baseMVA, bus, gen, branch, success, et] = runpf(mpc, mpopt);
t_is(gen(1:2, QG), [12.03; 12.03], 2, [t '2 gens, Qmin = Qmax for both']);

mpc.gen(1, [QMIN QMAX]) = [0 50];
mpc.gen(2, [QMIN QMAX]) = [0 100];
[baseMVA, bus, gen, branch, success, et] = runpf(mpc, mpopt);
t_is(gen(1:2, QG), [8.02; 16.05], 2, [t '2 gens, proportional']);

mpc.gen(1, [QMIN QMAX]) = [-50 0];
mpc.gen(2, [QMIN QMAX]) = [50 150];
[baseMVA, bus, gen, branch, success, et] = runpf(mpc, mpopt);
t_is(gen(1:2, QG), [-50+8.02; 50+16.05], 2, [t '2 gens, proportional']);

if have_fcn('octave')
    warning(s1.state, 'Octave:load-file-in-path');
end

t_end;
