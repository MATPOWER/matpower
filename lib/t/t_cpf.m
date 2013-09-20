function res = t_cpf(quiet)
%T_PF  Tests for continuation power flow.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2013 by Power System Engineering Research Center (PSERC)
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

t_begin(21, quiet);
plot_nose_curve = 0;

casefile = 't_case9_pfv2';
if quiet
    verbose = 0;
else
    verbose = 1;
end
if have_fcn('octave')
    s1 = warning('query', 'Octave:load-file-in-path');
    warning('off', 'Octave:load-file-in-path');
end
mpopt = mpoption('OUT_ALL', 0, 'VERBOSE', 0);
%mpopt = mpoption(mpopt, 'CPF_TRACE_FULL', 0);
mpopt = mpoption(mpopt, 'CPF_STEP', 0.02);
%mpopt = mpoption(mpopt, 'CPF_ADAPT_STEP', 1);
%mpopt = mpoption(mpopt, 'CPF_ERROR_TOL', 2e-5);
if plot_nose_curve
	mpopt = mpoption(mpopt, 'CPF_USER_CALLBACK', 1);
	mpopt = mpoption(mpopt, 'VERBOSE', 2);
end	

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% set up base and target cases
mpcb = loadcase(casefile);
%% add isolated bus to make sure int2ext works for V_p, V_c
mpcb.bus = [mpcb.bus(1:3, :); mpcb.bus(3, :); mpcb.bus(4:end, :)];
mpcb.bus(4, BUS_I) = 50;
mpcb.bus(4, BUS_TYPE) = NONE;
% r = runpf(mpcb, mpopt);
% mpcb.gen(1, [PG QG]) = r.gen(1, [PG QG]);	%% solved values for slack gen
mpct = mpcb;
factor = 2.5;
mpct.gen(:, [PG QG]) = mpct.gen(:, [PG QG]) * factor;
mpct.bus(:, [PD QD]) = mpct.bus(:, [PD QD]) * factor;

%% run CPF
t = 'Continuation PF to lambda = 0.7 : ';
mpopt = mpoption(mpopt, 'CPF_TRACE_LAMBDA', 0.7);
mpopt = mpoption(mpopt, 'CPF_TRACE_NOSE', 0);
r = runcpf(mpcb, mpct, mpopt);
iterations = 41;
t_ok(r.success, [t 'success']);
t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
t_is(r.cpf.max_lam, 0.707703, 4, [t 'max_lam']);
t_is(size(r.cpf.V_p), [10 iterations+1], 12, [t 'size(V_p)']);
t_is(size(r.cpf.V_c), [10 iterations+1], 12, [t 'size(V_c)']);
t_is(size(r.cpf.lam_p), [1 iterations+1], 12, [t 'size(lam_p)']);
t_is(size(r.cpf.lam_c), [1 iterations+1], 12, [t 'size(lam_c)']);

t = 'Continuation PF to nose point : ';
mpopt = mpoption(mpopt, 'CPF_TRACE_FULL', 0);
mpopt = mpoption(mpopt, 'CPF_TRACE_NOSE', 1);
mpopt = mpoption(mpopt, 'CPF_TRACE_LAMBDA', 0);
mpopt = mpoption(mpopt, 'CPF_ADAPT_STEP', 1);
r = runcpf(mpcb, mpct, mpopt);
iterations = 22;
t_ok(r.success, [t 'success']);
t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
t_is(r.cpf.max_lam, 0.99025, 3, [t 'max_lam']);
t_is(size(r.cpf.V_p), [10 iterations+1], 12, [t 'size(V_p)']);
t_is(size(r.cpf.V_c), [10 iterations+1], 12, [t 'size(V_c)']);
t_is(size(r.cpf.lam_p), [1 iterations+1], 12, [t 'size(lam_p)']);
t_is(size(r.cpf.lam_c), [1 iterations+1], 12, [t 'size(lam_c)']);

t = 'Continuation PF (full trace) : ';
mpopt = mpoption(mpopt, 'CPF_TRACE_FULL', 1);
mpopt = mpoption(mpopt, 'CPF_TRACE_NOSE', 0);
r = runcpf(mpcb, mpct, mpopt);
iterations = 45;
t_ok(r.success, [t 'success']);
t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
t_is(r.cpf.max_lam, 0.99025, 3, [t 'max_lam']);
t_is(size(r.cpf.V_p), [10 iterations+1], 12, [t 'size(V_p)']);
t_is(size(r.cpf.V_c), [10 iterations+1], 12, [t 'size(V_c)']);
t_is(size(r.cpf.lam_p), [1 iterations+1], 12, [t 'size(lam_p)']);
t_is(size(r.cpf.lam_c), [1 iterations+1], 12, [t 'size(lam_c)']);

if have_fcn('octave')
    warning(s1.state, 'Octave:load-file-in-path');
end

t_end;

if nargout
	res = r;
end