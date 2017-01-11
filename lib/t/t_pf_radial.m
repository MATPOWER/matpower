function t_pf_radial(quiet)
%T_PF_RADIAL  Tests for distribution power flow solvers.

if nargin < 1
    quiet = 0;
end

t_begin(84, quiet);

if quiet
    verbose = 0;
else
    verbose = 1;
end
if have_fcn('octave')
    if have_fcn('octave', 'vnum') >= 4
        file_in_path_warn_id = 'Octave:data-file-in-path';
    else
        file_in_path_warn_id = 'Octave:load-file-in-path';
    end
    s1 = warning('query', file_in_path_warn_id);
    warning('off', file_in_path_warn_id);
end
mpopt = mpoption('out.all', 0, 'verbose', verbose);

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% Test Distribution Power Flow
mpopt0 = mpoption(mpopt, 'pf.alg', 'NR');
mpopt1 = mpoption(mpopt, 'pf.alg', 'PQSUM');
mpopt2 = mpoption(mpopt, 'pf.alg', 'ISUM','pf.tol',1e-9);
mpopt3 = mpoption(mpopt, 'pf.alg', 'YSUM');
casefile = {
    'case4_dist'
    'case18'
    'case22'
    'case33bw'
    'case69'
    'case85'
    'case141'
    };
for i = 1:length(casefile)
    [baseMVA, bus0, gen0, branch0, success0, et0] = runpf(casefile{i}, mpopt0);
    [baseMVA, bus1, gen1, branch1, success1, et1] = runpf(casefile{i}, mpopt1);
    [baseMVA, bus2, gen2, branch2, success2, et2] = runpf(casefile{i}, mpopt2);
    [baseMVA, bus3, gen3, branch3, success3, et3] = runpf(casefile{i}, mpopt3);
    % Power Summation
    t = ['Power Summation, ' casefile{i} ' : '];
    t_ok(success1, [t 'success']);
    t_is(bus1, bus0, 6, [t 'bus']);
    t_is(gen1, gen0, 6, [t 'gen']);
    t_is(branch1, branch0, 6, [t 'branch']);
    % Current Summation
    t = ['Current Summation, ' casefile{i} ' : '];
    t_ok(success2, [t 'success']);
    t_is(bus2, bus0, 6, [t 'bus']);
    t_is(gen2, gen0, 6, [t 'gen']);
    t_is(branch2, branch0, 6, [t 'branch']);
    % Admittance Summation
    t = ['Admittance Summation, ' casefile{i} ' : '];
    t_ok(success3, [t 'success']);
    t_is(bus3, bus0, 6, [t 'bus']);
    t_is(gen3, gen0, 6, [t 'gen']);
    t_is(branch3, branch0, 6, [t 'branch']);
end

t_end;

if have_fcn('octave')
    warning(s1.state, file_in_path_warn_id);
end
