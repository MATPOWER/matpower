function t_pf_radial(quiet)
%T_PF_RADIAL  Tests for distribution power flow solvers.

if nargin < 1
    quiet = 0;
end

t_begin(329, quiet);

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
mpopt = mpoption(mpopt, 'pf.alg', 'NR', 'pf.radial.max_it', 100);
casefile = {
    'case4_dist'
    'case18'
    'case22'
    'case33bw'
    'case69'
    'case85'
    'case141'
    };
method = {
    'PQSUM' 'Power Summation'
    'ISUM'  'Current Summation'
    'YSUM'  'Admittance Summation'
    };

% Original test cases (no PV buses)
iter_expected = {
           'none'    'NR' 'PQSUM' 'ISUM' 'YSUM'
     'case4_dist'      3      12     12     11
         'case18'      4       9     12     10
         'case22'      3       3      6      6
       'case33bw'      3       5      8      8
         'case69'      4       5      8      8
         'case85'      4       5      9      9
        'case141'      3       4      7      7
        };
iter_nopv = zeros(size(casefile,1),size(method,1)+1);
for i = 1:length(casefile)
    % which row in iter_expected is for current casefile
    row = find(strcmp(casefile{i}, iter_expected(:,1)) == 1);
    mpc = loadcase(casefile{i});
    % solve it with NR
    r = runpf(mpc, mpopt);
    iter_nopv(i,1) = r.iterations;
    % which columng in iter_expected is for current method
    col = find(strcmp('NR', iter_expected(1,:)) == 1);
    t = ['Newton''s method, ' casefile{i} ' : '];
    t_is(r.iterations, iter_expected{row,col}, 6, [t 'iterations']);
    for j = 1:size(method,1)
        mpopt1 = mpoption(mpopt,'pf.alg',method{j,1});
        r1 = runpf(mpc,mpopt1);
        iter_nopv(i,j+1) = r1.iterations;
        % which columng in iter_expected is for current method
        col = find(strcmp(method{j,1}, iter_expected(1,:)) == 1);
        t = [method{j,2} ', ' casefile{i} ' : '];
        t_ok(r1.success, [t 'success']);
        t_is(r1.iterations, iter_expected{row,col}, 6, [t 'iterations']);
        t_is(r1.bus, r.bus, 6, [t 'bus']);
        t_is(r1.gen, r.gen, 6, [t 'gen']);
        t_is(r1.branch, r.branch, 6, [t 'branch']);
    end
end

% Test cases with added PV buses
clear iter_expected
iter_expected{1} = {
% pf.radial.vcorr = 0
           'none'    'NR' 'PQSUM' 'ISUM' 'YSUM'
     'case4_dist'      3      12     12     11
         'case18'      4      20     16     13
         'case22'      4      11     25     25
       'case33bw'      4       8     16     16
         'case69'      4      12     27     27
         'case85'      5       9     17     17
        'case141'      5       9     18     18
};
iter_expected{2} = {
% pf.radial.vcorr = 1
           'none'    'NR' 'PQSUM' 'ISUM' 'YSUM'
     'case4_dist'      3       7      6      7
         'case18'      4      10     12     12
         'case22'      4      12     14     14
       'case33bw'      4       8     11     11
         'case69'      4      13     17     17
         'case85'      5       9     13     13
        'case141'      5       9     12     12
};
iter_pv = zeros(size(casefile,1),size(method,1)+1,2);
for i = 1:length(casefile)
    % which row in iter_expected is for current casefile
    row = find(strcmp(casefile{i}, iter_expected{1}(:,1)) == 1);
    mpc = loadcase(casefile{i});
    pv = mpc.bus(:,BUS_TYPE) == PV;
    if all(~pv)
        %%% ADD PV BUSES %%%
        N = 5; % How many PV buses should I add?
        % Solve the case without PV buses with Newton method
        r = runpf(mpc, mpopt);
        % Find the last N buses with lowest voltage
        [Vm, B] = sort(r.bus(:,VM));
         B = B(1:N);
        Vm = Vm(1:N);
        % add PV generators at buses in B
        % set VG 0.05 pu bigger then voltage at buses in B
        mpc.gen = repmat(mpc.gen,N+1,1);
        mpc.gen(2:end,GEN_BUS) = B;
        mpc.gen(2:end,VG) = Vm + 0.05;
        mpc.bus(B,BUS_TYPE) = 2;
        if isfield(mpc,'gencost')
            mpc.gencost = repmat(mpc.gencost,N+1,1);
        end
        %%% END OF ADD PV BUSES %%%
    end
    r = runpf(mpc, mpopt);
    iter_pv(i,1,1) = r.iterations;
    iter_pv(i,1,2) = r.iterations;
    % which columng in iter_expected is for current method
    col = find(strcmp('NR', iter_expected{1}(1,:)) == 1);
    t = ['Newton''s method, ' casefile{i} ' : '];
    t_is(r.iterations, iter_expected{1}{row,col}, 6, [t 'iterations']);
    for j = 1:size(method,1)
        for vcorr = 0:1
            mpopt1 = mpoption(mpopt,'pf.alg',method{j,1},'pf.radial.vcorr',vcorr);
            r1 = runpf(mpc,mpopt1);
            iter_pv(i,j+1,vcorr+1) = r1.iterations;
            % which columng in iter_expected is for current method
            col = find(strcmp(method{j,1}, iter_expected{vcorr+1}(1,:)) == 1);
            t = [method{j,2} ', vcorr = ' num2str(vcorr) ', ' casefile{i} ' : '];
            t_ok(r1.success, [t 'success']);
            t_is(r1.iterations, iter_expected{vcorr+1}{row,col}, 6, [t 'iterations']);
            t_is(r1.bus, r.bus, 6, [t 'bus']);
            t_is(r1.gen, r.gen, 6, [t 'gen']);
            t_is(r1.branch, r.branch, 6, [t 'branch']);
        end
    end
end

t_end;

if verbose > 0
    fprintf('\nITERATIONS: Original test cases (no PV buses), pf.radial.vcorr = 0\n');
    print_iterations(casefile,method,iter_nopv)
    fprintf('\nITERATIONS: Test cases with added PV buses, pf.radial.vcorr = 0\n');
    print_iterations(casefile,method,iter_pv(:,:,1))
    fprintf('\nITERATIONS: Test cases with added PV buses, pf.radial.vcorr = 1\n');
    print_iterations(casefile,method,iter_pv(:,:,2))
end

if have_fcn('octave')
    warning(s1.state, file_in_path_warn_id);
end

function print_iterations(casefile,method,iterations)
fprintf('%22s','NR');
for j = 1:size(method,1)
    fprintf('%7s',method{j,1});
end
fprintf('\n');
for i = 1:length(casefile)
    fprintf('%15s',casefile{i});
    for j = 1:size(method,1)+1
        fprintf('%7i',iterations(i,j));
    end
    fprintf('\n');
end
