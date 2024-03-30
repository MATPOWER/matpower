function task = t_run_mp_3p(quiet)
% t_run_mp_3p - Tests for run_pf, run_cpf, run_opf for 3-phase and hybrid test cases.

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

casefiles = {
    't_case3p_a', ...
    't_case3p_b', ...
    't_case3p_c', ...
    't_case3p_d', ...
    't_case3p_e', ...
    't_case3p_f', ...
    't_case3p_g', ...
    't_case3p_h'  };

lam = [
    0.6540504347;
    0.6440398987;
    0.6540504347;
    0.6540504347;
    0.6540504347;
    0.6428930927;
    0.6466143248;
    0.6466969520;
];

vm1 = [
    1.0;
    1.0024291384;
    1.0024291384;
    1.0024291384;
    1.0;
    1.0193502204;
    1.0;
    1.0;
];

vm_end = [
    0.6158223964;
    0.6160907764;
    0.6158223965;
    0.6158223965;
    0.6158223965;
    0.6158556726;
    0.6160426241;
    0.6160460608;
];

%%  alg     name                            check   opts
cfg = {
    {'NR',  'Newton (polar-power)',         [],     {}  },
    {'NR',  'Newton (cartesian-power)',     [],     {'pf.v_cartesian', 1}  },
    {'NR',  'Newton (polar-current)',       [],     {'pf.current_balance', 1}  },
    {'NR',  'Newton (cartesian-current)',   [],     {'pf.v_cartesian', 1, 'pf.current_balance', 1}  },
};

%%  alg     name                            check   opts
cfg_cpf = {
    {'CPF',  'CPF (polar-power)',           [],     {}  },
    {'CPF',  'CPF (cartesian-power)',       [],     {'pf.v_cartesian', 1}  },
    {'CPF',  'CPF (polar-current)',         [],     {'pf.current_balance', 1}  },
    {'CPF',  'CPF (cartesian-current)',     [],     {'pf.v_cartesian', 1, 'pf.current_balance', 1}  },
};

%%  alg     name                            check   opts
cfg_opf = {
    {'MIPS',  'MIPS (polar-power)',         [],     {}  },
    {'MIPS',  'MIPS (cartesian-power)',     [],     {'opf.v_cartesian', 1}  },
    {'MIPS',  'MIPS (polar-current)',       [],     {'opf.current_balance', 1}  },
    {'MIPS',  'MIPS (cartesian-current)',   [],     {'opf.v_cartesian', 1, 'opf.current_balance', 1}  },
};

n_pf = 4 + 4 * length(casefiles)*length(cfg);
n_cpf = 4 * length(casefiles)*length(cfg_cpf);
n_opf = 4 * length(cfg_opf);

t_begin(n_pf + n_cpf + n_opf, quiet);

mpopt0 = mpoption('out.all', 0, 'verbose', 0, 'mips.comptol', 1e-9);

if have_feature('octave')
    sing_mat_warn_id = 'Octave:singular-matrix';
else
    sing_mat_warn_id = 'MATLAB:singularMatrix';
end
s = warning('query', sing_mat_warn_id);
warning('off', sing_mat_warn_id);
warning('off', 'update_z:multiple_nodes');

eva = 100 * [
     0.000000000033021  -1.199999999966978   1.200000000033022
    -0.001399911879082  -1.201847639219705   1.192648372538439
    -0.022580289275887  -1.236249929877042   1.147882346913731
    -0.041234097508356  -1.267980617541095   1.028457524322151 ];
evm = [
    7.199557856520505   7.199557856520505   7.199557856520503
    7.163739311071976   7.110502601043655   7.082050422574726
    2.305502785229222   2.254669273283851   2.202824420816035
    2.174983020357710   1.929841224321033   1.832702204181716 ];
epl = 1000 * [
   1.341424819811426   0.970523307120210   2.096102639009662   1.341414943637188   2.672344186177307   1.894590807224783  -1.337240567150252  -0.963372951103195  -2.074358800215885  -1.319145934210436  -2.652375608822653  -1.830632969045932
   1.323522914268509   0.881067034036327   2.043381640999108   1.133282979238420   2.598706416975550   1.508617820874939  -1.274999999960655  -0.790174031427826  -1.800000000002936  -0.871779788853622  -2.374999999391016  -0.780624749701872 ];

va_fields = {'va1', 'va2', 'va3'};
vm_fields = {'vm1', 'vm2', 'vm3'};
pl_fields = {'pl1_fr', 'ql1_fr', 'pl2_fr', 'ql2_fr', 'pl3_fr', 'ql3_fr', ...
             'pl1_to', 'ql1_to', 'pl2_to', 'ql2_to', 'pl3_to', 'ql3_to'};

%%-----  Power Flow  -----
%% Z-Gauss, 3-phase only
c = 1;
casefile = casefiles{c};
alg = 'ZG';
name = 'Z-Gauss (polar-power)';
opts = {};
t = sprintf('PF - %s : %s : ', casefile, name);
mpopt = mpoption(mpopt0, 'pf.alg', alg, opts{:});
pf = run_pf(casefile, mpopt, 'mpx', mp.xt_3p());

t_is(pf.success, 1, 12, [t 'success']);

va = cell2mat(cellfun(@(x)pf.dm.elements.bus3p.tab.(x), va_fields, 'UniformOutput', 0));
vm = cell2mat(cellfun(@(x)pf.dm.elements.bus3p.tab.(x), vm_fields, 'UniformOutput', 0)) .* ...
    (pf.dm.elements.bus3p.tab.base_kv*ones(1,3))/sqrt(3);
pl = cell2mat(cellfun(@(x)pf.dm.elements.line3p.tab.(x), pl_fields, 'UniformOutput', 0));

t_is(va, eva, 7, [t 'va']);
t_is(vm, evm, 7, [t 'vm']);
t_is(pl, epl, 4.8, [t 'pl']);

%% Newton-Raphson power flow, all cases
for k = 1:length(cfg)
    for c = 1:length(casefiles)
        casefile = casefiles{c};
        [alg, name, check, opts] = deal(cfg{k}{:});
        t = sprintf('PF - %s : %s : ', casefile, name);
        mpopt = mpoption(mpopt0, 'pf.alg', alg, opts{:});
        pf = run_pf(casefile, mpopt, 'mpx', mp.xt_3p());

        t_is(pf.success, 1, 12, [t 'success']);

        va = cell2mat(cellfun(@(x)pf.dm.elements.bus3p.tab.(x), va_fields, 'UniformOutput', 0));
        vm = cell2mat(cellfun(@(x)pf.dm.elements.bus3p.tab.(x), vm_fields, 'UniformOutput', 0)) .* ...
            (pf.dm.elements.bus3p.tab.base_kv*ones(1,3))/sqrt(3);
        pl = cell2mat(cellfun(@(x)pf.dm.elements.line3p.tab.(x), pl_fields, 'UniformOutput', 0));

        switch casefile
            case {'t_case3p_a', 't_case3p_b', 't_case3p_c', 't_case3p_d', 't_case3p_e'}
                t_is(va, eva, 7, [t 'va']);
                t_is(vm, evm, 7, [t 'vm']);
                t_is(pl, epl, 5, [t 'pl']);
            case 't_case3p_f'
                t_is(va, [eva; eva+3.570076720548507; eva+5.361141407429377], 7, [t 'va']);
                t_is(vm, [evm; evm; evm], 7, [t 'vm']);
                t_is(pl, [epl; epl; epl], 4.8, [t 'pl']);
            case 't_case3p_g'
                t_is(va, [eva; eva+3.228067320798516; eva+2.435814088796515], 7, [t 'va']);
                t_is(vm, [evm; evm; evm], 7, [t 'vm']);
                t_is(pl, [epl; epl; epl], 5, [t 'pl']);
            case 't_case3p_h'
                t_is(va, [eva; eva+4.194887190252387; eva+2.968367744185669], 7, [t 'va']);
                t_is(vm, [evm; evm; evm], 7, [t 'vm']);
                t_is(pl, [epl; epl; epl], 5, [t 'pl']);
        end
    end
end

%%-----  Continuation Power Flow  -----
for k = 1:length(cfg_cpf)
    for c = 1:length(casefiles)
        casefile = casefiles{c};
        [alg, name, check, opts] = deal(cfg_cpf{k}{:});
        t = sprintf('CPF - %s : %s : ', casefile, name);

        mpopt = mpoption(mpopt0, opts{:}, 'cpf.adapt_step', 1, ...
                        'cpf.nose_tol', 1e-8, 'cpf.stop_at', 'NOSE');
        mpc = loadcase(casefile);
        mpct = mpc;
        mpct.load3p(:, 4:6) = mpct.load3p(:, 4:6) * 1.2;
        cpf = run_cpf({mpc, mpct}, mpopt, 'mpx', mp.xt_3p());
        t_is(cpf.success, 1, 12, [t 'success']);
        t_is(cpf.mm.soln.output.max_lam, lam(c), 8, [t 'max_lam']);
        t_is(abs(cpf.mm.soln.output.V(1, end)), vm1(c), 8, [t '|v(1)|']);
        t_is(abs(cpf.mm.soln.output.V(end, end)), vm_end(c), 8, [t '|v(end)|']);
    end
end

%%-----  Optimal Power Flow  -----
for k = 1:length(cfg_opf)
    casefile = casefiles{7};
    [alg, name, check, opts] = deal(cfg_opf{k}{:});
    t = sprintf('OPF - %s : %s : ', casefile, name);
    mpopt = mpoption(mpopt0, 'opf.ac.solver', alg, opts{:});
    opf = run_opf(casefile, mpopt, 'mpx', {mp.xt_3p()});

    t_is(opf.success, 1, 12, [t 'success']);

    va = cell2mat(cellfun(@(x)opf.dm.elements.bus3p.tab.(x), va_fields, 'UniformOutput', 0));
    vm = cell2mat(cellfun(@(x)opf.dm.elements.bus3p.tab.(x), vm_fields, 'UniformOutput', 0)) .* ...
        (opf.dm.elements.bus3p.tab.base_kv*ones(1,3))/sqrt(3);
    pl = cell2mat(cellfun(@(x)opf.dm.elements.line3p.tab.(x), pl_fields, 'UniformOutput', 0));
    t_is(va, [eva; eva+3.656255521718970; eva+1.039772166853819], 7, [t 'va']);
    t_is(vm, [evm; evm; evm], 7, [t 'vm']);
    t_is(pl, [epl; epl; epl], 5, [t 'pl']);
end

warning('on', 'update_z:multiple_nodes');
warning(s.state, sing_mat_warn_id);

if nargout  %% set output arg
    task = pf;
end

t_end;
