function t_convert_1p_to_3p(quiet)
% t_convert_1p_to_3p - Tests of mp.case_utils.convert_1p_to_3p.
%
% Tests the conversion for many of the cases included with |MATPOWER|,
% ensuring that the power flow results of the original single-phase case
% and the resulting balanced, three-phase case match.

%   MATPOWER
%   Copyright (c) 2019-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 1
    quiet = 0;
end

cases = {'4_dist','4gs','5','6ww','9','10ba','12da','15da','15nbr','16ci',...
         '17me','18','18nbr','22','24_ieee_rts','28da','30','_ieee30', ...
         '33bw','33mg','34sa','38si','39','51ga','51he','59','60nordic','69',...
         '70da','_RTS_GMLC','74ds','85','94pi','118','118zh','136ma','141',...
         '145','_ACTIVSg200','300','_ACTIVSg500','533mt_hi','533mt_lo',...
         '1197', '_ACTIVSg2000','3012wp','3120sp'};

ntsucces = 1;
ntbus = 4;
ntgen = 4;
ntbranch = 8;
ntxfmr = 8;
ntshunt = 4;
ntbase = 9;
ntsave = 1;
warnings_off = quiet;
warnings_off = true;
diff_tool = 'bbdiff';
show_diff_on_fail = false;

nt = ntsucces+ntbus+ntgen+ntbranch+ntxfmr+ntshunt;

t_begin(nt*length(cases)+ntbase+ntsave, quiet);

if warnings_off     %% disable warning messages
    warn_id1 = 'MATPOWER:mp_case_utils_reactive_only_loads';
    warn_id2 = 'MATPOWER:mp_case_utils_remove_gen_q_lims';
    warn_id3 = 'MATPOWER:mp_case_utils_relocate_branch_shunts';
    s1 = warning('query', warn_id1);
    s2 = warning('query', warn_id2);
    s3 = warning('query', warn_id3);
    warning('off', warn_id1);
    warning('off', warn_id2);
    warning('off', warn_id3);
end
if have_feature('octave')
    if have_feature('octave', 'vnum') >= 4
        file_in_path_warn_id = 'Octave:data-file-in-path';
    else
        file_in_path_warn_id = 'Octave:load-file-in-path';
    end
    s4 = warning('query', file_in_path_warn_id);
    warning('off', file_in_path_warn_id);
end

mpopt = mpoption('verbose',0,'out.all',0,'pf.tol',1e-10,'pf.nr.max_it',20);
for c = 1:length(cases)
    skip_xfmr_flag = 0;
    skip_shunt_flag = 0;
    C = ['case' cases{c}];
    mpc = loadcase(C);
    mpc3p = mp.case_utils.convert_1p_to_3p(mpc);
    mpc = mp.case_utils.remove_gen_q_lims(mpc, C);
    mpc = mp.case_utils.to_consecutive_bus_numbers(mpc);
    mpc = mp.case_utils.relocate_branch_shunts(mpc);

    res1p = run_pf(mpc,mpopt);
    res3p = run_pf(mpc3p, mpopt, 'mpx', mp.xt_3p);
    soln3p = res3p.mm.soln;

    % 1) Check for convergence of mpc3p
    t = sprintf('3-phase %s - convergence - ', C);
    t_ok(soln3p.eflag==1, [t sprintf('%d iter', soln3p.output.iterations)])

    if soln3p.eflag
        t = sprintf('%s - 3-phase balanced results: ', C);
        % 2) Check for balanced results
        branch = res1p.dm.elements.branch.tab;
        bus3p = res3p.dm.elements.bus3p.tab;
        t_is(sum([bus3p.vm1 bus3p.vm2 bus3p.vm3],2)/3, bus3p.vm1, 6, [t 'buses (voltage magnitudes)']);
        t_is(sum(diff([bus3p.va1 bus3p.va2 bus3p.va3 bus3p.va1],1,2),2), zeros(size(bus3p,1),1), 6, [t 'buses (voltage angles)']);
        gen3p = res3p.dm.elements.gen3p.tab;
        t_is(sum([gen3p.pg1 gen3p.pg2 gen3p.pg3],2)/3, gen3p.pg1, 6, [t 'generators (active power)']);
        t_is(sum([gen3p.qg1 gen3p.qg2 gen3p.qg3],2)/3, gen3p.qg1, 6, [t 'generators (reactive power)']);
        line3p = res3p.dm.elements.line3p.tab;
        t_is(sum([line3p.pl1_fr line3p.pl2_fr line3p.pl3_fr],2)/3, line3p.pl1_fr, 6, [t 'lines (active power at from end)']);
        t_is(sum([line3p.pl1_to line3p.pl2_to line3p.pl3_to],2)/3, line3p.pl1_to, 6, [t 'lines (active power at to end)']);
        t_is(sum([line3p.ql1_fr line3p.ql2_fr line3p.ql3_fr],2)/3, line3p.ql1_fr, 6, [t 'lines (reactive power at from end)']);
        t_is(sum([line3p.ql1_to line3p.ql2_to line3p.ql3_to],2)/3, line3p.ql1_to, 6, [t 'lines (reactive power at to end)']);
        if ~isempty(mpc3p.xfmr3p)
            xfmr3p = res3p.dm.elements.xfmr3p.tab;
            t_is(sum([xfmr3p.pl1_fr xfmr3p.pl2_fr xfmr3p.pl3_fr],2)/3, xfmr3p.pl1_fr, 6, [t 'transformers (active power at from end)']);
            t_is(sum([xfmr3p.pl1_to xfmr3p.pl2_to xfmr3p.pl3_to],2)/3, xfmr3p.pl1_to, 6, [t 'transformers (active power at to end)']);
            t_is(sum([xfmr3p.ql1_fr xfmr3p.ql2_fr xfmr3p.ql3_fr],2)/3, xfmr3p.ql1_fr, 6, [t 'transformers (reactive power at from end)']);
            t_is(sum([xfmr3p.ql1_to xfmr3p.ql2_to xfmr3p.ql3_to],2)/3, xfmr3p.ql1_to, 6, [t 'transformers (reactive power at to end)']);
        else
            t_skip(ntxfmr, [t '3-phase case does not have transformers']);
            skip_xfmr_flag = 1;
        end
        if ~isempty(mpc3p.shunt3p)
            shunt3p = res3p.dm.elements.shunt3p.tab;
            t_is(sum([shunt3p.p1 shunt3p.p2 shunt3p.p3],2)/3, shunt3p.p1, 6, [t 'shunts (active power consumption)']);
            t_is(sum([shunt3p.q1 shunt3p.q2 shunt3p.q3],2)/3, shunt3p.q1, 6, [t 'shunts (reactive power consumption)']);
        else
            t_skip(ntshunt, [t '3-phase case does not have shunt elements']);
            skip_shunt_flag = 1;
        end


        t = sprintf('%s - comparison against 1-phase: ', C);
        % 3) Check for equivalent results with respect to single-phase case
        bus = res1p.dm.elements.bus.tab;
        t_is(bus3p.vm1, bus.vm, 6, [t 'buses (voltage magnitudes)'])
        t_is(bus3p.va1, bus.va, 6, [t 'buses (voltage angles)'])
        gen = res1p.dm.elements.gen.tab;
        t_is(3*gen3p.pg1/1e3, gen.pg, 6, [t 'generators (active power)']);
        t_is(3*gen3p.qg1/1e3, gen.qg, 6, [t 'generators (reactive power)']);
        id_line1p = find(branch.tm == 0);
        t_is(3*line3p.pl1_fr/1e3, branch.pl_fr(id_line1p), 6, [t 'lines (active power at from end)']);
        t_is(3*line3p.pl1_to/1e3, branch.pl_to(id_line1p), 6, [t 'lines (active power at to end)']);
        t_is(3*line3p.ql1_fr/1e3, branch.ql_fr(id_line1p), 6, [t 'lines (reactive power at from end)']);
        t_is(3*line3p.ql1_to/1e3, branch.ql_to(id_line1p), 6, [t 'lines (reactive power at to end)']);
        if ~skip_xfmr_flag
            id_xfmr1p = find(branch.tm ~= 0);
            t_is(3*xfmr3p.pl1_fr/1e3, branch.pl_fr(id_xfmr1p), 6, [t 'transformers (active power at from end)']);
            t_is(3*xfmr3p.pl1_to/1e3, branch.pl_to(id_xfmr1p), 6, [t 'transformers (active power at to end)']);
            t_is(3*xfmr3p.ql1_fr/1e3, branch.ql_fr(id_xfmr1p), 6, [t 'transformers (reactive power at from end)']);
            t_is(3*xfmr3p.ql1_to/1e3, branch.ql_to(id_xfmr1p), 6, [t 'transformers (reactive power at to end)']);
        end
        if ~skip_shunt_flag
            shunt = res1p.dm.elements.shunt.tab;
            t_is(3*shunt3p.p1/1e3, shunt.p, 6, [t 'shunts (active power consumption)']);
            t_is(3*shunt3p.q1/1e3, shunt.q, 6, [t 'shunts (reactive power consumption)']);
        end
    else
        t_skip(nt, [t '3-phase case did not converge!']);
    end
end

%% 4) Check for equal results with different bases
C = 'case141';
mpc = loadcase(C);
basekV_old = mpc.bus(:, 10);
basekV_new = basekV_old;
vmin = 0.98; vmax = 1.02;
for i = unique(basekV_old)'
    idvm = find(basekV_old==i);
    basekV_new(idvm) = i * (vmin+(vmax-vmin)*rand);
end
basekVA_new = 1000 * mpc.baseMVA * 1.5;
mpc3p_a = mp.case_utils.convert_1p_to_3p(mpc);
mpc3p_b = mp.case_utils.convert_1p_to_3p(mpc, basekVA_new, basekV_new);
mpc = mp.case_utils.remove_gen_q_lims(mpc, C);
mpc = mp.case_utils.to_consecutive_bus_numbers(mpc);
mpc = mp.case_utils.relocate_branch_shunts(mpc);

mpopt = mpoption('pf.tol', 1e-10, 'pf.nr.max_it', 1000, 'verbose', 0,'out.all', 0);
res3p_a = run_pf(mpc3p_a,mpopt,'mpx',mp.xt_3p);
res3p_b = run_pf(mpc3p_b,mpopt,'mpx',mp.xt_3p);

bus3p_a = res3p_a.dm.elements.bus3p.tab;
bus3p_b = res3p_b.dm.elements.bus3p.tab;
gen3p_a = res3p_a.dm.elements.gen3p.tab;
gen3p_b = res3p_b.dm.elements.gen3p.tab;
line3p_a = res3p_a.dm.elements.line3p.tab;
line3p_b = res3p_b.dm.elements.line3p.tab;

t = sprintf('%s - change of bases: ', C);
eflag = res3p_a.mm.soln.eflag == 1 & res3p_b.mm.soln.eflag == 1;
t_ok(eflag, [t 'convergence']);
t_is(bus3p_a.vm1.*basekV_old, bus3p_b.vm1.*basekV_new, 6, [t 'voltage magnitudes'])
t_is([bus3p_a.va1 bus3p_a.va2 bus3p_a.va3], [bus3p_b.va1 bus3p_b.va2 bus3p_b.va3], 6, [t 'voltage angles'])
t_is([gen3p_a.pg1 gen3p_a.pg2 gen3p_a.pg3], [gen3p_b.pg1 gen3p_b.pg2 gen3p_b.pg3], 6, [t 'active power generation'])
t_is([gen3p_a.qg1 gen3p_a.qg2 gen3p_a.qg3], [gen3p_b.qg1 gen3p_b.qg2 gen3p_b.qg3], 6, [t 'reactive power generation'])
t_is([line3p_a.pl1_fr line3p_a.pl2_fr line3p_a.pl3_fr], [line3p_b.pl1_fr line3p_b.pl2_fr line3p_b.pl3_fr], 5, [t 'line active power injection at from end'])
t_is([line3p_a.ql1_fr line3p_a.ql2_fr line3p_a.ql3_fr], [line3p_b.ql1_fr line3p_b.ql2_fr line3p_b.ql3_fr], 5, [t 'line reactive power injection at from end'])
t_is([line3p_a.pl1_to line3p_a.pl2_to line3p_a.pl3_to], [line3p_b.pl1_to line3p_b.pl2_to line3p_b.pl3_to], 5, [t 'line active power injection at to end'])
t_is([line3p_a.ql1_to line3p_a.ql2_to line3p_a.ql3_to], [line3p_b.ql1_to line3p_b.ql2_to line3p_b.ql3_to], 5, [t 'line reactive power injection at to end'])

t = 'save converted case: case24_ieee_rts';
C = 'case24_ieee_rts';
fname_e = 't_case24_ieee_rts_3p.m';
fname_g = sprintf('t_case24_ieee_rts_3p_%d.m', fix(1e9*rand));
mp.case_utils.convert_1p_to_3p(fname_g, C);
% equivalent to:
% savecase(fname_g, mp.case_utils.convert_1p_to_3p(C));
reps = {{'_(\d{6,9})', '', 1}};
if ~t_file_match(fname_g, fname_e, t, reps, 1);
    fprintf('  compare these 2 files:\n    %s\n    %s\n', fname_g, fname_e);
    if show_diff_on_fail
        cmd = sprintf('%s %s %s', diff_tool, fname_g, fname_e);
        [status, result] = system(cmd);
        keyboard
    end
end

if warnings_off     %% re-enable warning messages
    warning(s1.state, warn_id1);
    warning(s2.state, warn_id2);
    warning(s3.state, warn_id3);
end
if have_feature('octave')
    warning(s4.state, file_in_path_warn_id);
end

t_end;
