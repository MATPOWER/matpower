function t_convert_1p_to_3p(quiet)
% t_mpc2mpc3p - test of convert_1p_to_3p transformation function over many
%   of the cases included in the MATPOWER suite.

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
         '33bw','33mg','34sa','38si','39','51ga','51he','60nordic','69',...
         '70da','_RTS_GMLC','74ds','85','94pi','118','118zh','136ma','141',...
         '145','_ACTIVSg200','300','_ACTIVSg500','533mt_hi','533mt_lo',...
         '_ACTIVSg2000','3012wp','3120sp'};

ntsucces = 1;
ntbus = 4;
ntgen = 4;
ntbranch = 8;
ntxfmr = 8;
ntshunt = 4;
ntbase = 9;

nt = ntsucces+ntbus+ntgen+ntbranch+ntxfmr+ntshunt;

t_begin(nt*length(cases)+ntbase, quiet);

mpopt = mpoption('verbose',0,'out.all',0,'pf.tol',1e-10,'pf.nr.max_it',20);
for c = 1:length(cases)
    skip_xfmr_flag = 0;
    skip_shunt_flag = 0;
    C = ['case' cases{c}];
    mpc = ext2int(loadcase(C));
    mpc = check_mpc(mpc, C);
    mpc3p = convert_1p_to_3p(mpc);

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


        t = sprintf('%s - comparison againts 1-phase: ', C);
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
mpc = check_mpc(loadcase(C));
basekV_old = mpc.bus(:, 10);
basekV_new = basekV_old;
vmin = 0.98; vmax = 1.02;
for i = unique(basekV_old)'
    idvm = find(basekV_old==i);
    basekV_new(idvm) = i * (vmin+(vmax-vmin)*rand);
end
basekVA_new = 1000 * mpc.baseMVA * 1.5;
mpc3p_a = convert_1p_to_3p(mpc);
mpc3p_b = convert_1p_to_3p(mpc, basekVA_new, basekV_new);

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

t_end;

end

function mpc2 = check_mpc(mpc, case_name)
    [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
        VA, BASE_KV, ZONE, VMAX, VMIN] = idx_bus;
    [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
        TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
        ANGMIN, ANGMAX] = idx_brch;
    [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
        MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
        QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

    mpc2 = mpc;

    %% 1) Look for several generators connected at the same bus and ... (applies for case5 and case24_ieee_rts)
    id_bus_gen = mpc.gen(:,GEN_BUS);
    unique_id_bus_gen = unique(id_bus_gen);
    if length(unique_id_bus_gen) < length(id_bus_gen)
        for b = unique_id_bus_gen'
            id_b = find(id_bus_gen==b);
            if length(id_b) > 1     %% ... eliminate their reactive limits to avoid differences in reactive power allocation
                mpc2.gen(id_b, QMAX) = Inf(length(id_b),1);
                mpc2.gen(id_b, QMIN) = -Inf(length(id_b),1);
            end
        end
        warning(['t_convert_1p_to_3p: %s: Removing reactive power limits of generators connected to the same bus. This gens are: (%s).'], case_name, strjoin(cellstr(num2str(unique_id_bus_gen(:))), ', '));
    end

    %% 2) Look for general branches and ... ()
    general_branch_nom = mpc2.branch(:,BR_B) ~= 0 & mpc2.branch(:,TAP) ~= 0;
    if any(general_branch_nom)
        id_general_branch_nom = find(general_branch_nom);
        for b = id_general_branch_nom'
            branch = mpc2.branch(b,:);
            fbus = branch(F_BUS);
            tbus = branch(T_BUS);
            mpc2.bus(tbus,BS) = mpc2.bus(tbus,BS) + mpc2.baseMVA*(branch(BR_B)/2);                      % ... move right shunt element to receiving bus and ...
            mpc2.bus(fbus,BS) = mpc2.bus(fbus,BS) + mpc2.baseMVA*(branch(BR_B)/2)*(1/(branch(TAP)^2));  % ... move left shunt element to primary side of ideal transformer and ...
            branch(BR_B) = 0;                                                                           % ... zero-out branch shunt elements.
            mpc2.branch(b,:) = branch;
        end
        warning('t_convert_1p_to_3p: %s: Relocating branch shunt elements from pi-circuit model to sending/receiving buses in the following branches: %s', case_name, strjoin(cellstr(num2str(id_general_branch_nom(:))), ', '));
    end
end