function [mpc, success, iterations] = radial_pf(mpc,mpopt)
% radial_pf - Solves the power flow using a backward-forward sweep method.
% ::
%
%   [mpc, success, iterations] = radial_pf(mpc,mpopt)
%
%   Inputs:
%       mpc : MATPOWER case struct with internal bus numbering
%       mpopt : MATPOWER options struct to override default options
%           can be used to specify the solution algorithm, output options
%           termination tolerances, and more.
%
%   Outputs:
%       mpc : results struct with all fields from the input MATPOWER case,
%             with solved voltages, active and reactive power flows
%             and generator active and reactive power output.
%       success : success flag, 1 = succeeded, 0 = failed
%       iterations : number of iterations
%
% See also caseformat, loadcase, mpoption.

%   MATPOWER
%   Copyright (c) 2019-2024, Power Systems Engineering Research Center (PSERC)
%   by Mirko Todorovski
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% define named indices into bus, gen, branch matrices
define_constants;
%% branch ordering
mpc = order_radial(mpc);
if ~isempty(mpc.loop)
    fprintf('\nNumber of detected loops: %i\n', length(mpc.loop));
    if mpopt.verbose > 0
        fprintf('\nBranches forming loops\n')
        fprintf('LOOP# F_BUS T_BUS\n');
        for i = 1:length(mpc.loop)
            fprintf('%5i %5i %5i\n',i,mpc.order.bus.i2e(mpc.branch(mpc.loop(i),[F_BUS T_BUS])));
        end
    end
    error('radial_pf: power flow algorithm %s can only handle radial networks.', mpopt.pf.alg)
end
%% define vectors needed for backward-forward sweep method
% branch data
[f, t, Zb, Yb, Ysh] = ...
    deal(mpc.branch(:,F_BUS),mpc.branch(:,T_BUS), ...
         mpc.branch(:,BR_R)+1j*mpc.branch(:,BR_X),1j*mpc.branch(:,BR_B), ...
         mpc.bus(:,GS)+1j*mpc.bus(:,BS));
nl = size(mpc.branch,1);
nb = size(mpc.bus,1);
Ysh = Ysh/mpc.baseMVA;
tap = ones(nl, 1);        % default tap ratio = 1
i = find(mpc.branch(:, TAP)); % indices of non-zero tap ratios
tap(i) = mpc.branch(i, TAP);  % assign non-zero tap ratios
Ybf = Yb/2 + 1./tap .* (1./tap-1) ./ Zb;
Ybt = Yb/2 + (1-1./tap) ./ Zb;
[Ybf(mpc.br_reverse), Ybt(mpc.br_reverse)] = ...
    deal(Ybt(mpc.br_reverse), Ybf(mpc.br_reverse)); % reversed branches
Zb = Zb .* tap;
Yd = Ysh + (sparse(f, f, Ybf, nb, nb) + sparse(t, t, Ybt, nb, nb)) * ones(nb,1);
% vector of complex bus power injections
Sbus = makeSbus(mpc.baseMVA, mpc.bus, mpc.gen);
% generator data (other than the slack bus)
[ref, pv, pq] = bustypes(mpc.bus, mpc.gen);
on = find(mpc.gen(:, GEN_STATUS) > 0);  % which generators are on?
gbus = mpc.gen(on, GEN_BUS);  % what buses are they at?
vcb = zeros(nb, 1);  % mask of voltage-controlled buses
vcb(pv) = 1;  % include PV buses
k = find(vcb(gbus));  % in-service gens at v-c buses
Vg = ones(nb, 1);
Vg(gbus(k)) = mpc.gen(on(k), VG);  % vector of gen set voltages
Vg = Vg(pv);
Pg = real(Sbus(pv))*mpc.baseMVA + mpc.bus(pv,PD);
% load data
Sd = zeros(nb, 1);
Sd(pq) = -Sbus(pq);
Sd(pv) = (mpc.bus(pv,PD) + 1j*mpc.bus(pv, QD))/mpc.baseMVA;
%% calculate voltages and power flows
Vslack = mpc.gen(1,VG);
switch upper(mpopt.pf.alg)
    case 'PQSUM'
        [V, Qpv, Sf, St, Sslack, iterations, success] = calc_v_pq_sum(Vslack,nb,nl,f,Zb,Ybf,Ybt,Yd,Sd,pv,Pg,Vg,mpopt);
    case 'ISUM'
        [V, Qpv, Sf, St, Sslack, iterations, success] = calc_v_i_sum(Vslack,nb,nl,f,Zb,Ybf,Ybt,Yd,Sd,pv,Pg,Vg,mpopt);
    case 'YSUM'
        [V, Qpv, Sf, St, Sslack, iterations, success] = calc_v_y_sum(Vslack,nb,nl,f,Zb,Ybf,Ybt,Yd,Sd,pv,Pg,Vg,mpopt);
end
mpc.success = success;
%% update data matrices with solution
mpc.bus(:,VM) = abs(V);
mpc.bus(:,VA) = angle(V)/pi*180;
mpc.branch(:,PF) = real(Sf)*mpc.baseMVA;
mpc.branch(:,QF) = imag(Sf)*mpc.baseMVA;
mpc.branch(:,PT) = -real(St)*mpc.baseMVA;
mpc.branch(:,QT) = -imag(St)*mpc.baseMVA;
[mpc.branch(mpc.br_reverse,PF), mpc.branch(mpc.br_reverse,PT)] = ...
    deal(mpc.branch(mpc.br_reverse,PT), mpc.branch(mpc.br_reverse,PF));
[mpc.branch(mpc.br_reverse,QF), mpc.branch(mpc.br_reverse,QT)] = ...
    deal(mpc.branch(mpc.br_reverse,QT), mpc.branch(mpc.br_reverse,QF));
mpc.gen(1,PG) = real(Sslack)*mpc.baseMVA;
mpc.gen(1,QG) = imag(Sslack)*mpc.baseMVA;
if ~isempty(pv)
    % C = gbus == pv';  % map gbus to pv (does not work on MATLAB <R2016b)
    C = gbus*ones(1,length(pv)) == ones(length(gbus),1)*pv';
    Q = C*Qpv;  % Q for generators at gbus
    k = find(Q);  % update QG in gen for buses with non-zero Qpv
    mpc.gen(k,QG) = Q(k)*mpc.baseMVA;
end
% At this point any buses with more than one generator will have
% the total Q dispatch for the bus assigned to each generator. This
% must be split between them using pfsoln.
temp = mpc;
% Check whether branches need to be reversed before creating Ybus
k = mpc.br_reverse;
if any(k)
    % We need this in case there are transformers
    temp.branch(k,BR_R) = tap(k).^2.*temp.branch(k,BR_R);
    temp.branch(k,BR_X) = tap(k).^2.*temp.branch(k,BR_X);
    temp.branch(k,TAP) = 1./tap(k);
end
[Ybus, Yf, Yt] = makeYbus(temp);
baseMVA = mpc.baseMVA; bus0 = mpc.bus; gen0 = mpc.gen; branch0 = mpc.branch;
[~, gen] = pfsoln(baseMVA, bus0, gen0, branch0, Ybus, Yf, Yt, V, ref, pv, pq);
mpc.gen = gen;
%% reverse bus and branch ordering
mpc.bus = mpc.bus(mpc.bus_order_inv,:);
mpc.bus(:,BUS_I) = mpc.bus_order(mpc.bus(:,BUS_I));
[f, t] = deal(mpc.branch(:,F_BUS),mpc.branch(:,T_BUS));
[f(mpc.br_reverse), t(mpc.br_reverse)] =  deal(t(mpc.br_reverse), f(mpc.br_reverse));
mpc.branch(:,[F_BUS T_BUS]) = [mpc.bus_order(f) mpc.bus_order(t)];
mpc.branch = mpc.branch(mpc.branch_order_inv,:);
mpc.gen(:,GEN_BUS) = mpc.bus_order(mpc.gen(:,GEN_BUS));
