function printpf(baseMVA, bus, gen, branch, f, success, et, fd, mpopt)
%PRINTPF   Prints power flow results.
%   printpf(baseMVA, bus, gen, branch, f, success, et, fd, mpopt) prints
%   powerflow results to fd (a file descriptor which defaults to STDOUT).
%   mpopt is a MATPOWER options vector (see 'help mpoption' for details).
%   Uses default options if this parameter is not given. The objective
%   function value is given in f and the elapsed time (seconds to compute
%   opf) in et.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2003 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%%----- initialization -----
%% default arguments
if nargin < 9
	mpopt = mpoption;	%% use default options
	if nargin < 8
		fd = 1;			%% print to stdio by default
	end
end
if isempty(f)
	isOPF = 0;		%% have only simple PF data
else
	isOPF = 1;		%% have OPF data
end

%% options
dc				= mpopt(10);		%% use DC formulation?
OUT_ALL			= mpopt(32);
OUT_ANY			= OUT_ALL == 1;		%% set to true if any pretty output is to be generated
OUT_SYS_SUM		= OUT_ALL == 1 | (OUT_ALL == -1 & mpopt(33));
OUT_AREA_SUM	= OUT_ALL == 1 | (OUT_ALL == -1 & mpopt(34));
OUT_BUS			= OUT_ALL == 1 | (OUT_ALL == -1 & mpopt(35));
OUT_BRANCH		= OUT_ALL == 1 | (OUT_ALL == -1 & mpopt(36));
OUT_GEN			= OUT_ALL == 1 | (OUT_ALL == -1 & mpopt(37));
OUT_ANY			= OUT_ANY | (OUT_ALL == -1 & (OUT_SYS_SUM | OUT_AREA_SUM | OUT_BUS | OUT_BRANCH | OUT_GEN));
if OUT_ALL == -1
	OUT_ALL_LIM	= mpopt(38);
elseif OUT_ALL == 1
	OUT_ALL_LIM	= 2;
else
	OUT_ALL_LIM = 0;
end
OUT_ANY			= OUT_ANY | OUT_ALL_LIM >= 1;
if OUT_ALL_LIM == -1
	OUT_V_LIM		= mpopt(39);
	OUT_LINE_LIM	= mpopt(40);
	OUT_PG_LIM		= mpopt(41);
	OUT_QG_LIM		= mpopt(42);
else
	OUT_V_LIM		= OUT_ALL_LIM;
	OUT_LINE_LIM	= OUT_ALL_LIM;
	OUT_PG_LIM		= OUT_ALL_LIM;
	OUT_QG_LIM		= OUT_ALL_LIM;
end
OUT_ANY			= OUT_ANY | (OUT_ALL_LIM == -1 & (OUT_V_LIM | OUT_LINE_LIM | OUT_PG_LIM | OUT_QG_LIM));
OUT_RAW			= mpopt(43);

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
	VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, ...
	GEN_STATUS, PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
	RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST] = idx_brch;

%% constant
j = sqrt(-1);

%% internal bus number
e2i = zeros(max(bus(:, BUS_I)), 1);		%% need internal bus numbering for a second
e2i(bus(:, BUS_I)) = [1:size(bus, 1)]';

%% sizes of things
nb = size(bus, 1);		%% number of buses
nl = size(branch, 1);	%% number of branches
ng = size(gen, 1);		%% number of generators

%% zero out some data to make printout consistent for DC case
if dc
	bus(:, [QD, BS])			= zeros(nb, 2);
	gen(:, [VG, QMAX, QMIN])	= zeros(ng, 3);
	branch(:, [BR_R, BR_B])		= zeros(nl, 2);
end

%% parameters
ties = find(bus(e2i(branch(:, F_BUS)), BUS_AREA) ~= bus(e2i(branch(:, T_BUS)), BUS_AREA));
						%% area inter-ties
tap = ones(nl, 1);								%% default tap ratio = 1 for lines
xfmr = find(branch(:, TAP));					%% indices of transformers
tap(xfmr) = branch(xfmr, TAP);					%% include transformer tap ratios
tap = tap .* exp(-j*pi/180 * branch(:, SHIFT));	%% add phase shifters
nzld = find(bus(:, PD) | bus(:, QD));
sorted_areas = sort(bus(:, BUS_AREA));
s_areas = sorted_areas([1; find(diff(sorted_areas))+1]);	%% area numbers
na = length(s_areas);							%% number of areas
on = find(gen(:, GEN_STATUS) > 0);
V = bus(:, VM) .* exp(sqrt(-1) * pi/180 * bus(:, VA));
out = find(branch(:, BR_STATUS) == 0);			%% out-of-service branches
nout = length(out);
if dc
	loss = zeros(nl,1);
else
	loss = baseMVA * abs(V(e2i(branch(:, F_BUS))) ./ tap - V(e2i(branch(:, T_BUS)))) .^ 2 ./ ...
				(branch(:, BR_R) - j * branch(:, BR_X));
end
fchg = abs(V(e2i(branch(:, F_BUS))) ./ tap) .^ 2 .* branch(:, BR_B) * baseMVA / 2;
tchg = abs(V(e2i(branch(:, T_BUS)))       ) .^ 2 .* branch(:, BR_B) * baseMVA / 2;
loss(out) = zeros(nout, 1);
fchg(out) = zeros(nout, 1);
tchg(out) = zeros(nout, 1);

%%----- print the stuff -----
if OUT_ANY
	%% convergence & elapsed time
	if success
		fprintf(fd, '\nConverged in %.2f seconds', et);
	else
		fprintf(fd, '\nDid not converge (%.2f seconds)\n', et);
	end
	
	%% objective function value
	if isOPF
		fprintf(fd, '\nObjective Function Value = %.2f $/hr', f);
	end
end
if OUT_SYS_SUM
	fprintf(fd, '\n================================================================================');
	fprintf(fd, '\n|     System Summary                                                           |');
	fprintf(fd, '\n================================================================================');
	fprintf(fd, '\n\nHow many?                  How much?              P (MW)       Q (MVAr)');
	fprintf(fd, '\n---------------------      --------------------  --------  -----------------');
	fprintf(fd, '\nBuses          %5d       Total Gen Capacity    %7.1f  %7.1f to +%.1f', nb, sum(gen(:, PMAX)), sum(gen(:, QMIN)), sum(gen(:, QMAX)));
	fprintf(fd, '\nGenerators     %5d       On-line Capacity      %7.1f  %7.1f to +%.1f', ng, sum(gen(on, PMAX)), sum(gen(on, QMIN)), sum(gen(on, QMAX)));
	fprintf(fd, '\nCommited Gens  %5d       Generation (current)  %7.1f      %7.1f', length(on), sum(gen(on, PG)), sum(gen(on, QG)));
	fprintf(fd, '\nLoads          %5d       Load                  %7.1f      %7.1f', length(nzld), sum(bus(nzld, PD)), sum(bus(nzld, QD)));
	fprintf(fd, '\nBranches       %5d       Losses (I^2 * Z)      %8.2f     %8.2f', nl, sum(real(loss)), sum(imag(loss)) );
	fprintf(fd, '\nTransformers   %5d       Branch Charging (inj)      -       %7.1f', length(xfmr), sum(fchg) + sum(tchg) );
	fprintf(fd, '\nAreas          %5d       Shunt (inj)           %7.1f      %7.1f', length(s_areas), ...
		-sum(bus(:, VM) .^ 2 .* bus(:, GS)), sum(bus(:, VM) .^ 2 .* bus(:, BS))	);
	fprintf(fd, '\nInter-ties     %5d       Total Inter-tie Flow  %7.1f      %7.1f', length(ties), sum(abs(branch(ties, PF)-branch(ties, PT))) / 2, sum(abs(branch(ties, QF)-branch(ties, QT))) / 2);
	fprintf(fd, '\n');
	fprintf(fd, '\n                          Minimum                      Maximum');
	fprintf(fd, '\n                 -------------------------  --------------------------------');
	[minv, mini] = min(bus(:, VM));
	[maxv, maxi] = max(bus(:, VM));
	fprintf(fd, '\nVoltage Magnitude %7.3f p.u. @ bus %-4d     %7.3f p.u. @ bus %-4d', minv, bus(mini, BUS_I), maxv, bus(maxi, BUS_I));
	[minv, mini] = min(bus(:, VA));
	[maxv, maxi] = max(bus(:, VA));
	fprintf(fd, '\nVoltage Angle   %8.2f deg   @ bus %-4d   %8.2f deg   @ bus %-4d', minv, bus(mini, BUS_I), maxv, bus(maxi, BUS_I));
	[maxv, maxi] = max(real(loss));
	fprintf(fd, '\nP Losses (I^2*R)             -              %8.2f MW    @ line %d-%d', maxv, branch(maxi, F_BUS), branch(maxi, T_BUS));
	[maxv, maxi] = max(imag(loss));
	fprintf(fd, '\nQ Losses (I^2*X)             -              %8.2f MVAr  @ line %d-%d', maxv, branch(maxi, F_BUS), branch(maxi, T_BUS));
	if isOPF
		[minv, mini] = min(bus(:, LAM_P));
		[maxv, maxi] = max(bus(:, LAM_P));
		fprintf(fd, '\nLambda P        %8.2f $/MWh @ bus %-4d   %8.2f $/MWh @ bus %-4d', minv, bus(mini, BUS_I), maxv, bus(maxi, BUS_I));
		[minv, mini] = min(bus(:, LAM_Q));
		[maxv, maxi] = max(bus(:, LAM_Q));
		fprintf(fd, '\nLambda Q        %8.2f $/MWh @ bus %-4d   %8.2f $/MWh @ bus %-4d', minv, bus(mini, BUS_I), maxv, bus(maxi, BUS_I));
	end
	fprintf(fd, '\n');
end

if OUT_AREA_SUM
	fprintf(fd, '\n================================================================================');
	fprintf(fd, '\n|     Area Summary                                                             |');
	fprintf(fd, '\n================================================================================');
	fprintf(fd, '\nArea  # of   # of   Gens   # of   # of   # of   # of    Total Gen    On-line');
	fprintf(fd, '\n Num  Buses  Gens  Online  Loads  Brchs  Xfmrs  Ties     Capacity    Capacity');
	fprintf(fd, '\n----  -----  ----  ------  -----  -----  -----  ----   -----------  -----------');
	for i=1:length(s_areas)
		a = s_areas(i);
		ib = find(bus(:, BUS_AREA) == a);
		ig = find(bus(e2i(gen(:, GEN_BUS)), BUS_AREA) == a);
		igon = find(bus(e2i(gen(:, GEN_BUS)), BUS_AREA) == a & gen(:, GEN_STATUS) > 0);
		inzld = find(bus(:, BUS_AREA) == a & (bus(:, PD) | bus(:, QD)));
		ibrch = find(bus(e2i(branch(:, F_BUS)), BUS_AREA) == a & bus(e2i(branch(:, T_BUS)), BUS_AREA) == a);
		in_tie = find(bus(e2i(branch(:, F_BUS)), BUS_AREA) == a & bus(e2i(branch(:, T_BUS)), BUS_AREA) ~= a);
		out_tie = find(bus(e2i(branch(:, F_BUS)), BUS_AREA) ~= a & bus(e2i(branch(:, T_BUS)), BUS_AREA) == a);
		if length(xfmr)
			nxfmr = length(find(bus(e2i(branch(xfmr, F_BUS)), BUS_AREA) == a & bus(e2i(branch(xfmr, T_BUS)), BUS_AREA) == a));
		else
			nxfmr = 0;
		end
		fprintf(fd, '\n%3d   %5d %5d  %5d  %5d  %5d  %5d %5d    %7.1f MW   %7.1f MW', ...
			a, length(ib), length(ig), length(igon), length(inzld), length(ibrch), nxfmr, ...
			length(in_tie)+length(out_tie), sum(gen(ig, PMAX)), sum(gen(igon, PMAX)));
	end
	fprintf(fd, '\n----  -----  ----  ------  -----  -----  -----  ----   -----------  -----------');
	fprintf(fd, '\nTot:  %5d %5d  %5d  %5d  %5d  %5d %5d    %7.1f MW   %7.1f MW', ...
		nb, ng, length(on), length(nzld), nl, length(xfmr), ...
		length(ties), sum(gen(:, PMAX)), sum(gen(on, PMAX)));
	fprintf(fd, '\n');
	fprintf(fd, '\nArea   Generation       Load          Losses       Net Export    Brnch   Shunt');
	fprintf(fd, '\n Num   MW    MVAr     MW    MVAr     MW    MVAr    MW     MVAr   Chrgng   MVAr');
	fprintf(fd, '\n---- ------ ------  ------ ------  ------ ------  ------ ------  ------  ------');
	for i=1:length(s_areas)
		a = s_areas(i);
		ib = find(bus(:, BUS_AREA) == a);
		ig = find(bus(e2i(gen(:, GEN_BUS)), BUS_AREA) == a);
		igon = find(bus(e2i(gen(:, GEN_BUS)), BUS_AREA) == a & gen(:, GEN_STATUS) > 0);
		inzld = find(bus(:, BUS_AREA) == a & (bus(:, PD) | bus(:, QD)));
		ibrch = find(bus(e2i(branch(:, F_BUS)), BUS_AREA) == a & bus(e2i(branch(:, T_BUS)), BUS_AREA) == a & branch(:, BR_STATUS));
		in_tie = find(bus(e2i(branch(:, F_BUS)), BUS_AREA) ~= a & bus(e2i(branch(:, T_BUS)), BUS_AREA) == a & branch(:, BR_STATUS));
		out_tie = find(bus(e2i(branch(:, F_BUS)), BUS_AREA) == a & bus(e2i(branch(:, T_BUS)), BUS_AREA) ~= a & branch(:, BR_STATUS));
		fprintf(fd, '\n%3d %7.1f%7.1f %7.1f%7.1f%7.2f %7.2f %7.1f%7.1f %7.1f %7.1f', ...
			a, sum(gen(igon, PG)), sum(gen(igon, QG)), sum(bus(inzld, PD)), sum(bus(inzld, QD)), ...
			sum(real(loss(ibrch))) + sum(real(loss([in_tie; out_tie]))) / 2, ...
			sum(imag(loss(ibrch))) + sum(imag(loss([in_tie; out_tie]))) / 2, ...
			sum(branch(in_tie, PT))+sum(branch(out_tie, PF)), ...
			sum(branch(in_tie, QT))+sum(branch(out_tie, QF)), ...
			sum(fchg(ibrch)) + sum(tchg(ibrch)) + sum(fchg(out_tie)) + sum(tchg(in_tie)), ...
			sum(bus(ib, VM) .^ 2 .* bus(ib, BS))	);
	end
	fprintf(fd, '\n---- ------ ------  ------ ------  ------ ------  ------ ------  ------  ------');
	fprintf(fd, '\nTot:%7.1f%7.1f %7.1f%7.1f%7.2f %7.2f %7s%7s %7.1f %7.1f', ...
		sum(gen(on, PG)), sum(gen(on, QG)), sum(bus(nzld, PD)), sum(bus(nzld, QD)), ...
		sum(real(loss)), sum(imag(loss)), '-  ', '-  ', ...
		sum(fchg) + sum(tchg), sum(bus(:, VM) .^ 2 .* bus(:, BS))	);
	fprintf(fd, '\n');
end

%% generator data
if OUT_GEN
	if isOPF
		genlamP = bus(e2i(gen(:, GEN_BUS)), LAM_P);
		genlamQ = bus(e2i(gen(:, GEN_BUS)), LAM_Q);
	end
	fprintf(fd, '\n================================================================================');
	fprintf(fd, '\n|     Generator Data                                                           |');
	fprintf(fd, '\n================================================================================');
	fprintf(fd, '\nGen  Bus     Pg        Qg   ');
	if isOPF, fprintf(fd, '   Lambda ($/MVA-hr)'); end
	fprintf(fd, '\n #    #     (MW)     (MVAr) ');
	if isOPF, fprintf(fd, '     P         Q    '); end
	fprintf(fd, '\n---  ---  --------  --------');
	if isOPF, fprintf(fd, '  --------  --------'); end
	for i = 1:ng
		if gen(i, PG) | gen(i, QG)
			fprintf(fd, '\n%3d%5d%9.2f%10.2f%10.2f%10.2f', i, gen(i, GEN_BUS), gen(i, PG), gen(i, QG));
		else
			fprintf(fd, '\n%3d%5d      -         -  ', i, gen(i, GEN_BUS));
		end
		if isOPF, fprintf(fd, '%10.2f%10.2f', genlamP(i), genlamQ(i)); end
	end
	fprintf(fd, '\n          --------  --------');
	fprintf(fd, '\n  Total:%9.2f%10.2f', sum(gen(:, PG)), sum(gen(:, QG)));
	fprintf(fd, '\n');
end
		
%% bus data
if OUT_BUS
	fprintf(fd, '\n================================================================================');
	fprintf(fd, '\n|     Bus Data                                                                 |');
	fprintf(fd, '\n================================================================================');
	fprintf(fd, '\nBus      Voltage           Generation             Load        ');
	if isOPF, fprintf(fd, '  Lambda($/MVA-hr)'); end
	fprintf(fd, '\n #   Mag(pu)  Ang(deg)   P (MW)   Q (MVAr)   P (MW)   Q (MVAr)');
	if isOPF, fprintf(fd, '     P        Q   '); end
	fprintf(fd, '\n---  -------  --------  --------  --------  --------  --------');
	if isOPF, fprintf(fd, '  -------  -------'); end
	for i = 1:nb
		fprintf(fd, '\n%3d%8.3f%10.3f', bus(i, [BUS_I, VM, VA]));
		g = find(gen(:, GEN_BUS) == bus(i, BUS_I));
		if any(gen(g, GEN_STATUS) > 0)
			fprintf(fd, '%10.2f%10.2f', sum(gen(g, PG)), sum(gen(g, QG)));
		else
			fprintf(fd, '       -         -  ');
		end
		if bus(i, PD) | bus(i, QD)
			fprintf(fd, '%10.2f%10.2f', bus(i, [PD, QD]));
		else
			fprintf(fd, '       -         -  ');
		end
		if isOPF
			fprintf(fd, '%10.3f', bus(i, LAM_P));
			if abs(bus(i, LAM_Q)) > 1e-6
				fprintf(fd, '%8.3f', bus(i, LAM_Q));
			else
				fprintf(fd, '     -   ');
			end
		end
	end
	fprintf(fd, '\n                        --------  --------  --------  --------');
	fprintf(fd, '\n               Total: %9.2f %9.2f %9.2f %9.2f', ...
		sum(gen(:, PG)), sum(gen(:, QG)), sum(bus(nzld, PD)), sum(bus(nzld, QD)));
	fprintf(fd, '\n');
end

%% branch data
if OUT_BRANCH
	fprintf(fd, '\n================================================================================');
	fprintf(fd, '\n|     Branch Data                                                              |');
	fprintf(fd, '\n================================================================================');
	fprintf(fd, '\nBrnch   From   To    From Bus Injection   To Bus Injection     Loss (I^2 * Z)  ');
	fprintf(fd, '\n  #     Bus    Bus    P (MW)   Q (MVAr)   P (MW)   Q (MVAr)   P (MW)   Q (MVAr)');
	fprintf(fd, '\n-----  -----  -----  --------  --------  --------  --------  --------  --------');
	fprintf(fd, '\n%4d%7d%7d%10.2f%10.2f%10.2f%10.2f%10.3f%10.2f', ...
			[	[1:nl]', branch(:, [F_BUS, T_BUS]), ...
				branch(:, [PF, QF]), branch(:, [PT, QT]), ...
				real(loss), imag(loss) ...
			]');
	fprintf(fd, '\n                                                             --------  --------');
	fprintf(fd, '\n                                                    Total:%10.3f%10.2f', ...
			sum(real(loss)), sum(imag(loss)));
	fprintf(fd, '\n');
end
	
%%-----  constraint data  -----
if isOPF
	%% voltage constraints
	if OUT_V_LIM == 2 | (OUT_V_LIM == 1 & ...
						 (any(bus(:, MU_VMIN) > 1e-6) | any(bus(:, MU_VMAX) > 1e-6)))
		fprintf(fd, '\n================================================================================');
		fprintf(fd, '\n|     Voltage Constraints                                                      |');
		fprintf(fd, '\n================================================================================');
		fprintf(fd, '\nBus');
		fprintf(fd, '\n #   Vmin mu  Vmin    |V|   Vmax   Vmax mu');
		fprintf(fd, '\n---  -------  -----  -----  -----  -------');
		for i = 1:nb
			if OUT_V_LIM == 2 | (OUT_V_LIM == 1 & ...
						 (bus(i, MU_VMIN) > 1e-6 | bus(i, MU_VMAX) > 1e-6))
				fprintf(fd, '\n%3d', bus(i, BUS_I));
				if bus(i, MU_VMIN) > 1e-6
					fprintf(fd, '%8.3f', bus(i, MU_VMIN));
				else
					fprintf(fd, '     -  ');
				end
				fprintf(fd, '%8.3f%7.3f%7.3f', bus(i, [VMIN, VM, VMAX]));
				if bus(i, MU_VMAX) > 1e-6
					fprintf(fd, '%8.3f', bus(i, MU_VMAX));
				else
					fprintf(fd, '     -   ');
				end
			end
		end
		fprintf(fd, '\n');
	end
		
	%% generator P constraints
	if OUT_PG_LIM == 2 | OUT_QG_LIM == 2 | ...
			(OUT_PG_LIM == 1 & (any(gen(:, MU_PMIN) > 1e-6) | any(gen(:, MU_PMAX) > 1e-6))) | ...
			(OUT_QG_LIM == 1 & (any(gen(:, MU_QMIN) > 1e-6) | any(gen(:, MU_QMAX) > 1e-6)))
		fprintf(fd, '\n================================================================================');
		fprintf(fd, '\n|     Generation Constraints                                                   |');
		fprintf(fd, '\n================================================================================');
	end
	if OUT_PG_LIM == 2 | (OUT_PG_LIM == 1 & ...
						 (any(gen(:, MU_PMIN) > 1e-6) | any(gen(:, MU_PMAX) > 1e-6)))
		fprintf(fd, '\nGen  Bus               Active Power Limits');
		fprintf(fd, '\n #    #   Pmin mu    Pmin       P        Pmax    Pmax mu');
		fprintf(fd, '\n---  ---  -------  --------  --------  --------  -------');
		for i = 1:ng
			if OUT_PG_LIM == 2 | (OUT_PG_LIM == 1 & (gen(i, MU_PMIN) > 1e-6 | gen(i, MU_PMAX) > 1e-6))
				fprintf(fd, '\n%3d%5d', i, gen(i, GEN_BUS));
				if gen(i, MU_PMIN) > 1e-6
					fprintf(fd, '%8.3f', gen(i, MU_PMIN));
				else
					fprintf(fd, '     -  ');
				end
				if gen(i, PG)
					fprintf(fd, '%10.2f%10.2f%10.2f', gen(i, [PMIN, PG, PMAX]));
				else
					fprintf(fd, '%10.2f       -  %10.2f', gen(i, [PMIN, PMAX]));
				end
				if gen(i, MU_PMAX) > 1e-6
					fprintf(fd, '%9.3f', gen(i, MU_PMAX));
				else
					fprintf(fd, '      -  ');
				end
			end
		end
		fprintf(fd, '\n');
	end
		
	%% generator Q constraints
	if OUT_QG_LIM == 2 | (OUT_QG_LIM == 1 & ...
						 (any(gen(:, MU_QMIN) > 1e-6) | any(gen(:, MU_QMAX) > 1e-6)))
		fprintf(fd, '\nGen  Bus              Reactive Power Limits');
		fprintf(fd, '\n #    #   Qmin mu    Qmin       Q        Qmax    Qmax mu');
		fprintf(fd, '\n---  ---  -------  --------  --------  --------  -------');
		for i = 1:ng
			if OUT_QG_LIM == 2 | (OUT_QG_LIM == 1 & (gen(i, MU_QMIN) > 1e-6 | gen(i, MU_QMAX) > 1e-6))
				fprintf(fd, '\n%3d%5d', i, gen(i, GEN_BUS));
				if gen(i, MU_QMIN) > 1e-6
					fprintf(fd, '%8.3f', gen(i, MU_QMIN));
				else
					fprintf(fd, '     -  ');
				end
				if gen(i, QG)
					fprintf(fd, '%10.2f%10.2f%10.2f', gen(i, [QMIN, QG, QMAX]));
				else
					fprintf(fd, '%10.2f       -  %10.2f', gen(i, [QMIN, QMAX]));
				end
				if gen(i, MU_QMAX) > 1e-6
					fprintf(fd, '%9.3f', gen(i, MU_QMAX));
				else
					fprintf(fd, '      -  ');
				end
			end
		end
		fprintf(fd, '\n');
	end
		
	%% line flow constraints
	if OUT_LINE_LIM == 2 | (OUT_LINE_LIM == 1 & ...
						 (any(branch(:, MU_SF) > 1e-6) | any(branch(:, MU_ST) > 1e-6)))
		fprintf(fd, '\n================================================================================');
		fprintf(fd, '\n|     Branch Flow Constraints                                                  |');
		fprintf(fd, '\n================================================================================');
		fprintf(fd, '\nBrnch   From     "From" End        Limit       "To" End        To');
		fprintf(fd, '\n  #     Bus   |Sf| mu    |Sf|     |Smax|     |St|    |St| mu   Bus');
		fprintf(fd, '\n-----  -----  -------  --------  --------  --------  -------  -----');
		for i = 1:nl
			if OUT_LINE_LIM == 2 | (OUT_LINE_LIM == 1 & ...
						 (branch(i, MU_SF) > 1e-6 | branch(i, MU_ST) > 1e-6))
				fprintf(fd, '\n%4d%7d', i, branch(i, F_BUS));
				if branch(i, MU_SF) > 1e-6
					fprintf(fd, '%10.3f', branch(i, MU_SF));
				else
					fprintf(fd, '      -   ');
				end
				fprintf(fd, '%9.2f%10.2f%10.2f', ...
					[abs(branch(i, PF) + j * branch(i, QF)), ...
					branch(i, RATE_A), abs(branch(i, PT) + j * branch(i, QT))]);
				if branch(i, MU_ST) > 1e-6
					fprintf(fd, '%10.3f', branch(i, MU_ST));
				else
					fprintf(fd, '      -   ');
				end
				fprintf(fd, '%6d', branch(i, T_BUS));
			end
		end
		fprintf(fd, '\n');
	end
end

%% print raw data for Perl database interface
if OUT_RAW
	fprintf(fd, '----------  raw PB::Soln data below  ----------\n');
	fprintf(fd, 'bus\n');
	if isOPF
		fprintf(fd, '%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\n', ...
					bus(:, [BUS_I, BUS_TYPE, VM, VA, LAM_P, LAM_Q, MU_VMAX, MU_VMIN])');
	
		fprintf(fd, 'branch\n');
		fprintf(fd, '%d\t%g\t%g\t%g\t%g\t%g\t%g\n', ...
					[[1:nl]' branch(:, [PF, QF, PT, QT, MU_SF, MU_ST])]');
	
		fprintf(fd, 'gen\n');
		fprintf(fd, '%d\t%g\t%g\t%g\t%d\t%g\t%g\t%g\t%g\n', ...
					[[1:ng]' gen(:, [PG, QG, VG, GEN_STATUS, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN])]');
	else
		fprintf(fd, '%d\t%d\t%f\t%f\t%d\t%d\t%d\t%d\n', ...
					[bus(:, [BUS_I, BUS_TYPE, VM, VA]) zeros(nb, 4)]');
	
		fprintf(fd, 'branch\n');
		fprintf(fd, '%d\t%f\t%f\t%f\t%f\t%d\t%d\n', ...
					[[1:nl]' branch(:, [PF, QF, PT, QT]) zeros(nl, 2)]');
	
		fprintf(fd, 'gen\n');
		fprintf(fd, '%d\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\n', ...
					[[1:ng]' gen(:, [PG, QG, VG, GEN_STATUS]) zeros(ng, 4)]');
	end
	fprintf(fd, 'end\n');
	fprintf(fd, '----------  raw PB::Soln data above  ----------\n');
end

return;
