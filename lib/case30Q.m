function [baseMVA, bus, gen, branch, area, gencost] = case
%CASE    Defines the power flow data in a format similar to PTI.
%   [baseMVA, bus, gen, branch, area, gencost] = case
%   The format for the data is similar to PTI format except where noted.
%   An item marked with (+) indicates that it is included in this data
%   but is not part of the PTI format. An item marked with (-) is one that
%   is in the PTI format but is not included here.
%
%   Bus Data Format
%       1   bus number (1 to 29997)
%       2   bus type
%               PQ bus          = 1
%               PV bus          = 2
%               reference bus   = 3
%               isolated bus    = 4
%       3   Pd, real power demand (MW)
%       4   Qd, reactive power demand (MVAr)
%       5   Gs, shunt conductance (MW (demanded?) at V = 1.0 p.u.)
%       6   Bs, shunt susceptance (MVAr (injected?) at V = 1.0 p.u.)
%       7   area number, 1-100
%       8   Vm, voltage magnitude (p.u.)
%       9   Va, voltage angle (degrees)
%   (-)     (bus name)
%       10  baseKV, base voltage (kV)
%       11  zone, loss zone (1-999)
%   (+) 12  maxVm, maximum voltage magnitude (p.u.)
%   (+) 13  minVm, minimum voltage magnitude (p.u.)
%
%   Generator Data Format
%       1   bus number
%   (-)     (machine identifier, 0-9, A-Z)
%       2   Pg, real power output (MW)
%       3   Qg, reactive power output (MVAr)
%       4   Qmax, maximum reactive power output (MVAr)
%       5   Qmin, minimum reactive power output (MVAr)
%       6   Vg, voltage magnitude setpoint (p.u.)
%   (-)     (remote controlled bus index)
%       7   mBase, total MVA base of this machine, defaults to baseMVA
%   (-)     (machine impedance, p.u. on mBase)
%   (-)     (step up transformer impedance, p.u. on mBase)
%   (-)     (step up transformer off nominal turns ratio)
%       8   status,  >  0 - machine in service
%                    <= 0 - machine out of service
%   (-)     (% of total VAr's to come from this gen in order to hold V at
%               remote bus controlled by several generators)
%       9   Pmax, maximum real power output (MW)
%       10  Pmin, minimum real power output (MW)
%
%   Branch Data Format
%       1   f, from bus number
%       2   t, to bus number
%   (-)     (circuit identifier)
%       3   r, resistance (p.u.)
%       4   x, reactance (p.u.)
%       5   b, total line charging susceptance (p.u.)
%       6   rateA, MVA rating A (long term rating)
%       7   rateB, MVA rating B (short term rating)
%       8   rateC, MVA rating C (emergency rating)
%       9   ratio, transformer off nominal turns ratio ( = 0 for lines )
%           (taps at 'from' bus, impedance at 'to' bus, i.e. ratio = Vf / Vt)
%       10  angle, transformer phase shift angle (degrees)
%   (-)     (Gf, shunt conductance at from bus p.u.)
%   (-)     (Bf, shunt susceptance at from bus p.u.)
%   (-)     (Gt, shunt conductance at to bus p.u.)
%   (-)     (Bt, shunt susceptance at to bus p.u.)
%       11  initial branch status, 1 - in service, 0 - out of service
%
% (+) Area Data Format
%       1   i, area number
%       2   price_ref_bus, reference bus for that area
% 
% (+) Generator Cost Data Format
%       NOTE: If gen has n rows, then the first n rows of gencost contain
%       the cost for active power produced by the corresponding generators.
%       If gencost has 2*n rows then rows n+1 to 2*n contain the reactive
%       power costs in the same format.
%       1   model, 1 - piecewise linear, 2 - polynomial
%       2   startup, startup cost in US dollars
%       3   shutdown, shutdown cost in US dollars
%       4   n, number of cost coefficients to follow for polynomial
%           (or data points for piecewise linear) total cost function
%       5 and following, cost data, piecewise linear data as:
%                   x0, y0, x1, y1, x2, y2, ...
%           and polynomial data as, e.g.:
%                   c2, c1, c0
%           where the polynomial is c0 + c1*P + c2*P^2
%
% << this file created [97-Aug-26 12:29:15] by PB::System version 1.3 >>

%%-----  Power Flow Data  -----%%
%% system MVA base
baseMVA = 100.0000;

%% bus data
bus = [
	1	3	0.0	0.0	0.0	0.0	1	1.0000	0.0000	135.0000	1	1.0500	0.9500;
	2	2	21.7000	12.7000	0.0	0.0	1	1.0000	0.0000	135.0000	1	1.1000	0.9500;
	3	1	2.4000	1.2000	0.0	0.0	1	1.0000	0.0000	135.0000	1	1.0500	0.9500;
	4	1	7.6000	1.6000	0.0	0.0	1	1.0000	0.0000	135.0000	1	1.0500	0.9500;
	5	1	0.0	0.0	0.0000	0.1900	1	1.0000	0.0000	135.0000	1	1.0500	0.9500;
	6	1	0.0	0.0	0.0	0.0	1	1.0000	0.0000	135.0000	1	1.0500	0.9500;
	7	1	22.8000	10.9000	0.0	0.0	1	1.0000	0.0000	135.0000	1	1.0500	0.9500;
	8	1	30.0000	30.0000	0.0	0.0	1	1.0000	0.0000	135.0000	1	1.0500	0.9500;
	9	1	0.0	0.0	0.0	0.0	1	1.0000	0.0000	135.0000	1	1.0500	0.9500;
	10	1	5.8000	2.0000	0.0	0.0	3	1.0000	0.0000	135.0000	1	1.0500	0.9500;
	11	1	0.0	0.0	0.0	0.0	1	1.0000	0.0000	135.0000	1	1.0500	0.9500;
	12	1	11.2000	7.5000	0.0	0.0	2	1.0000	0.0000	135.0000	1	1.0500	0.9500;
	13	2	0.0	0.0	0.0	0.0	2	1.0000	0.0000	135.0000	1	1.1000	0.9500;
	14	1	6.2000	1.6000	0.0	0.0	2	1.0000	0.0000	135.0000	1	1.0500	0.9500;
	15	1	8.2000	2.5000	0.0	0.0	2	1.0000	0.0000	135.0000	1	1.0500	0.9500;
	16	1	3.5000	1.8000	0.0	0.0	2	1.0000	0.0000	135.0000	1	1.0500	0.9500;
	17	1	9.0000	5.8000	0.0	0.0	2	1.0000	0.0000	135.0000	1	1.0500	0.9500;
	18	1	3.2000	0.9000	0.0	0.0	2	1.0000	0.0000	135.0000	1	1.0500	0.9500;
	19	1	9.5000	3.4000	0.0	0.0	2	1.0000	0.0000	135.0000	1	1.0500	0.9500;
	20	1	2.2000	0.7000	0.0	0.0	2	1.0000	0.0000	135.0000	1	1.0500	0.9500;
	21	1	17.5000	11.2000	0.0	0.0	3	1.0000	0.0000	135.0000	1	1.0500	0.9500;
	22	2	0.0	0.0	0.0	0.0	3	1.0000	0.0000	135.0000	1	1.1000	0.9500;
	23	2	3.2000	1.6000	0.0	0.0	2	1.0000	0.0000	135.0000	1	1.1000	0.9500;
	24	1	8.7000	6.7000	0.0000	0.0400	3	1.0000	0.0000	135.0000	1	1.0500	0.9500;
	25	1	0.0	0.0	0.0	0.0	3	1.0000	0.0000	135.0000	1	1.0500	0.9500;
	26	1	3.5000	2.3000	0.0	0.0	3	1.0000	0.0000	135.0000	1	1.0500	0.9500;
	27	2	0.0	0.0	0.0	0.0	3	1.0000	0.0000	135.0000	1	1.1000	0.9500;
	28	1	0.0	0.0	0.0	0.0	1	1.0000	0.0000	135.0000	1	1.0500	0.9500;
	29	1	2.4000	0.9000	0.0	0.0	3	1.0000	0.0000	135.0000	1	1.0500	0.9500;
	30	1	10.6000	1.9000	0.0	0.0	3	1.0000	0.0000	135.0000	1	1.0500	0.9500;
];

%% generator data
gen = [
	1	23.5400	0.0000	150.0000	-20.0000	1.0000	100.0000	1	80.0000	0.0000;
	2	60.9700	0.0000	60.0000	-20.0000	1.0000	100.0000	1	80.0000	0.0000;
	22	21.5900	0.0000	62.5000	-15.0000	1.0000	100.0000	1	50.0000	0.0000;
	27	26.9100	0.0000	48.7000	-15.0000	1.0000	100.0000	1	55.0000	0.0000;
	23	19.2000	0.0000	40.0000	-10.0000	1.0000	100.0000	1	30.0000	0.0000;
	13	37.0000	0.0000	44.7000	-15.0000	1.0000	100.0000	1	40.0000	0.0000;
];

%% branch data
branch = [
	1	2	0.0200	0.0600	0.0300	130.0000	130.0000	130.0000	0.0000	0.0000	1;
	1	3	0.0500	0.1900	0.0200	130.0000	130.0000	130.0000	0.0000	0.0000	1;
	2	4	0.0600	0.1700	0.0200	65.0000	65.0000	65.0000	0.0000	0.0000	1;
	3	4	0.0100	0.0400	0.0000	130.0000	130.0000	130.0000	0.0000	0.0000	1;
	2	5	0.0500	0.2000	0.0200	130.0000	130.0000	130.0000	0.0000	0.0000	1;
	2	6	0.0600	0.1800	0.0200	65.0000	65.0000	65.0000	0.0000	0.0000	1;
	4	6	0.0100	0.0400	0.0000	90.0000	90.0000	90.0000	0.0000	0.0000	1;
	5	7	0.0500	0.1200	0.0100	70.0000	70.0000	70.0000	0.0000	0.0000	1;
	6	7	0.0300	0.0800	0.0100	130.0000	130.0000	130.0000	0.0000	0.0000	1;
	6	8	0.0100	0.0400	0.0000	32.0000	32.0000	32.0000	0.0000	0.0000	1;
	6	9	0.0000	0.2100	0.0000	65.0000	65.0000	65.0000	0.0000	0.0000	1;
	6	10	0.0000	0.5600	0.0000	32.0000	32.0000	32.0000	0.0000	0.0000	1;
	9	11	0.0000	0.2100	0.0000	65.0000	65.0000	65.0000	0.0000	0.0000	1;
	9	10	0.0000	0.1100	0.0000	65.0000	65.0000	65.0000	0.0000	0.0000	1;
	4	12	0.0000	0.2600	0.0000	65.0000	65.0000	65.0000	0.0000	0.0000	1;
	12	13	0.0000	0.1400	0.0000	65.0000	65.0000	65.0000	0.0000	0.0000	1;
	12	14	0.1200	0.2600	0.0000	32.0000	32.0000	32.0000	0.0000	0.0000	1;
	12	15	0.0700	0.1300	0.0000	32.0000	32.0000	32.0000	0.0000	0.0000	1;
	12	16	0.0900	0.2000	0.0000	32.0000	32.0000	32.0000	0.0000	0.0000	1;
	14	15	0.2200	0.2000	0.0000	16.0000	16.0000	16.0000	0.0000	0.0000	1;
	16	17	0.0800	0.1900	0.0000	16.0000	16.0000	16.0000	0.0000	0.0000	1;
	15	18	0.1100	0.2200	0.0000	16.0000	16.0000	16.0000	0.0000	0.0000	1;
	18	19	0.0600	0.1300	0.0000	16.0000	16.0000	16.0000	0.0000	0.0000	1;
	19	20	0.0300	0.0700	0.0000	32.0000	32.0000	32.0000	0.0000	0.0000	1;
	10	20	0.0900	0.2100	0.0000	32.0000	32.0000	32.0000	0.0000	0.0000	1;
	10	17	0.0300	0.0800	0.0000	32.0000	32.0000	32.0000	0.0000	0.0000	1;
	10	21	0.0300	0.0700	0.0000	32.0000	32.0000	32.0000	0.0000	0.0000	1;
	10	22	0.0700	0.1500	0.0000	32.0000	32.0000	32.0000	0.0000	0.0000	1;
	21	22	0.0100	0.0200	0.0000	32.0000	32.0000	32.0000	0.0000	0.0000	1;
	15	23	0.1000	0.2000	0.0000	16.0000	16.0000	16.0000	0.0000	0.0000	1;
	22	24	0.1200	0.1800	0.0000	16.0000	16.0000	16.0000	0.0000	0.0000	1;
	23	24	0.1300	0.2700	0.0000	16.0000	16.0000	16.0000	0.0000	0.0000	1;
	24	25	0.1900	0.3300	0.0000	16.0000	16.0000	16.0000	0.0000	0.0000	1;
	25	26	0.2500	0.3800	0.0000	16.0000	16.0000	16.0000	0.0000	0.0000	1;
	25	27	0.1100	0.2100	0.0000	16.0000	16.0000	16.0000	0.0000	0.0000	1;
	28	27	0.0000	0.4000	0.0000	65.0000	65.0000	65.0000	0.0000	0.0000	1;
	27	29	0.2200	0.4200	0.0000	16.0000	16.0000	16.0000	0.0000	0.0000	1;
	27	30	0.3200	0.6000	0.0000	16.0000	16.0000	16.0000	0.0000	0.0000	1;
	29	30	0.2400	0.4500	0.0000	16.0000	16.0000	16.0000	0.0000	0.0000	1;
	8	28	0.0600	0.2000	0.0200	32.0000	32.0000	32.0000	0.0000	0.0000	1;
	6	28	0.0200	0.0600	0.0100	32.0000	32.0000	32.0000	0.0000	0.0000	1;
];

%%-----  OPF Data  -----%%
%% area data
area = [
	1	8;
	2	23;
	3	26;
];

%% generator cost data
gencost = [
	2	0.00	0.00	3	0.02	2	0;
	2	0.00	0.00	3	0.0175	1.75	0;
	2	0.00	0.00	3	0.0625	1	0;
	2	0.00	0.00	3	0.0083	3.25	0;
	2	0.00	0.00	3	0.025	3	0;
	2	0.00	0.00	3	0.025	3	0;
	2	0.00	0.00	3	0.02	0	0;
	2	0.00	0.00	3	0.0175	0	0;
	2	0.00	0.00	3	0.0625	0	0;
	2	0.00	0.00	3	0.0083	0	0;
	2	0.00	0.00	3	0.025	0	0;
	2	0.00	0.00	3	0.025	0	0;
];

return;
