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
	1	3	0.0	0.0	0.0	0.0	1	1.0000	0.0000	345.0000	1	1.1000	0.9000;
	2	2	0.0	0.0	0.0	0.0	1	1.0000	0.0000	345.0000	1	1.1000	0.9000;
	3	2	0.0	0.0	0.0	0.0	1	1.0000	0.0000	345.0000	1	1.1000	0.9000;
	4	1	0.0	0.0	0.0	0.0	1	1.0000	0.0000	345.0000	1	1.1000	0.9000;
	5	1	90.0000	30.0000	0.0	0.0	1	1.0000	0.0000	345.0000	1	1.1000	0.9000;
	6	1	0.0	0.0	0.0	0.0	1	1.0000	0.0000	345.0000	1	1.1000	0.9000;
	7	1	100.0000	35.0000	0.0	0.0	1	1.0000	0.0000	345.0000	1	1.1000	0.9000;
	8	1	0.0	0.0	0.0	0.0	1	1.0000	0.0000	345.0000	1	1.1000	0.9000;
	9	1	125.0000	50.0000	0.0	0.0	1	1.0000	0.0000	345.0000	1	1.1000	0.9000;
];

%% generator data
gen = [
	1	0.0000	0.0000	300.0000	-300.0000	1.0000	100.0000	1	250.0000	10.0000;
	2	163.0000	0.0000	300.0000	-300.0000	1.0000	100.0000	1	300.0000	10.0000;
	3	85.0000	0.0000	300.0000	-300.0000	1.0000	100.0000	1	270.0000	10.0000;
];

%% branch data
branch = [
	1	4	0.0000	0.0576	0.0000	250.0000	250.0000	250.0000	0.0000	0.0000	1;
	4	5	0.0170	0.0920	0.1580	250.0000	250.0000	250.0000	0.0000	0.0000	1;
	5	6	0.0390	0.1700	0.3580	150.0000	150.0000	150.0000	0.0000	0.0000	1;
	3	6	0.0000	0.0586	0.0000	300.0000	300.0000	300.0000	0.0000	0.0000	1;
	6	7	0.0119	0.1008	0.2090	150.0000	150.0000	150.0000	0.0000	0.0000	1;
	7	8	0.0085	0.0720	0.1490	250.0000	250.0000	250.0000	0.0000	0.0000	1;
	8	2	0.0000	0.0625	0.0000	250.0000	250.0000	250.0000	0.0000	0.0000	1;
	8	9	0.0320	0.1610	0.3060	250.0000	250.0000	250.0000	0.0000	0.0000	1;
	9	4	0.0100	0.0850	0.1760	250.0000	250.0000	250.0000	0.0000	0.0000	1;
];

%%-----  OPF Data  -----%%
%% area data
area = [
	1	5;
];

%% generator cost data
gencost = [
	2	1500.00	0.00	3	0.11	5	150;
	2	2000.00	0.00	3	0.085	1.2	600;
	2	3000.00	0.00	3	0.1225	1	335;
];

return;
