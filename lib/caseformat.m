%CASEFORMAT    Defines the MATPOWER case file format.
%   A MATPOWER case file is an M-file or MAT-file which defines the variables
%   baseMVA, bus, gen, branch, areas, and gencost. With the exception of
%   baseMVA, a scalar, each data variable is a matrix, where a row corresponds
%   to a single bus, branch, gen, etc. The format of the data is similar to
%   the PTI format described in
%       http://www.ee.washington.edu/research/pstca/formats/pti.txt
%   except where noted. An item marked with (+) indicates that it is included
%   in this data but is not part of the PTI format. An item marked with (-) is
%   one that is in the PTI format but is not included here. Those marked with
%   (2) were added for version 2 of the case file format. The columns for
%   each data matrix are given below.
%
%   MATPOWER Case Version Information:
%   A version 1 case file defined the data matrices directly. The last two,
%   areas and gencost, were optional since they were not needed for running
%   a simple power flow. In version 2, each of the data matrices are stored
%   as fields in a struct. It is this struct, rather than the individual
%   matrices that is returned by a version 2 M-casefile. Likewise a version 2
%   MAT-casefile stores a struct named 'mpc' (for MATPOWER case). The struct
%   also contains a 'version' field so MATPOWER knows how to interpret the
%   data. Any case file which does not return a struct, or any struct which
%   does not have a 'version' field is considered to be in version 1 format.
%
%   See also IDX_BUS, IDX_BRCH, IDX_GEN, IDX_AREA and IDX_COST regarding
%   constants which can be used as named column indices for the data matrices.
%   Also described in the first three are additional columns that are added
%   to the bus, branch and gen matrices by the power flow and OPF solvers.
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
%       5   Gs, shunt conductance (MW (demanded) at V = 1.0 p.u.)
%       6   Bs, shunt susceptance (MVAr (injected) at V = 1.0 p.u.)
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
%   (2) 11  Qmax2, maximum reactive power output at Pmax (MVAr)
%   (2) 12  Qmin2, minimum reactive power output at Pmax (MVAr)
%   (2) 13  ramp rate for load following/AGC (MW/min)
%   (2) 14  ramp rate for 10 minute reserves (MW)
%   (2) 15  ramp rate for 30 minute reserves (MW)
%   (2) 16  ramp rate for reactive power (2 sec timescale) (MVAr/min)
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
%   (2) 12  minimum angle difference, angle(Vf) - angle(Vt) (degrees)
%   (2) 13  maximum angle difference, angle(Vf) - angle(Vt) (degrees)
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
%           cost function, or number of data points for piecewise linear
%       5 and following, cost data defining total cost function
%           For polynomial cost:
%                   c2, c1, c0
%           where the polynomial is c0 + c1*P + c2*P^2
%           For piecewise linear cost:
%                   x0, y0, x1, y1, x2, y2, ...
%           where x0 < x1 < x2 < ... and the points (x0,y0), (x1,y1),
%           (x2,y2), ... are the end- and break-points of the cost function.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2005 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.
