%CASEFORMAT    Defines the MATPOWER case file format.
%   A MATPOWER case file is an M-file or MAT-file that defines or returns
%   a struct named mpc, referred to as a "MATPOWER case struct". The fields
%   of this struct are baseMVA, bus, gen, branch, and (optional) gencost. With
%   the exception of baseMVA, a scalar, each data variable is a matrix, where
%   a row corresponds to a single bus, branch, gen, etc. The format of the
%   data is similar to the PTI format described in
%       http://www.ee.washington.edu/research/pstca/formats/pti.txt
%   except where noted. An item marked with (+) indicates that it is included
%   in this data but is not part of the PTI format. An item marked with (-) is
%   one that is in the PTI format but is not included here. Those marked with
%   (2) were added for version 2 of the case file format. The columns for
%   each data matrix are given below.
%
%   MATPOWER Case Version Information:
%   There are two versions of the MATPOWER case file format. The current
%   version of MATPOWER uses version 2 of the MATPOWER case format
%   internally, and includes a 'version' field with a value of '2' to make
%   the version explicit. Earlier versions of MATPOWER used the version 1
%   case format, which defined the data matrices as individual variables,
%   as opposed to fields of a struct. Case files in version 1 format with
%   OPF data also included an (unused) 'areas' variable. While the version 1
%   format has now been deprecated, it is still be handled automatically by
%   LOADCASE and SAVECASE which are able to load and save case files in both
%   version 1 and version 2 formats.
%
%   See also IDX_BUS, IDX_BRCH, IDX_GEN, IDX_AREA and IDX_COST regarding
%   constants which can be used as named column indices for the data matrices.
%   Also described in the first three are additional results columns that
%   are added to the bus, branch and gen matrices by the power flow and OPF
%   solvers.
%
%   The case struct also allows for additional fields to be included.
%   The OPF is designed to recognize fields named A, l, u, H, Cw, N,
%   fparm, z0, zl and zu as parameters used to directly extend the OPF
%   formulation (see OPF for details). Additional standard optional fields
%   include bus_name, gentype and genfuel. Other user-defined fields may also
%   be included and will be automatically loaded by the LOADCASE function
%   and, given an appropriate 'savecase' callback function (see
%   ADD_USERFCN), saved by the SAVECASE function.
%
%   Bus Data Format
%       1   bus number (positive integer)
%       2   bus type
%               PQ bus          = 1
%               PV bus          = 2
%               reference bus   = 3
%               isolated bus    = 4
%       3   Pd, real power demand (MW)
%       4   Qd, reactive power demand (MVAr)
%       5   Gs, shunt conductance (MW demanded at V = 1.0 p.u.)
%       6   Bs, shunt susceptance (MVAr injected at V = 1.0 p.u.)
%       7   area number, (positive integer)
%       8   Vm, voltage magnitude (p.u.)
%       9   Va, voltage angle (degrees)
%   (-)     (bus name)
%       10  baseKV, base voltage (kV)
%       11  zone, loss zone (positive integer)
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
%   (2) 11  Pc1, lower real power output of PQ capability curve (MW)
%   (2) 12  Pc2, upper real power output of PQ capability curve (MW)
%   (2) 13  Qc1min, minimum reactive power output at Pc1 (MVAr)
%   (2) 14  Qc1max, maximum reactive power output at Pc1 (MVAr)
%   (2) 15  Qc2min, minimum reactive power output at Pc2 (MVAr)
%   (2) 16  Qc2max, maximum reactive power output at Pc2 (MVAr)
%   (2) 17  ramp rate for load following/AGC (MW/min)
%   (2) 18  ramp rate for 10 minute reserves (MW)
%   (2) 19  ramp rate for 30 minute reserves (MW)
%   (2) 20  ramp rate for reactive power (2 sec timescale) (MVAr/min)
%   (2) 21  APF, area participation factor
%
%   Branch Data Format
%       1   f, from bus number
%       2   t, to bus number
%   (-)     (circuit identifier)
%       3   r, resistance (p.u.)
%       4   x, reactance (p.u.)
%       5   b, total line charging susceptance (p.u.)
%       6   rateA, MVA rating A (long term rating), set to 0 for unlimited
%       7   rateB, MVA rating B (short term rating), set to 0 for unlimited
%       8   rateC, MVA rating C (emergency rating), set to 0 for unlimited
%       9   tap, transformer off nominal turns ratio, if non-zero
%           (taps at "from" bus, impedance at "to" bus, i.e. if r = x = b = 0,
%            then tap = Vf / Vt; tap = 0 used to indicate transmission
%           line rather than transformer, i.e. mathematically equivalent to
%           transformer with tap = 1)
%       10  shift, transformer phase shift angle (degrees), positive => delay
%   (-)     (Gf, shunt conductance at from bus p.u.)
%   (-)     (Bf, shunt susceptance at from bus p.u.)
%   (-)     (Gt, shunt conductance at to bus p.u.)
%   (-)     (Bt, shunt susceptance at to bus p.u.)
%       11  initial branch status, 1 - in service, 0 - out of service
%   (2) 12  minimum angle difference, angle(Vf) - angle(Vt) (degrees)
%   (2) 13  maximum angle difference, angle(Vf) - angle(Vt) (degrees)
%           (The voltage angle difference is taken to be unbounded below
%            if ANGMIN < -360 and unbounded above if ANGMAX > 360.
%            If both parameters are zero, it is unconstrained.)
%
% (+) Generator Cost Data Format
%       NOTE: If gen has ng rows, then the first ng rows of gencost contain
%       the cost for active power produced by the corresponding generators.
%       If gencost has 2*ng rows then rows ng+1 to 2*ng contain the reactive
%       power costs in the same format.
%       1   model, 1 - piecewise linear, 2 - polynomial
%       2   startup, startup cost in US dollars
%       3   shutdown, shutdown cost in US dollars
%       4   N, number of cost coefficients to follow for polynomial
%           cost function, or number of data points for piecewise linear
%       5 and following, parameters defining total cost function f(p),
%           units of f and p are $/hr and MW (or MVAr), respectively.
%           (MODEL = 1) : p0, f0, p1, f1, ..., pn, fn
%               where p0 < p1 < ... < pn and the cost f(p) is defined by
%               the coordinates (p0,f0), (p1,f1), ..., (pn,fn) of the
%               end/break-points of the piecewise linear cost function
%           (MODEL = 2) : cn, ..., c1, c0
%               n+1 coefficients of an n-th order polynomial cost function,
%               starting with highest order, where cost is
%               f(p) = cn*p^n + ... + c1*p + c0
% 
% (+) Area Data Format (deprecated)
%     (this data is not used by MATPOWER and is no longer necessary for
%      version 2 case files with OPF data).
%       1   i, area number
%       2   price_ref_bus, reference bus for that area
%
%   See also LOADCASE, SAVECASE, IDX_BUS, IDX_BRCH, IDX_GEN, IDX_AREA
%   and IDX_COST.

%   MATPOWER
%   Copyright (c) 1996-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.
