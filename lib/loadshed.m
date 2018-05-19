function shed = loadshed(gen, ild)
%LOADSHED  Returns a vector of curtailments of dispatchable loads.
%   SHED = LOADSHED(GEN)
%   SHED = LOADSHED(GEN, ILD)
%
%   Returns a column vector of MW curtailments of dispatchable loads.
%
%   Inputs:
%       GEN - MATPOWER generator matrix
%       ILD - (optional) NLD x 1 vector of generator indices corresponding
%             to the dispatchable loads of interest, default is all
%             dispatchable loads as determined by the ISLOAD() function.
%
%   Output:
%       SHED - NLD x 1 vector of the MW curtailment for each dispatchable
%              load of interest
%
%   Example:
%       total_load_shed = max(loadshed(mpc.gen));

%   MATPOWER
%   Copyright (c) 2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% default input
if nargin < 2
    ild = find(isload(gen));
end

tol = 1e-5;
shed = gen(ild, PG) - gen(ild, PMIN);

%% zero out anything less than tolerance
shed(shed < tol) = 0;
