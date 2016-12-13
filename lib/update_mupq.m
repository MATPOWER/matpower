function gen = update_mupq(baseMVA, gen, mu_PQh, mu_PQl, data)
%UPDATE_MUPQ  Updates values of generator limit shadow prices.
%   GEN = UPDATE_MUPQ(BASEMVA, GEN, MU_PQH, MU_PQL, DATA)
%
%   Updates the values of MU_PMIN, MU_PMAX, MU_QMIN, MU_QMAX based
%   on any shadow prices on the sloped portions of the generator
%   capability curve constraints.
%
%   MU_PQH - shadow prices on upper sloped portion of capability curves
%   MU_PQL - shadow prices on lower sloped portion of capability curves
%   DATA   - "data" struct returned by MAKEAPQ
%
%   See also MAKEAPQ.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% extract the constraint parameters
[ipqh, ipql, Apqhdata, Apqldata] = ...
    deal(data.ipqh, data.ipql, data.h, data.l);

%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% combine original limit multipliers into single value
muP = gen(:, MU_PMAX) - gen(:, MU_PMIN);
muQ = gen(:, MU_QMAX) - gen(:, MU_QMIN);

%% add P and Q components of multipliers on upper sloped constraint
if ~isempty(ipqh)
    muP(ipqh) = muP(ipqh) - mu_PQh .* Apqhdata(:,1)/baseMVA;
    muQ(ipqh) = muQ(ipqh) - mu_PQh .* Apqhdata(:,2)/baseMVA;
end

%% add P and Q components of multipliers on lower sloped constraint
if ~isempty(ipql)
    muP(ipql) = muP(ipql) - mu_PQl .* Apqldata(:,1)/baseMVA;
    muQ(ipql) = muQ(ipql) - mu_PQl .* Apqldata(:,2)/baseMVA;
end

%% split back into upper and lower multipliers based on sign
gen(:, MU_PMAX) = (muP > 0) .*  muP;
gen(:, MU_PMIN) = (muP < 0) .* -muP;
gen(:, MU_QMAX) = (muQ > 0) .*  muQ;
gen(:, MU_QMIN) = (muQ < 0) .* -muQ;
