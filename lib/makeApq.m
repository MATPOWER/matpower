function [Apqh, ubpqh, Apql, ubpql, data] = makeApq(baseMVA, gen)
%MAKEAPQ Construct linear constraints for generator capability curves.
%   [APQH, UBPQH, APQL, UBPQL, DATA] = MAKEAPQ(BASEMVA, GEN)
%
%   Constructs the parameters for the following linear constraints
%   implementing trapezoidal generator capability curves, where
%   Pg and Qg are the real and reactive generator injections.
%
%   APQH * [Pg; Qg] <= UBPQH
%   APQL * [Pg; Qg] <= UBPQL
%
%   DATA constains additional information as shown below.
%
%   Example:
%       [Apqh, ubpqh, Apql, ubpql, data] = makeApq(baseMVA, gen);
%
%       data.h      [QC1MAX-QC2MAX, PC2-PC1]
%       data.l      [QC2MIN-QC1MIN, PC1-PC2]
%       data.ipqh   indices of gens with general PQ cap curves (upper)
%       data.ipql   indices of gens with general PQ cap curves (lower)

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% data dimensions
ng = size(gen, 1);      %% number of dispatchable injections

%% which generators require additional linear constraints
%% (in addition to simple box constraints) on (Pg,Qg) to correctly
%% model their PQ capability curves
ipqh = find( hasPQcap(gen, 'U') );
ipql = find( hasPQcap(gen, 'L') );
npqh = size(ipqh, 1);   %% number of general PQ capability curves (upper)
npql = size(ipql, 1);   %% number of general PQ capability curves (lower)

%% make Apqh if there is a need to add general PQ capability curves;
%% use normalized coefficient rows so multipliers have right scaling
%% in $$/pu
if npqh > 0
  data.h = [gen(ipqh,QC1MAX)-gen(ipqh,QC2MAX), gen(ipqh,PC2)-gen(ipqh,PC1)];
  ubpqh = data.h(:, 1) .* gen(ipqh,PC1) + data.h(:, 2) .* gen(ipqh,QC1MAX);
  for i=1:npqh,
    tmp = norm(data.h(i,:));
    data.h(i,:) = data.h(i, :) / tmp;
    ubpqh(i) = ubpqh(i) / tmp;
  end
  Apqh = sparse([1:npqh, 1:npqh]', [ipqh; ipqh+ng], ...
                data.h(:), npqh, 2*ng);
  ubpqh = ubpqh / baseMVA;
else
  data.h = [];
  Apqh  = sparse(0, 2*ng);
  ubpqh = [];
end

%% similarly Apql
if npql > 0
  data.l = [gen(ipql,QC2MIN)-gen(ipql,QC1MIN), gen(ipql,PC1)-gen(ipql,PC2)];
  ubpql= data.l(:, 1) .* gen(ipql,PC1) + data.l(:, 2) .* gen(ipql,QC1MIN) ;
  for i=1:npql,
    tmp = norm(data.l(i, : ));
    data.l(i, :) = data.l(i, :) / tmp;
    ubpql(i) = ubpql(i) / tmp;
  end
  Apql = sparse([1:npql, 1:npql]', [ipql; ipql+ng], ...
                data.l(:), npql, 2*ng);
  ubpql = ubpql / baseMVA;
else
  data.l = [];
  Apql  = sparse(0, 2*ng);
  ubpql = [];
end

data.ipql = ipql;
data.ipqh = ipqh;
