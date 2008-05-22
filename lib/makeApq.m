function [Apqh, lbpqh, ubpqh, Apql, lbpql, ubpql, data] = makeApq(baseMVA, gen)
%   MAKEAPQ

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Autonoma de Manizales
%   Copyright (c) 1996-2008 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% data dimensions
ng = size(gen, 1);      %% number of dispatchable injections

% Find out which generators require additional linear constraints
% (as opposed to simple box constraints) on (Pg,Qg) to correctly
% model their PQ capability curves
ipqh = find( hasPQcap(gen, 'U') );
ipql = find( hasPQcap(gen, 'L') );

npqh = size(ipqh, 1);   %% number of general PQ capability curves (upper)
npql = size(ipql, 1);   %% number of general PQ capability curves (lower)


% Make Apqh if there is a need to add general PQ capability curves;
% use normalized coefficient rows so multipliers have right scaling
% in $$/pu
if npqh > 0
  data.h = [gen(ipqh,QC1MAX)-gen(ipqh,QC2MAX), gen(ipqh,PC2)-gen(ipqh,PC1)];
  ubpqh = (gen(ipqh,QC1MAX)-gen(ipqh,QC2MAX)) .* gen(ipqh,PC1) ...
         + (gen(ipqh,PC2)-gen(ipqh,PC1)) .* gen(ipqh,QC1MAX);
  for i=1:npqh,
    tmp = norm(data.h(i,:));
    data.h(i,:) = data.h(i, :) / tmp;
    ubpqh(i) = ubpqh(i) / tmp;
  end
  Apqh = sparse([1:npqh, 1:npqh]', [ipqh; ipqh+ng], ...
                data.h(:), npqh, 2*ng);
  ubpqh = ubpqh / baseMVA;
  lbpqh = -1e10*ones(npqh,1);
else
  data.h = [];
  Apqh = [];
  ubpqh = [];
  lbpqh = [];
end

% similarly Apql
if npql > 0
  data.l = [gen(ipql,QC2MIN)-gen(ipql,QC1MIN), gen(ipql,PC1)-gen(ipql,PC2)];
  ubpql= (gen(ipql,QC2MIN)-gen(ipql,QC1MIN)) .* gen(ipql,PC1) ...
         - (gen(ipql,PC2)-gen(ipql,PC1)) .* gen(ipql,QC1MIN) ;
  for i=1:npql,
    tmp = norm(data.l(i, : ));
    data.l(i, :) = data.l(i, :) / tmp;
    ubpql(i) = ubpql(i) / tmp;
  end
  Apql = sparse([1:npql, 1:npql]', [ipql; ipql+ng], ...
                data.l(:), npql, 2*ng);
  ubpql = ubpql / baseMVA;
  lbpql = -1e10*ones(npql,1);
else
  data.l = [];
  Apql = [];
  ubpql = [];
  lbpql = [];
end

data.ipql = ipql;
data.ipqh = ipqh;

return;
