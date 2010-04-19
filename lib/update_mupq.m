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
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Autonoma de Manizales
%   Copyright (c) 1996-2010 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as M-files and MEX-files) available in a
%   Matlab (or compatible) environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.


[ipqh, ipql, Apqhdata, Apqldata] = ...
    deal(data.ipqh, data.ipql, data.h, data.l);

%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

% If we succeeded and there were generators with general PQ curve
% characteristics, this is the time to re-compute the multipliers,
% splitting any nonzero multiplier on one of the linear bounds among the
% Pmax, Pmin, Qmax or Qmin limits, producing one multiplier for a P limit and
% another for a Q limit. For upper Q limit, if we are neither at Pmin nor at 
% Pmax, the limit is taken as Pmin if the Qmax line's normal has a negative P
% component, Pmax if it has a positive P component. Messy but there really
% are many cases.
muPmax = gen(:, MU_PMAX);
muPmin = gen(:, MU_PMIN);
if ~isempty(mu_PQh)
%   gen(:, [MU_PMIN MU_PMAX MU_QMIN MU_QMAX])
  k = 1;
  for i = ipqh'
    if muPmax(i) > 0
      gen(i,MU_PMAX)=gen(i,MU_PMAX)-mu_PQh(k)*Apqhdata(k,1)/baseMVA;
    elseif muPmin(i) > 0
      gen(i,MU_PMIN)=gen(i,MU_PMIN)+mu_PQh(k)*Apqhdata(k,1)/baseMVA;
    else
      if Apqhdata(k, 1) >= 0
         gen(i,MU_PMAX)=gen(i,MU_PMAX)-mu_PQh(k)*Apqhdata(k,1)/baseMVA;
      else
         gen(i,MU_PMIN)=gen(i,MU_PMIN)+mu_PQh(k)*Apqhdata(k,1)/baseMVA;
      end
    end
    gen(i,MU_QMAX)=gen(i,MU_QMAX)-mu_PQh(k)*Apqhdata(k,2)/baseMVA;
    k = k + 1;
  end
end

if ~isempty(mu_PQl)
%   gen(:, [MU_PMIN MU_PMAX MU_QMIN MU_QMAX])
  k = 1;
  for i = ipql'
    if muPmax(i) > 0
      gen(i,MU_PMAX)=gen(i,MU_PMAX)-mu_PQl(k)*Apqldata(k,1)/baseMVA;
    elseif muPmin(i) > 0
      gen(i,MU_PMIN)=gen(i,MU_PMIN)+mu_PQl(k)*Apqldata(k,1)/baseMVA;
    else
      if Apqldata(k,1) >= 0
        gen(i,MU_PMAX)=gen(i,MU_PMAX)-mu_PQl(k)*Apqldata(k,1)/baseMVA;
      else
        gen(i,MU_PMIN)=gen(i,MU_PMIN)+mu_PQl(k)*Apqldata(k,1)/baseMVA;
      end
    end
    gen(i,MU_QMIN)=gen(i,MU_QMIN)+mu_PQl(k)*Apqldata(k,2)/baseMVA;
    k = k + 1;
  end
%   gen(:, [MU_PMIN MU_PMAX MU_QMIN MU_QMAX])
%   -[ mu_PQl(1:2) mu_PQh(1:2) ]/baseMVA
%   -[ mu_PQl(1:2).*Apqldata(1:2,1) mu_PQh(1:2).*Apqhdata(1:2,1) ]/baseMVA
%   -[ mu_PQl(1:2).*Apqldata(1:2,2) mu_PQh(1:2).*Apqhdata(1:2,2) ]/baseMVA
end
