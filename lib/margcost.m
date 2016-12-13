function marginalcost = margcost(gencost, Pg)
%MARGCOST    Computes marginal cost for generators at given output level.
%   MARGINALCOST = MARGCOST(GENCOST, PG) computes marginal cost for generators
%   given a matrix in gencost format and a column vector of generation levels.
%   The return value has the same dimensions as PG. Each row of GENCOST is
%   used to evaluate the cost at the points specified in the corresponding row
%   of PG.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   & Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

[ng, m] = size(gencost);
marginalcost = zeros(ng, 1);

if ~isempty(gencost)
  ipwl = find(gencost(:, MODEL) == PW_LINEAR);
  ipol = find(gencost(:, MODEL) == POLYNOMIAL);
  if ~isempty(ipwl)
    x = gencost(:, COST:2:(m-1));
    y = gencost(:, (COST+1):2:m);
    for i = ipwl'
      if gencost(i, NCOST) > 0
        c = diff(y(i,:)) ./ diff(x(i,:));
        k = find(Pg(i,:) <= x(i,:), 1);
        if isempty(k)
          marginalcost(i,:) = c(end);
        elseif k == 1
          marginalcost(i,:) = c(1);
        else
          marginalcost(i,:) = c(k-1);
        end
      end
    end
  end
  for i = ipol'
    marginalcost(i,:) = polyval(polyder(gencost(i, COST:(COST+gencost(i, NCOST)-1) )), Pg(i,:) );
  end
end
