function totalcost = totcost(gencost, Pg)
%TOTCOST    Computes total cost for generators at given output level.
%   TOTALCOST = TOTCOST(GENCOST, PG) computes total cost for generators given
%   a matrix in gencost format and a column vector or matrix of generation
%   levels. The return value has the same dimensions as PG. Each row
%   of GENCOST is used to evaluate the cost at the points specified in the
%   corresponding row of PG.

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
totalcost = zeros(ng, size(Pg, 2));

if ~isempty(gencost)
  ipwl = find(gencost(:, MODEL) == PW_LINEAR);
  ipol = find(gencost(:, MODEL) == POLYNOMIAL);
  if ~isempty(ipwl)
    x = gencost(:, COST:2:(m-1));
    y = gencost(:, (COST+1):2:m);
    for i = ipwl'
      if gencost(i, NCOST) > 0
%         totalcost(i,:) = interp1(x(i, 1:gencost(i, NCOST)), ...
%                                  y(i, 1:gencost(i, NCOST)), ...
%                                  Pg(i,:), 'linear', 'extrap');
        j1 = 1:(gencost(i, NCOST) - 1);    j2 = 2:gencost(i, NCOST);
        pp = mkpp(x(i, 1:gencost(i, NCOST))', [(y(i,j2) - y(i,j1)) ./ (x(i,j2) - x(i,j1));  y(i,j1)]');
        totalcost(i,:) = ppval(pp, Pg(i,:));
      end
    end
  end
  for i = 1:size(totalcost, 2)
    totalcost(ipol, i) = polycost(gencost(ipol, :), Pg(ipol, i));
  end
end
