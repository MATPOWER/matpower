function totalcost = totcost(gencost, Pg)
%TOTCOST    Computes total cost for generators at given output level.
%   totalcost = totcost(gencost, Pg) computes total cost for generators given
%   a matrix in gencost format and a column vector or matrix of generation
%   levels. The return value has the same dimensions as Pg. Each row
%   of gencost is used to evaluate the cost at the points specified in the
%   corresponding row of Pg.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   & Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Autonoma de Manizales
%   Copyright (c) 1996-2006 by Power System Engineering Research Center (PSERC)
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
        j1 = 1:(gencost(i, NCOST) - 1);    j2 = 2:gencost(i, NCOST);
        pp = mkpp(x(i, 1:gencost(i, NCOST))', [(y(i,j2) - y(i,j1)) ./ (x(i,j2) - x(i,j1));  y(i,j1)]');
        totalcost(i,:) = ppval(pp, Pg(i,:));
      end
    end
  end
  tc = totalcost;
  for i = 1:size(totalcost, 2)
  	totalcost(ipol, i) = polycost(gencost(ipol, :), Pg(ipol, i));
  end
end
