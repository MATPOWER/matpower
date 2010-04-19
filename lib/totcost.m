function totalcost = totcost(gencost, Pg)
%TOTCOST    Computes total cost for generators at given output level.
%   TOTALCOST = TOTCOST(GENCOST, PG) computes total cost for generators given
%   a matrix in gencost format and a column vector or matrix of generation
%   levels. The return value has the same dimensions as PG. Each row
%   of GENCOST is used to evaluate the cost at the points specified in the
%   corresponding row of PG.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   & Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Autonoma de Manizales
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
  for i = 1:size(totalcost, 2)
    totalcost(ipol, i) = polycost(gencost(ipol, :), Pg(ipol, i));
  end
end
