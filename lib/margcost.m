function marginalcost = margcost(gencost, Pg)
%MARGCOST    Computes marginal cost for generators at given output level.
%   MARGINALCOST = MARGCOST(GENCOST, PG) computes marginal cost for generators
%   given a matrix in gencost format and a column vector of generation levels.
%   The return value has the same dimensions as PG. Each row of GENCOST is
%   used to evaluate the cost at the points specified in the corresponding row
%   of PG.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   & Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Autonoma de Manizales
%   Copyright (c) 1996-2012 by Power System Engineering Research Center (PSERC)
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
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

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
