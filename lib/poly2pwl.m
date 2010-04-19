function pwlcost = poly2pwl(polycost, Pmin, Pmax, npts)
%POLY2PWL  Converts polynomial cost variable to piecewise linear.
%   PWLCOST = POLY2PWL(POLYCOST, PMIN, PMAX, NPTS) converts the polynomial
%   cost variable POLYCOST into a piece-wise linear cost by evaluating at
%   zero and then at NPTS evenly spaced points between PMIN and PMAX. If
%   PMIN <= 0 (such as for reactive power, where P really means Q) it just
%   uses NPTS evenly spaced points between PMIN and PMAX.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
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

pwlcost = polycost;
[m, n] = size(polycost);                                %% size of piece being changed
pwlcost(:, MODEL)  = PW_LINEAR * ones(m, 1);            %% change cost model
pwlcost(:, COST:n) = zeros(size(pwlcost(:, COST:n)));   %% zero out old data
pwlcost(:, NCOST)  = npts * ones(m, 1);                 %% change number of data points

for i = 1:m
    if Pmin(i) == 0
        step = (Pmax(i) - Pmin(i)) / (npts - 1);
        xx = (Pmin(i):step:Pmax(i));
    elseif Pmin(i) > 0
        step = (Pmax(i) - Pmin(i)) / (npts - 2);
        xx = [0 Pmin(i):step:Pmax(i)];
    elseif Pmin(i) < 0 && Pmax(i) > 0        %% for when P really means Q
        step = (Pmax(i) - Pmin(i)) / (npts - 1);
        xx = (Pmin(i):step:Pmax(i));
    end
    yy = totcost(polycost(i, :), xx);
    pwlcost(i,      COST:2:(COST + 2*(npts-1)    )) = xx;
    pwlcost(i,  (COST+1):2:(COST + 2*(npts-1) + 1)) = yy;
end
