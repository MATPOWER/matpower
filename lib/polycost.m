function f = polycost(gencost, Pg, der)
%POLYCOST  Evaluates polynomial generator cost & derivatives.
%   F = POLYCOST(GENCOST, PG) returns the vector of costs evaluated at PG
%
%   DF = POLYCOST(GENCOST, PG, 1) returns the vector of first derivatives
%   of costs evaluated at PG
%
%   D2F = POLYCOST(GENCOST, PG, 2) returns the vector of second derivatives
%   of costs evaluated at PG
%
%   GENCOST must contain only polynomial costs
%   PG is in MW, not p.u. (works for QG too)
%
%   This is a more effecient implementation that what can be done with
%   MATLAB's built-in POLYVAL and POLYDER functions.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2009-2010 by Power System Engineering Research Center (PSERC)
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

%%----- initialize -----
%% define named indices into data matrices
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

if nargin < 3
    der = 0;
end

if any(gencost(:, MODEL) == PW_LINEAR)
    error('polycost: all costs must be polynomial');
end

ng = length(Pg);
maxN = max(gencost(:, NCOST));
minN = min(gencost(:, NCOST));

%% form coefficient matrix where 1st column is constant term, 2nd linear, etc.
c = zeros(ng, maxN);
for n = minN:maxN
    k = find(gencost(:, NCOST) == n);   %% cost with n coefficients
    c(k, 1:n) = gencost(k, (COST+n-1):-1:COST);
end

%% do derivatives
for d = 1:der
    if size(c, 2) >= 2
        c = c(:, 2:maxN-d+1);
    else
        c = zeros(ng, 1);
        break;
    end
    for k = 2:maxN-d
        c(:, k) = k * c(:, k);
    end
end

%% evaluate polynomial
if isempty(c)
    f = zeros(size(Pg));
else
    f = c(:, 1);        %% constant term
    for k = 2:size(c, 2)
        f = f + c(:, k) .* Pg .^ (k-1);
    end
end
