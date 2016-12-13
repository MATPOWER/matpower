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
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

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
