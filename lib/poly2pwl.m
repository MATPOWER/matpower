function pwlcost = poly2pwl(polycost, Pmin, Pmax, npts)
%POLY2PWL  Converts polynomial cost variable to piecewise linear.
%   PWLCOST = POLY2PWL(POLYCOST, PMIN, PMAX, NPTS) converts the polynomial
%   cost variable POLYCOST into a piece-wise linear cost by evaluating at
%   NPTS evenly spaced points between PMIN and PMAX. If the range does not
%   include 0, then it is evaluated at 0 and NPTS-1 evenly spaced points
%   between PMIN and PMAX.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

pwlcost = polycost;
[m, n] = size(polycost);            %% size of piece being changed
pwlcost(:, MODEL)  = PW_LINEAR;     %% change cost model
pwlcost(:, COST:n) = 0;             %% zero out old data
pwlcost(:, NCOST)  = npts;          %% change number of data points
pwlcost(1, COST+2*(npts-1)+1) = 0;  %% expand as needed

for i = 1:m
    if Pmin(i) > 0
        step = (Pmax(i) - Pmin(i)) / (npts - 2);
        xx = [0 Pmin(i):step:Pmax(i)];
    elseif Pmax(i) < 0
        step = (Pmax(i) - Pmin(i)) / (npts - 2);
        xx = [Pmin(i):step:Pmax(i) 0];
    else
        step = (Pmax(i) - Pmin(i)) / (npts - 1);
        xx = (Pmin(i):step:Pmax(i));
    end
    yy = totcost(polycost(i, :), xx);
    pwlcost(i,      COST:2:(COST + 2*(npts-1)    )) = xx;
    pwlcost(i,  (COST+1):2:(COST + 2*(npts-1) + 1)) = yy;
end
