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
%   See http://www.pserc.cornell.edu/matpower/ for more info.

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
