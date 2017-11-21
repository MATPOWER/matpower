function [Ay, by]  = makeAy(baseMVA, ng, gencost, pgbas, qgbas, ybas)
%MAKEAY  Make the A matrix and RHS for the CCV formulation.
%   [AY, BY]  = MAKEAY(BASEMVA, NG, GENCOST, PGBAS, QGBAS, YBAS)
%
%   Constructs the parameters for linear "basin constraints" on Pg, Qg
%   and Y used by the CCV cost formulation, expressed as
%
%       AY * X <= BY
%
%   where X is the vector of optimization variables. The starting index
%   within the X vector for the active, reactive sources and the Y
%   variables should be provided in arguments PGBAS, QGBAS, YBAS. The
%   number of generators is NG.
%
%   Assumptions: All generators are in-service.  Filter any generators
%   that are offline from the GENCOST matrix before calling MAKEAY.
%   Efficiency depends on Qg variables being after Pg variables, and
%   the Y variables must be the last variables within the vector X for
%   the dimensions of the resulting AY to be conformable with X.
%
%   Example:
%       [Ay, by]  = makeAy(baseMVA, ng, gencost, pgbas, qgbas, ybas);

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

% find all pwl cost rows in gencost, either real or reactive
iycost = find(gencost(:, MODEL) == PW_LINEAR);

% this is the number of extra "y" variables needed to model those costs
ny = size(iycost, 1);

if ny == 0
   Ay = sparse([], [], [], 0, ybas+ny-1, 0);
   by = [];
   return
end

% if p(i),p(i+1),c(i),c(i+1) define one of the cost segments, then
% the corresponding constraint on Pg (or Qg) and Y is
%                                             c(i+1) - c(i)
%  Y   >=   c(i) + m * (Pg - p(i)),      m = ---------------
%                                             p(i+1) - p(i)
%
% this becomes   m * Pg - Y   <=   m*p(i) - c(i)

% Form A matrix.  Use two different loops, one for the Pg/Qg coeffs,
% then another for the y coeffs so that everything is filled in the
% same order as the compressed column sparse format used by MATLAB;
% this should be the quickest.

m = sum(gencost(iycost, NCOST));  % total number of cost points
Ay = sparse([], [], [], m-ny, ybas+ny-1, 2*(m-ny)); 
by = [];
% First fill the Pg or Qg coefficients (since their columns come first)
% and the rhs
k = 1;
for i=iycost'
   ns = gencost(i, NCOST);                % # of cost points; segments = ns-1
   p = gencost(i, COST:2:COST+2*ns-1) / baseMVA;
   c = gencost(i, COST+1:2:COST+2*ns);
   m = diff(c) ./ diff(p);                % slopes for Pg (or Qg)
   if any(diff(p) == 0)
     fprintf('\nmakeAy: bad x axis data in row %i of gencost matrix\n',i);
   end
   b = m .* p(1:ns-1) - c(1:ns-1);        % and rhs
   by = [by;  b'];
   if i > ng
     sidx = qgbas + (i-ng) - 1;           % this was for a q cost
   else
     sidx = pgbas + i - 1;                % this was for a p cost
   end
   Ay(k:k+ns-2, sidx) = m';
   k = k + ns - 1;
end
% Now fill the y columns with -1's
k = 1;
j = 1;
for i=iycost'
   ns = gencost(i, NCOST);
   Ay(k:k+ns-2, ybas+j-1) = -ones(ns-1,1);
   k = k + ns - 1;
   j = j + 1;
end
