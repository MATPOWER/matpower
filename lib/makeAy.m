function [Ay, by]  = makeAy(baseMVA, ng, gencost, pgbas, qgbas, ybas)
% makeAy:  Make the A matrix and RHS for the CCV formulation.
%
% [Ay, by]  = makeAy(ng, gencost, pgbas, qgbas, ybas) constructs
% a matrix Ay and vector by such that the "basin constraints" on Pg,
% Qg and Y that the CCV cost formulation uses are expressed as
%  Ay * x <= by
% where x are the optimization variables.  The starting index within the x
% vector for the active, reactive sources and the y variables should be 
% provided in arguments pgbas, qgbas, ybas. The number of generators is ng.
%
% Assumptions: all generators are committed.  Filter any generators
% that are offline from the gencost() matrix before calling makeAy.
% Efficiency depends on Qg variables being after Pg variables,
% and the Y variables must be the last variables within the vector x for
% the dimensions of the resulting Ay to be conformable with x .

%   MATPOWER
%   $Id$
%   by Carlos Murillo-Sanchez, PSERC Cornell
%   Copyright (c) 1996-2003 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

iycost = find(gencost(:, MODEL) == PW_LINEAR);
ny = size(iycost, 1);

if ny == 0
   Ay = sparse([], [], [], 0, 0, 0);;
   by = [];
   return
end

% if p(i),p(i+1),c(i),c(i+1) define one of the cost segments, then
% the corresponding constraint on Pg and Y is
%                                             c(i+1) - c(i)
%  Y   >=   c(i) + m * (Pg - p(i)),      m = ---------------
%                                             p(i+1) - p(i)
%
% this becomes   m * Pg - Y   <=   m*p(i) - c(i)

% Form A matrix.  Use two different loops, one for the PG/Qg coefs,
% then another for the y coefs so that everything is filled in the
% same order as the compressed column sparse format used by matlab;
% this should be the quickest.

m = sum(gencost(iycost, NCOST));
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
     sidx = qgbas + (i-ng) - 1;
   else
     sidx = pgbas + i - 1;
   end
%   if (length(k:k+ns-2) ~= length(m)) | isempty(m)
%     keyboard
%   end
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

