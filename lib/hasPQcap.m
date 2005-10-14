function TorF = hasPQcap(gen, hilo)
%HASPQCAP  Checks for P-Q capability curve constraints.
%   TorF = hasPQcap(gen, hilo) returns a column vector of 1's and 0's. The 1's
%   correspond to rows of the gen matrix which correspond to generators which
%   have defined a capability curve (with sloped upper and/or lower bound on
%   Q).
%
%   The gen matrix in version 2 of the MATPOWER case format includes columns
%   for specifying a P-Q capability curve for a generator defined as the
%   intersection of two half-planes and the box constraints on P and Q. The
%   two half planes are defined respectively as the area below the line
%   connecting (Pc1, Qc1max) and (Pc2, Qc2max) and the area above the line
%   connecting (Pc1, Qc1min) and (Pc2, Qc2min).
%
%   If the optional second argument is 'U' this function returns true only for
%   rows with a sloped upper portion. If it is 'L', it only looks at the lower
%   portion. If it is absent or has any other value it returns true for gens
%   specifying either upper or lower sloped capability curves.
%
%   It is smart enough to return true only if the corresponding linear
%   constraint is not redundant w.r.t the box constraints.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2005 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% default value
if nargin < 2
    hilo = 'B';     %% look at both top and bottom by default
end

%% check for errors capability curve data
if any( gen(:, PC1) > gen(:, PC2) )
    error('hasPQcap: Pc1 > Pc2');
end
if any( gen(:, QC2MAX) > gen(:, QC1MAX) )
    error('hasPQcap: Qc2max > Qc1max');
end
if any( gen(:, QC2MIN) < gen(:, QC1MIN) )
    error('hasPQcap: Qc2min < Qc1min');
end

L = zeros(size(gen, 1), 1);
U = zeros(size(gen, 1), 1);
k = find( gen(:, PC1) ~= gen(:, PC2) );

if ~strcmp(hilo, 'U')       %% include lower constraint
    Qmin_at_Pmax = gen(k, QC1MIN) + (gen(k, PMAX) - gen(k, PC1)) .* ...
        (gen(k, QC2MIN) - gen(k, QC1MIN)) ./ (gen(k, PC2) - gen(k, PC1));
    L(k) = Qmin_at_Pmax > gen(k, QMIN);
end

if ~strcmp(hilo, 'L')       %% include upper constraint
    Qmax_at_Pmax = gen(k, QC1MAX) + (gen(k, PMAX) - gen(k, PC1)) .* ...
        (gen(k, QC2MAX) - gen(k, QC1MAX)) ./ (gen(k, PC2) - gen(k, PC1));
    U(k) = Qmax_at_Pmax < gen(k, QMAX);
end

TorF = L | U;

return;
