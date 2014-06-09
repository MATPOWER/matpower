function TorF = hasPQcap(gen, hilo)
%HASPQCAP  Checks for P-Q capability curve constraints.
%   TORF = HASPQCAP(GEN, HILO) returns a column vector of 1's and 0's. The 1's
%   correspond to rows of the GEN matrix which correspond to generators which
%   have defined a capability curve (with sloped upper and/or lower bound on
%   Q) and require that additional linear constraints be added to the OPF.
%
%   The GEN matrix in version 2 of the MATPOWER case format includes columns
%   for specifying a P-Q capability curve for a generator defined as the
%   intersection of two half-planes and the box constraints on P and Q. The
%   two half planes are defined respectively as the area below the line
%   connecting (Pc1, Qc1max) and (Pc2, Qc2max) and the area above the line
%   connecting (Pc1, Qc1min) and (Pc2, Qc2min).
%
%   If the optional 2nd argument is 'U' this function returns true only for
%   rows corresponding to generators that require the upper constraint on Q.
%   If it is 'L', only for those requiring the lower constraint. If the 2nd
%   argument is not specified or has any other value it returns true for rows
%   corresponding to gens that require either or both of the constraints.
%
%   It is smart enough to return true only if the corresponding linear
%   constraint is not redundant w.r.t the box constraints.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2005-2014 by Power System Engineering Research Center (PSERC)
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

[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% default value
if nargin < 2
    hilo = 'B';     %% look at both top and bottom by default
end

%% check for errors capability curve data
if any( gen(:, PC1) > gen(:, PC2) )
    error('hasPQcap: must have Pc1 < Pc2');
end
if any( gen(:, QC2MAX) < gen(:, QC2MIN) & gen(:, QC1MAX) < gen(:, QC1MIN) )
    error('hasPQcap: capability curve defines and empty set');
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
