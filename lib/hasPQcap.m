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
%   Copyright (c) 2005-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% default value
if nargin < 2
    hilo = 'B';     %% look at both top and bottom by default
end

%% for which gens is it specified
k = find( gen(:, PC1) | gen(:, PC2) );
ng = size(gen, 1);

if isempty(k)
    TorF = zeros(ng, 1);
else
    %% eliminate cases where QMIN = QMAX = QC
    kk = find( gen(k, QMIN) == gen(k, QMAX) & ...
                gen(k, QMIN) == gen(k, QC1MAX) & ...
                gen(k, QMIN) == gen(k, QC1MIN) & ...
                gen(k, QMIN) == gen(k, QC2MAX) & ...
                gen(k, QMIN) == gen(k, QC2MIN) );
    k(kk) = [];

    %% check for errors in capability curve data
    if any( gen(k, PC1) >= gen(k, PC2) )
        error('hasPQcap: must have Pc1 < Pc2');
    end
    if any( gen(k, QC2MAX) <= gen(k, QC2MIN) & gen(k, QC1MAX) <= gen(k, QC1MIN) )
        error('hasPQcap: capability curve defines an empty set');
    end

    %% for which gens is it specified
    k = find( gen(:, PC1) ~= gen(:, PC2) );
    L = zeros(ng, 1);
    U = zeros(ng, 1);
    dPc = gen(k, PC2) - gen(k, PC1);

    if ~strcmp(hilo, 'U')       %% include lower constraint
        dQc = gen(k, QC2MIN) - gen(k, QC1MIN);
        Qmin_at_Pmin = gen(k, QC1MIN) + (gen(k, PMIN) - gen(k, PC1)) .* ...
            dQc ./ dPc;
        Qmin_at_Pmax = gen(k, QC1MIN) + (gen(k, PMAX) - gen(k, PC1)) .* ...
            dQc ./ dPc;
        L(k) = Qmin_at_Pmin > gen(k, QMIN) | Qmin_at_Pmax > gen(k, QMIN);
    end

    if ~strcmp(hilo, 'L')       %% include upper constraint
        dQc = gen(k, QC2MAX) - gen(k, QC1MAX);
        Qmax_at_Pmin = gen(k, QC1MAX) + (gen(k, PMIN) - gen(k, PC1)) .* ...
            dQc ./ dPc;
        Qmax_at_Pmax = gen(k, QC1MAX) + (gen(k, PMAX) - gen(k, PC1)) .* ...
            dQc ./ dPc;
        U(k) = Qmax_at_Pmin < gen(k, QMAX) | Qmax_at_Pmax < gen(k, QMAX);
    end

    TorF = L | U;
end
