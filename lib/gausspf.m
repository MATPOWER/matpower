function [V, converged, i] = gausspf(Ybus, Sbus, V0, ref, pv, pq, mpopt)
%GAUSSPF  Solves the power flow using a Gauss-Seidel method.
%   [V, converged, i] = gausspf(Ybus, Sbus, V0, ref, pv, pq, mpopt)
%   solves for bus voltages given the full system admittance matrix (for
%   all buses), the complex bus power injection vector (for all buses),
%   the initial vector of complex bus voltages, and column vectors with
%   the lists of bus indices for the swing bus, PV buses, and PQ buses,
%   respectively. The bus voltage vector contains the set point for
%   generator (including ref bus) buses, and the reference angle of the
%   swing bus, as well as an initial guess for remaining magnitudes and
%   angles. mpopt is a MATPOWER options vector which can be used to 
%   set the termination tolerance, maximum number of iterations, and 
%   output options (see 'help mpoption' for details). Uses default
%   options if this parameter is not given. Returns the final complex
%   voltages, a flag which indicates whether it converged or not, and
%   the number of iterations performed.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   and Alberto Borghetti, University of Bologna, Italy
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% default arguments
if nargin < 7
    mpopt = mpoption;
end

%% options
tol     = mpopt(2);
max_it  = mpopt(5);
verbose = mpopt(31);

%% initialize
j = sqrt(-1);
converged = 0;
i = 0;
V = V0;
Va = angle(V);
Vm = abs(V);

%% set up indexing for updating V
npv = length(pv);
npq = length(pq);
j1 = 1;         j2 = npv;           %% j1:j2 - V angle of pv buses
j3 = j2 + 1;    j4 = j2 + npq;      %% j3:j4 - V angle of pq buses
j5 = j4 + 1;    j6 = j4 + npq;      %% j5:j6 - V mag of pq buses

%% evaluate F(x0)
mis = V .* conj(Ybus * V) - Sbus;
F = [   real(mis([pv; pq]));
        imag(mis(pq))   ];

%% check tolerance
normF = norm(F, inf);
if verbose > 1
    fprintf('\n it    max P & Q mismatch (p.u.)');
    fprintf('\n----  ---------------------------');
    fprintf('\n%3d        %10.3e', i, normF);
end
if normF < tol
    converged = 1;
    if verbose > 1
        fprintf('\nConverged!\n');
    end
end

%% do Gauss-Seidel iterations
while (~converged && i < max_it)
    %% update iteration counter
    i = i + 1;

    %% update voltage
    %% at PQ buses
    for k = pq(1:npq)'
        V(k) =  V(k) + (conj(Sbus(k) / V(k)) - Ybus(k,:) * V ) / Ybus(k,k);
    end

    %% at PV buses
    if npv
        for k = pv(1:npv)'
            Sbus(k) = real(Sbus(k)) + j * imag( V(k) .* conj(Ybus(k,:) * V));
            V(k) =  V(k) + (conj(Sbus(k) / V(k)) - Ybus(k,:) * V ) / Ybus(k,k);
%           V(k) = Vm(k) * V(k) / abs(V(k));
        end
        V(pv) = Vm(pv) .* V(pv) ./ abs(V(pv));
    end

    %% evalute F(x)
    mis = V .* conj(Ybus * V) - Sbus;
    F = [   real(mis(pv));
            real(mis(pq));
            imag(mis(pq))   ];

    %% check for convergence
    normF = norm(F, inf);
    if verbose > 1
        fprintf('\n%3d        %10.3e', i, normF);
    end
    if normF < tol
        converged = 1;
        if verbose
            fprintf('\nGauss-Seidel power flow converged in %d iterations.\n', i);
        end
    end
end

if verbose
    if ~converged
        fprintf('\nGauss-Seidel power did not converge in %d iterations.\n', i);
    end
end
