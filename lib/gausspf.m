function [V, converged, i] = gausspf(Ybus, Sbus, V0, ref, pv, pq, mpopt)
%GAUSSPF  Solves the power flow using a Gauss-Seidel method.
%   [V, CONVERGED, I] = GAUSSPF(YBUS, SBUS, V0, REF, PV, PQ, MPOPT)
%   solves for bus voltages given the full system admittance matrix (for
%   all buses), the complex bus power injection vector (for all buses),
%   the initial vector of complex bus voltages, and column vectors with
%   the lists of bus indices for the swing bus, PV buses, and PQ buses,
%   respectively. The bus voltage vector contains the set point for
%   generator (including ref bus) buses, and the reference angle of the
%   swing bus, as well as an initial guess for remaining magnitudes and
%   angles. MPOPT is a MATPOWER options struct which can be used to 
%   set the termination tolerance, maximum number of iterations, and 
%   output options (see MPOPTION for details). Uses default options
%   if this parameter is not given. Returns the final complex voltages,
%   a flag which indicates whether it converged or not, and the number
%   of iterations performed.
%
%   See also RUNPF.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   and Alberto Borghetti, University of Bologna, Italy
%   Copyright (c) 1996-2011 by Power System Engineering Research Center (PSERC)
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

%% default arguments
if nargin < 7
    mpopt = mpoption;
end

%% options
tol     = mpopt.pf.tol;
max_it  = mpopt.pf.gs.max_it;

%% initialize
converged = 0;
i = 0;
V = V0;
Vm = abs(V);

%% set up indexing for updating V
npv = length(pv);
npq = length(pq);

%% evaluate F(x0)
mis = V .* conj(Ybus * V) - Sbus;
F = [   real(mis([pv; pq]));
        imag(mis(pq))   ];

%% check tolerance
normF = norm(F, inf);
if mpopt.verbose > 1
    fprintf('\n it    max P & Q mismatch (p.u.)');
    fprintf('\n----  ---------------------------');
    fprintf('\n%3d        %10.3e', i, normF);
end
if normF < tol
    converged = 1;
    if mpopt.verbose > 1
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
            Sbus(k) = real(Sbus(k)) + 1j * imag( V(k) .* conj(Ybus(k,:) * V));
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
    if mpopt.verbose > 1
        fprintf('\n%3d        %10.3e', i, normF);
    end
    if normF < tol
        converged = 1;
        if mpopt.verbose
            fprintf('\nGauss-Seidel power flow converged in %d iterations.\n', i);
        end
    end
end

if mpopt.verbose
    if ~converged
        fprintf('\nGauss-Seidel power flow did not converge in %d iterations.\n', i);
    end
end
