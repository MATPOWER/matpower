function [V, converged, i, lam] = cpf_corrector(Ybus, Sbus, V0, ref, pv, pq, ...
                                  lam0, Sxfr, Vprv, lamprv, z, step, mpopt)
%CPF_CORRECTOR  Solves the corrector step of a continuation power flow using a
%   full Newton method with pseudo-arclength parameterization.
%   [V, CONVERGED, I, LAM] = CPF_CORRECTOR(YBUS, SBUS, V0, REF, PV, PQ, ...
%                            LAM0, SXFR, VPRV, LPRV, Z, STEP, MPOPT)
%   solves for bus voltages and lambda given the full system admittance
%   matrix (for all buses), the complex bus power injection vector (for
%   all buses), the initial vector of complex bus voltages, and column
%   vectors with the lists of bus indices for the swing bus, PV buses, and
%   PQ buses, respectively. The bus voltage vector contains the set point
%   for generator (including ref bus) buses, and the reference angle of the
%   swing bus, as well as an initial guess for remaining magnitudes and
%   angles. MPOPT is a MATPOWER options vector which can be used to
%   set the termination tolerance, maximum number of iterations, and
%   output options (see MPOPTION for details). Uses default options if
%   this parameter is not given. Returns the final complex voltages, a
%   flag which indicates whether it converged or not, the number
%   of iterations performed, and the final lambda.
%
%   The extra continuation inputs are LAM0 (initial predicted lambda),
%   SXFR ([delP+j*delQ] transfer/loading vector for all buses), VPRV
%   (final complex V corrector solution from previous continuation step),
%   LAMPRV (final lambda corrector solution from previous continuation step),
%   Z (normalized predictor for all buses), and STEP (continuation step size).
%   The extra continuation output is LAM (final corrector lambda).
%
%   See also RUNCPF.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell,
%   Shrirang Abhyankar, Argonne National Laboratory,
%   and Alexander Flueck, IIT
%   Copyright (c) 1996-2013 by Power System Engineering Research Center (PSERC)
%
%   Modified by Alexander J. Flueck, Illinois Institute of Technology
%   2001.02.22 - corrector.m (ver 1.0) based on newtonpf.m (MATPOWER 2.0)
%
%   Modified by Shrirang Abhyankar, Argonne National Laboratory
%   Updated to be compatible with MATPOWER version 4.1)
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
if nargin < 13
    mpopt = mpoption;
end

%% options
tol     = mpopt(2);
max_it  = mpopt(3);
verbose = mpopt(31);

%% initialize
converged = 0;
i = 0;
V = V0;
Va = angle(V);
Vm = abs(V);
lam = lam0;             %% set lam to initial lam0
Vprva = angle(Vprv);    %% previous step's voltage angles
Vprvm = abs(Vprv);      %% previous step's voltage magnitudes

%% set up indexing for updating V
npv = length(pv);
npq = length(pq);
nb = length(V);         %% number of buses
j1 = 1;         j2 = npv;           %% j1:j2 - V angle of pv buses
j3 = j2 + 1;    j4 = j2 + npq;      %% j3:j4 - V angle of pq buses
j5 = j4 + 1;    j6 = j4 + npq;      %% j5:j6 - V mag of pq buses
j7 = j6 + 1;    j8 = j6 + 1;        %% j7:j8 - lambda

%% evaluate F(x0, lam0), including Sxfr transfer/loading
mis = V .* conj(Ybus * V) - Sbus - lam*Sxfr;
F = [   real(mis([pv; pq]));
        imag(mis(pq))   ];

%% evaluate P(x0, lambda0)
P = z([pv; pq; nb+pq; 2*nb+1])' * ...
    ( [Va([pv; pq]); Vm(pq); lam] - [Vprva([pv; pq]); Vprvm(pq); lamprv] )...
    - step;

%% augment F(x,lambda) with P(x,lambda)
F = [ F; 
      P ];

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

%% do Newton iterations
while (~converged && i < max_it)
    %% update iteration counter
    i = i + 1;
    
    %% evaluate Jacobian
    [dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, V);
    
    j11 = real(dSbus_dVa([pv; pq], [pv; pq]));
    j12 = real(dSbus_dVm([pv; pq], pq));
    j21 = imag(dSbus_dVa(pq, [pv; pq]));
    j22 = imag(dSbus_dVm(pq, pq));
    
    J = [   j11 j12;
            j21 j22;    ];

    %% augment J with real/imag -Sxfr and z^T
    J = [ J, -[real(Sxfr([pv; pq])); imag(Sxfr(pq))]; 
                  z([pv; pq; nb+pq; 2*nb+1])'        ];

    %% compute update step
    dx = -(J \ F);

    %% update voltage
    if npv
        Va(pv) = Va(pv) + dx(j1:j2);
    end
    if npq
        Va(pq) = Va(pq) + dx(j3:j4);
        Vm(pq) = Vm(pq) + dx(j5:j6);
    end
    V = Vm .* exp(1j * Va);
    Vm = abs(V);            %% update Vm and Va again in case
    Va = angle(V);          %% we wrapped around with a negative Vm

    %% update lambda
    lam = lam + dx(j7:j8);

    %% evalute F(x, lam)
    mis = V .* conj(Ybus * V) - Sbus - lam*Sxfr;
    F = [   real(mis(pv));
            real(mis(pq));
            imag(mis(pq))   ];

    %% evaluate P(x, lambda)
    P = z([pv; pq; nb+pq; 2*nb+1])' * ...
        ([Va([pv; pq]); Vm(pq); lam] - [Vprva([pv; pq]); Vprvm(pq); lamprv]) ...
        - step;

    %% augment F(x,lambda) with P(x,lambda)
    F = [ F; 
          P ];

    %% check for convergence
    normF = norm(F, inf);
    if verbose > 1
        fprintf('\n%3d        %10.3e', i, normF);
    end
    if normF < tol
        converged = 1;
        if verbose
            fprintf('\nNewton''s method corrector converged in %d iterations.\n', i);
        end
    end
end

if verbose
    if ~converged
        fprintf('\nNewton''s method corrector did not converge in %d iterations.\n', i);
    end
end
