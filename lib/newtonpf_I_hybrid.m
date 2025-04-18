function [V, converged, i] = newtonpf_I_hybrid(Ybus, Sbus, V0, ref, pv, pq, mpopt)
% newtonpf_I_hybrid - Solves power flow using full Newton's method (current/hybrid).
% ::
%
%   [V, CONVERGED, I] = NEWTONPF_I_HYBRID(YBUS, SBUS, V0, REF, PV, PQ, MPOPT)
%
%   Solves for bus voltages using a full Newton-Raphson method, using nodal
%   current balance equations and a hybrid representation of voltages, where
%   a polar update is computed using a cartesian Jacobian, given the
%   following inputs:
%       YBUS  - full system admittance matrix (for all buses)
%       SBUS  - handle to function that returns the complex bus power
%               injection vector (for all buses), given the bus voltage
%               magnitude vector (for all buses)
%       V0    - initial vector of complex bus voltages
%       REF   - bus index of reference bus (voltage ang reference & gen slack)
%       PV    - vector of bus indices for PV buses
%       PQ    - vector of bus indices for PQ buses
%       MPOPT - (optional) MATPOWER option struct, used to set the
%               termination tolerance, maximum number of iterations, and
%               output options (see MPOPTION for details).
%
%   The bus voltage vector contains the set point for generator
%   (including ref bus) buses, and the reference angle of the swing
%   bus, as well as an initial guess for remaining magnitudes and
%   angles.
%
%   Returns the final complex voltages, a flag which indicates whether it
%   converged or not, and the number of iterations performed.
%
% See also runpf, newtonpf, newtonpf_S_cart, newtonpf_I_polar.

%   MATPOWER
%   Copyright (c) 1996-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Baljinnyam Sereeter, Delft University of Technology
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% default arguments
if nargin < 7
    mpopt = mpoption;
end

%% options
tol         = mpopt.pf.tol;
max_it      = mpopt.pf.nr.max_it;
lin_solver  = mpopt.pf.nr.lin_solver;

%% initialize
converged = 0;
i = 0;
V = V0;
Va = angle(V);
Vm = abs(V);
n = length(V0);

%% set up indexing for updating V
npv = length(pv);
npq = length(pq);
j1 = 1;         j2 = npv;           %% j1:j2 - Q of pv buses
j3 = j2 + 1;    j4 = j2 + npq;      %% j3:j4 - Vr of pq buses
j5 = j4 + 1;    j6 = j4 + npv;      %% j5:j6 - Vi of pv buses
j7 = j6 + 1;    j8 = j6 + npq;      %% j7:j9 - Vi of pq buses

%% evaluate F(x0)
Sb = Sbus(Vm);
Sb(pv) = real(Sb(pv)) + 1j * imag(V(pv) .* conj(Ybus(pv, :) * V));
mis = Ybus * V - conj(Sb ./ V);
F = [   real(mis([pv; pq]));
        imag(mis([pv; pq]))   ];

%% check tolerance
normF = norm(F, inf);
if mpopt.verbose > 1
    fprintf('\n it   max Ir & Ii mismatch (p.u.)');
    fprintf('\n----  ---------------------------');
    fprintf('\n%3d        %10.3e', i, normF);
end
if normF < tol
    converged = 1;
    if mpopt.verbose > 1
        fprintf('\nConverged!\n');
    end
end

%% attempt to pick fastest linear solver, if not specified
if isempty(lin_solver)
    nx = length(F);
    if nx <= 10 || have_feature('octave')
        lin_solver = '\';       %% default \ operator
    else    %% MATLAB and nx > 10 or Octave and nx > 2000
        lin_solver = 'LU3';     %% LU decomp with 3 output args, AMD ordering
    end
end

%% do Newton iterations
while (~converged && i < max_it)
    %% update iteration counter
    i = i + 1;

    %% evaluate Jacobian
    dImis_dQ = sparse(1:n, 1:n, 1j./conj(V), n, n);
    [dImis_dVr, dImis_dVi] = dImis_dV(Sb, Ybus, V, 1);
    dImis_dVi(:, pv) = ...
        dImis_dVi(:, pv) * sparse(1:npv, 1:npv, real(V(pv)), npv, npv) - ...
        dImis_dVr(:, pv) * sparse(1:npv, 1:npv, imag(V(pv)), npv, npv);
    dImis_dVr(:, pv) = dImis_dQ(:, pv);

    %% handling of derivatives for voltage dependent loads
    %% (not yet implemented) goes here

    j11 = real(dImis_dVr([pv; pq], [pv; pq]));
    j12 = real(dImis_dVi([pv; pq], [pv; pq]));
    j21 = imag(dImis_dVr([pv; pq], [pv; pq]));
    j22 = imag(dImis_dVi([pv; pq], [pv; pq]));

    J = [   j11 j12;
            j21 j22;    ];

    %% compute update step
    dx = mplinsolve(J, -F, lin_solver);

    %% update voltage
    if npv
        Va(pv) = Va(pv) + dx(j5:j6);
        Sb(pv) = real(Sb(pv)) + 1j * (imag(Sb(pv)) + dx(j1:j2));
    end
    if npq
        Vm(pq) = Vm(pq) + (real(V(pq))./Vm(pq)) .* dx(j3:j4) ...
                        + (imag(V(pq))./Vm(pq)) .* dx(j7:j8);
        Va(pq) = Va(pq) + (real(V(pq))./(Vm(pq).^2)) .* dx(j7:j8) ...
                        - (imag(V(pq))./(Vm(pq).^2)) .* dx(j3:j4);
    end
    V = Vm .* exp(1j * Va);
    Vm = abs(V);            %% update Vm and Va again in case
    Va = angle(V);          %% we wrapped around with a negative Vm

    %% evalute F(x)
    mis = Ybus * V - conj(Sb ./ V);
    F = [   real(mis([pv; pq]));
            imag(mis([pv; pq]))   ];

    %% check for convergence
    normF = norm(F, inf);
    if mpopt.verbose > 1
        fprintf('\n%3d        %10.3e', i, normF);
    end
    if normF < tol
        converged = 1;
        if mpopt.verbose
            fprintf('\nNewton''s method power flow (current balance, hybrid) converged in %d iterations.\n', i);
        end
    end
end

if mpopt.verbose
    if ~converged
        fprintf('\nNewton''s method power flow (current balance, hybrid) did not converge in %d iterations.\n', i);
    end
end
