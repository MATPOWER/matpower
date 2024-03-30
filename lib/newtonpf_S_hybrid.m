function [V, converged, i] = newtonpf_S_hybrid(Ybus, Sbus, V0, ref, pv, pq, mpopt)
% newtonpf_S_hybrid - Solves power flow using full Newton's method (power/hybrid).
% ::
%
%   [V, CONVERGED, I] = NEWTONPF_S_HYBRID(YBUS, SBUS, V0, REF, PV, PQ, MPOPT)
%
%   Solves for bus voltages using a full Newton-Raphson method, using nodal
%   power balance equations and a hybrid representation of voltages, where
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
% See also runpf, newtonpf, newtonpf_I_polar, newtonpf_I_cart.

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

%% set up indexing for updating V
npv = length(pv);
npq = length(pq);
j1 = 1;         j2 = npq;           %% j1:j2 - Vr of pq buses
j3 = j2 + 1;    j4 = j2 + npv;      %% j3:j4 - Vi of pv buses
j5 = j4 + 1;    j6 = j4 + npq;      %% j5:j6 - Vi of pq buses

%% evaluate F(x0)
mis = V .* conj(Ybus * V) - Sbus(Vm);
F = [   real(mis([pq; pv]));
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
    [dSbus_dVr, dSbus_dVi] = dSbus_dV(Ybus, V, 1);
    dSbus_dVi(:, pv) = ...
        dSbus_dVi(:, pv) * sparse(1:npv, 1:npv, real(V(pv)), npv, npv) - ...
        dSbus_dVr(:, pv) * sparse(1:npv, 1:npv, imag(V(pv)), npv, npv);

    %% handling of derivatives for voltage dependent loads
    %% (not yet implemented) goes here

    j11 = real(dSbus_dVr([pq; pv], pq));
    j12 = real(dSbus_dVi([pq; pv], [pv; pq]));
    j21 = imag(dSbus_dVr(pq, pq));
    j22 = imag(dSbus_dVi(pq, [pv; pq]));

    J = [   j11 j12;
            j21 j22;    ];

    %% compute update step
    dx = mplinsolve(J, -F, lin_solver);

    %% update voltage
    if npv
        Va(pv) = Va(pv) + dx(j3:j4);
    end
    if npq
        Vm(pq) = Vm(pq) + (real(V(pq))./Vm(pq)) .* dx(j1:j2) ...
                        + (imag(V(pq))./Vm(pq)) .* dx(j5:j6);
        Va(pq) = Va(pq) + (real(V(pq))./(Vm(pq).^2)) .* dx(j5:j6) ...
                        - (imag(V(pq))./(Vm(pq).^2)) .* dx(j1:j2);
    end
    V = Vm .* exp(1j * Va);
    Vm = abs(V);            %% update Vm and Va again in case
    Va = angle(V);          %% we wrapped around with a negative Vm

    %% evalute F(x)
    mis = V .* conj(Ybus * V) - Sbus(Vm);
    F = [   real(mis([pq; pv]));
            imag(mis(pq))   ];

    %% check for convergence
    normF = norm(F, inf);
    if mpopt.verbose > 1
        fprintf('\n%3d        %10.3e', i, normF);
    end
    if normF < tol
        converged = 1;
        if mpopt.verbose
            fprintf('\nNewton''s method power flow (power balance, hybrid) converged in %d iterations.\n', i);
        end
    end
end

if mpopt.verbose
    if ~converged
        fprintf('\nNewton''s method power flow (power balance, hybrid) did not converge in %d iterations.\n', i);
    end
end
