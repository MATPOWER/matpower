function [V, converged, i] = newtonpf(Ybus, Sbus, V0, ref, pv, pq, mpopt)
%NEWTONPF  Solves the power flow using a full Newton's method.
%   [V, CONVERGED, I] = NEWTONPF(YBUS, SBUS, V0, REF, PV, PQ, MPOPT)
%
%   Solves for bus voltages using a full Newton-Raphson method, given the
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
%   See also RUNPF.

%   MATPOWER
%   Copyright (c) 1996-2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% default arguments
if nargin < 7
    mpopt = mpoption;
end

%% options
tol         = mpopt.pf.tol;
max_it      = mpopt.pf.nr.max_it;
lin_solver  = mpopt.pf.nr.lin_solver;
v_cartesian = mpopt.pf.v_cartesian;
current_balance = mpopt.pf.current_balance;

if current_balance
    variants = 'current balance in';
else
    variants = 'power balance in';
end
if v_cartesian
    variants = strcat(variants,' Cartesian coordinates');
else
    variants = strcat(variants,' Polar coordinates');
end
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
j1 = 1;         j2 = npv;               %% j1:j2 - Va/Vi of pv buses
j3 = j2 + 1;    j4 = j2 + npq;          %% j3:j4 - Va/Vi of pq buses
if current_balance
    j7 = j4 + 1;    j8 = j4 + npv;      %% j7:j8 - Q of pv buses
    j5 = j8 + 1;    j6 = j8 + npq;      %% j5:j6 - Vm/Vr of pq buses    
else
    j5 = j4 + 1;    j6 = j4 + npq;      %% j5:j6 - Vm/Vr of pq buses
end
%% evaluate F(x0)
if current_balance
    Sbus1 = Sbus(Vm);
    Sbus1(pv) = real(Sbus1(pv)) + 1j*imag(V(pv).* conj(Ybus(pv,:)*V));  
    mis = Ybus*V - conj(Sbus1./V);
    F = [   real(mis([pv; pq]));
            imag(mis([pv; pq]))   ];
else
    mis = V .* conj(Ybus * V) - Sbus(Vm);
    F = [   real(mis([pv; pq]));
            imag(mis(pq))   ];
end

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
    if nx <= 10 || have_fcn('octave')
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
    if current_balance    
        dImis_dQ  = sparse(1:n, 1:n, 1j./conj(V), n, n); 
        [dImis_dV1, dImis_dV2] = dImis_dV(Sbus1, Ybus, V, mpopt.pf.v_cartesian);
        if v_cartesian
            dImis_dV2(:,pv) = dImis_dV2(:,pv) - dImis_dV1(:,pv)*sparse(1:npv, 1:npv, imag(V(pv))./real(V(pv)), npv, npv);  
            dImis_dV1(:,pv) = dImis_dQ(:,pv);        
            j11 = real(dImis_dV2([pv; pq], [pv; pq]));
            j12 = real(dImis_dV1([pv; pq], [pv; pq]));
            j21 = imag(dImis_dV2([pv; pq], [pv; pq]));
            j22 = imag(dImis_dV1([pv; pq], [pv; pq]));             
        else
            dImis_dV2(:,pv) = dImis_dQ(:,pv); 
            j11 = real(dImis_dV1([pv; pq], [pv; pq]));
            j12 = real(dImis_dV2([pv; pq], [pv; pq]));
            j21 = imag(dImis_dV1([pv; pq], [pv; pq]));
            j22 = imag(dImis_dV2([pv; pq], [pv; pq]));            
        end
    else
        [dSbus_dV1, dSbus_dV2] = dSbus_dV(Ybus, V, mpopt.pf.v_cartesian);        
        if v_cartesian
            dSbus_dV2(:,pv) = dSbus_dV2(:,pv) - dSbus_dV1(:,pv)*sparse(1:npv, 1:npv, imag(V(pv))./real(V(pv)), npv, npv);
            j11 = real(dSbus_dV2([pv; pq], [pv; pq]));
            j12 = real(dSbus_dV1([pv; pq], pq));
            j21 = imag(dSbus_dV2(pq, [pv; pq]));
            j22 = imag(dSbus_dV1(pq, pq));            
        else
            [dummy, neg_dSd_dVm] = Sbus(Vm);
            dSbus_dV2 = dSbus_dV2 - neg_dSd_dVm;            
            j11 = real(dSbus_dV1([pv; pq], [pv; pq]));
            j12 = real(dSbus_dV2([pv; pq], pq));
            j21 = imag(dSbus_dV1(pq, [pv; pq]));
            j22 = imag(dSbus_dV2(pq, pq));            
        end
    end
    J = [   j11 j12;
            j21 j22;    ];

    %% compute update step
    dx = mplinsolve(J, -F, lin_solver);

    %% update voltage
    if npv
        if v_cartesian
            Va(pv) = Va (pv) + dx(j1:j2)./real(V(pv));
        else
            Va(pv) = Va(pv) + dx(j1:j2);            
        end
        if current_balance
            Sbus1(pv) = real(Sbus1(pv)) + 1j*(imag(Sbus1(pv)) + dx(j7:j8));
        end
    end
    if npq
        if v_cartesian
            Vm(pq) = Vm(pq) + (real(V(pq))./Vm(pq)).*dx(j5:j6) + (imag(V(pq))./Vm(pq)).*dx(j3:j4);           
            Va(pq) = Va(pq) + (real(V(pq))./(Vm(pq).^2)).*dx(j3:j4) - (imag(V(pq))./(Vm(pq).^2)).*dx(j5:j6);                 
        else
            Va(pq) = Va(pq) + dx(j3:j4);
            Vm(pq) = Vm(pq) + dx(j5:j6);                            
        end
    end
    V = Vm .* exp(1j * Va);
    Vm = abs(V);            %% update Vm and Va again in case
    Va = angle(V);          %% we wrapped around with a negative Vm

    %% evalute F(x)
    if current_balance
        mis = Ybus*V - conj(Sbus1./V);
        F = [   real(mis([pv; pq]));
                imag(mis([pv; pq]))   ];
    else
        mis = V .* conj(Ybus * V) - Sbus(Vm);
        F = [   real(mis([pv; pq]));
                imag(mis(pq))   ];
    end

    %% check for convergence
    normF = norm(F, inf);
    if mpopt.verbose > 1
        fprintf('\n%3d        %10.3e', i, normF);
    end
    if normF < tol
        converged = 1;
        if mpopt.verbose
            fprintf('\nNewton''s method power flow using %s converged in %d iterations.\n', variants, i);
        end
    end
end

if mpopt.verbose
    if ~converged
        fprintf('\nNewton''s method power flow using %s did not converge in %d iterations.\n', variants, i);
    end
end
