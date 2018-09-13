function [V, converged, i] = newtonpf_I_cart(Ybus, Sbus, V0, ref, pv, pq, mpopt)
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
j1 = 1;         j2 = npv;           %% j1:j2 - Vi of pv buses
j3 = j2 + 1;    j4 = j2 + npq;      %% j3:j4 - Vi of pq buses
j5 = j4 + 1;    j6 = j4 + npv;      %% j5:j6 - Q of pv buses
j7 = j6 + 1;    j8 = j6 + npq;      %% j7:j8 - Vr of pq buses
%% evaluate F(x0)
Sbus1 = Sbus(Vm);
Sbus1(pv) = real(Sbus1(pv)) + 1j*imag(V(pv).* conj(Ybus(pv,:)*V));
mis = Ybus*V - conj(Sbus1./V);
F = [   real(mis([pv; pq]));
        imag(mis([pv; pq]))   ];
%% check tolerance
normF = norm(F, inf);
if mpopt.verbose > 1
    fprintf('\n it    max Ir & Ii mismatch (p.u.)');
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
    dImis_dQ  = sparse(1:n, 1:n, 1j./conj(V), n, n); 
    [dImis_dVr, dImis_dVi] = dImis_dV(Sbus1, Ybus, V, 1);
        dImis_dVi(:,pv) = dImis_dVi(:,pv) - dImis_dVr(:,pv)*sparse(1:npv, 1:npv, imag(V(pv))./real(V(pv)), npv, npv);  
        dImis_dVr(:,pv) = dImis_dQ(:,pv);
%    [dummy, neg_dSd_dVm] = Sbus(Vm);
%    dSbus_dVm = dSbus_dVm - neg_dSd_dVm;
    j11 = real(dImis_dVi([pv; pq], [pv; pq]));
    j12 = real(dImis_dVr([pv; pq], [pv; pq]));
    j21 = imag(dImis_dVi([pv; pq], [pv; pq]));
    j22 = imag(dImis_dVr([pv; pq], [pv; pq]));
    
    J = [   j11 j12;
            j21 j22;    ];
        
    %% compute update step
    dx = mplinsolve(J, -F, lin_solver);

    %% update voltage
    if npv
        Va(pv) = Va(pv) + dx(j1:j2)./real(V(pv));
        Sbus1(pv) = real(Sbus1(pv)) + 1j*(imag(Sbus1(pv)) + dx(j5:j6));
    end
    if npq
        Vm(pq) = Vm(pq) + (real(V(pq))./Vm(pq)).*dx(j7:j8) + (imag(V(pq))./Vm(pq)).*dx(j3:j4);
        Va(pq) = Va(pq) + (real(V(pq))./(Vm(pq).^2)).*dx(j3:j4) - (imag(V(pq))./(Vm(pq).^2)).*dx(j7:j8);
    end
    V = Vm .* exp(1j * Va);
    Vm = abs(V);            %% update Vm and Va again in case
    Va = angle(V);          %% we wrapped around with a negative Vm

    %% evalute F(x)
    mis = Ybus*V - conj(Sbus1./V);
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
            fprintf('\nNewton''s method power flow using current balance in Cartesian coordinates converged in %d iterations.\n', i);
        end
    end
end

if mpopt.verbose
    if ~converged
        fprintf('\nNewton''s method power flow using current balance in Cartesianr coordinates did not converge in %d iterations.\n', i);
    end
end
