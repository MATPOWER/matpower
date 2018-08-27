function [V, converged, i] = newtonpf_PC(Ybus, Sbus, V0, ref, pv, pq, mpopt)
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
j1 = 1;         j2 = npv;           %% j1:j2 - Vi of pv buses
j3 = j2 + 1;    j4 = j2 + npq;      %% j3:j4 - Vi of pq buses
j5 = j4 + 1;    j6 = j4 + npq;      %% j5:j6 - Vr of pq buses
%% evaluate F(x0)
mis = V .* conj(Ybus * V) - Sbus(Vm);
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
    [dSbus_dVr, dSbus_dVi] = dSbus_dV(Ybus, V, 1);
        dSbus_dVi(:,pv) = dSbus_dVi(:,pv) - dSbus_dVr(:,pv)*sparse(1:npv, 1:npv, imag(V(pv))./real(V(pv)), npv, npv);
%    [dummy, neg_dSd_dVm] = Sbus(Vm);
%    dSbus_dVm = dSbus_dVm - neg_dSd_dVm;
    j11 = real(dSbus_dVi([pv; pq], [pv; pq]));
    j12 = real(dSbus_dVr([pv; pq], pq));
    j21 = imag(dSbus_dVi(pq, [pv; pq]));
    j22 = imag(dSbus_dVr(pq, pq));
    
    J = [   j11 j12;
            j21 j22;    ];
        
    %% compute update step
    dx = mplinsolve(J, -F, lin_solver);

    %% update voltage
    if npv
        Va(pv) = Va(pv) + dx(j1:j2)./real(V(pv));
    end
    if npq
        Vm(pq) = Vm(pq) + (real(V(pq))./Vm(pq)).*dx(j5:j6) + (imag(V(pq))./Vm(pq)).*dx(j3:j4);           
        Va(pq) = Va(pq) + (real(V(pq))./(Vm(pq).^2)).*dx(j3:j4) - (imag(V(pq))./(Vm(pq).^2)).*dx(j5:j6);
    end
    V = Vm .* exp(1j * Va);
    Vm = abs(V);            %% update Vm and Va again in case
    Va = angle(V);          %% we wrapped around with a negative Vm

    %% evalute F(x)
    mis = V .* conj(Ybus * V) - Sbus(Vm);
    F = [   real(mis([pv; pq]));
            imag(mis(pq))   ];

    %% check for convergence
    normF = norm(F, inf);
    if mpopt.verbose > 1
        fprintf('\n%3d        %10.3e', i, normF);
    end
    if normF < tol
        converged = 1;
        if mpopt.verbose
            fprintf('\nNewton''s method power flow using power balance in Cartesian coordinates converged in %d iterations.\n', i);
        end
    end
end

if mpopt.verbose
    if ~converged
        fprintf('\nNewton''s method power flow using power balance in Cartesian coordinates did not converge in %d iterations.\n', i);
    end
end
