function [V, converged, i, lam] = cpf_corrector(Ybus, Sbusb, V_hat, ref, pv, pq, ...
                lam_hat, Sbust, Vprv, lamprv, z, step, parameterization, mpopt)
%CPF_CORRECTOR  Solves the corrector step of a continuation power flow
%   [V, CONVERGED, I, LAM] = CPF_CORRECTOR(YBUS, SBUSB, V_HAT, REF, PV, PQ, ...
%                                       LAM_HAT, SBUST, VPRV, LPRV, Z, ...
%                                       STEP, PARAMETERIZATION, MPOPT)
%
%   Computes the corrector step of a continuation power flow using a
%   full Newton method with selected parameterization scheme.
%
%   Inputs:
%       YBUS : complex bus admittance matrix
%       SBUSB : handle of function returning nb x 1 vector of complex
%               base case injections in p.u. and derivatives w.r.t. |V|
%       V_HAT :  predicted complex bus voltage vector
%       REF : vector of indices for REF buses
%       PV : vector of indices of PV buses
%       PQ : vector of indices of PQ buses
%       LAM_HAT : predicted scalar lambda
%       SBUST : handle of function returning nb x 1 vector of complex
%               target case injections in p.u. and derivatives w.r.t. |V|
%       VPRV : complex bus voltage vector at previous solution
%       LAMPRV : scalar lambda value at previous solution
%       STEP : continuation step length
%       Z : normalized tangent prediction vector
%       STEP : continuation step size
%       PARAMETERIZATION : Value of cpf.parameterization option.
%       MPOPT : Options struct
%
%   Outputs:
%       V : complex bus voltage solution vector
%       CONVERGED : Newton iteration count
%       I : Newton iteration count
%       LAM : lambda continuation parameter
%
%   See also RUNCPF.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell,
%   Shrirang Abhyankar, Argonne National Laboratory,
%   and Alexander Flueck, IIT
%
%   Modified by Alexander J. Flueck, Illinois Institute of Technology
%   2001.02.22 - corrector.m (ver 1.0) based on newtonpf.m (MATPOWER 2.0)
%
%   Modified by Shrirang Abhyankar, Argonne National Laboratory
%   (Updated to be compatible with MATPOWER version 4.1)
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% default arguments
if nargin < 14
    mpopt = mpoption;
end

%% options
tol     = mpopt.pf.tol;
max_it  = mpopt.pf.nr.max_it;

%% initialize
converged = 0;
i = 0;
V = V_hat;          %% initialize V with predicted V
Va = angle(V);
Vm = abs(V);
lam = lam_hat;      %% initialize lam with predicted lam

%% set up indexing for updating V
npv = length(pv);
npq = length(pq);
nb = length(V);         %% number of buses
j1 = 1;         j2 = npv;           %% j1:j2 - V angle of pv buses
j3 = j2 + 1;    j4 = j2 + npq;      %% j3:j4 - V angle of pq buses
j5 = j4 + 1;    j6 = j4 + npq;      %% j5:j6 - V mag of pq buses
j7 = j6 + 1;    j8 = j6 + 1;        %% j7:j8 - lambda

%% evaluate F(x0, lam_hat), including transfer/loading
Sb = Sbusb(Vm);
St = Sbust(Vm);
mis = V .* conj(Ybus * V) - Sb - lam * (St - Sb);
F = [   real(mis([pv; pq]));
        imag(mis(pq))   ];

%% evaluate P(x0, lambda0)
P = cpf_p(parameterization, step, z, V, lam, Vprv, lamprv, pv, pq);

%% augment F(x,lambda) with P(x,lambda)
F = [ F; 
      P ];

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

%% do Newton iterations
while (~converged && i < max_it)
    %% update iteration counter
    i = i + 1;
    
    %% evaluate Jacobian
    [dSbus_dVa, dSbus_dVm] = dSbus_dV(Ybus, V);
    [dummy, neg_dSdb_dVm] = Sbusb(Vm);
    [dummy, neg_dSdt_dVm] = Sbust(Vm);
    dSbus_dVm = dSbus_dVm - neg_dSdb_dVm - lam * (neg_dSdt_dVm - neg_dSdb_dVm);
    
    j11 = real(dSbus_dVa([pv; pq], [pv; pq]));
    j12 = real(dSbus_dVm([pv; pq], pq));
    j21 = imag(dSbus_dVa(pq, [pv; pq]));
    j22 = imag(dSbus_dVm(pq, pq));
    
    J = [   j11 j12;
            j21 j22;    ];

    Sxf = St - Sb;
    dF_dlam = -[real(Sxf([pv; pq])); imag(Sxf(pq))];
    [dP_dV, dP_dlam] = cpf_p_jac(parameterization, z, V, lam, Vprv, lamprv, pv, pq);

    %% augment J with real/imag -Sxfr and z^T
    J = [   J   dF_dlam; 
          dP_dV dP_dlam ];

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
    Sb = Sbusb(Vm);
    St = Sbust(Vm);
    mis = V .* conj(Ybus * V) - Sb - lam * (St - Sb);
    F = [   real(mis(pv));
            real(mis(pq));
            imag(mis(pq))   ];

    %% evaluate P(x, lambda)
    P = cpf_p(parameterization, step, z, V, lam, Vprv, lamprv, pv, pq);

    %% augment F(x,lambda) with P(x,lambda)
    F = [ F; 
          P ];

    %% check for convergence
    normF = norm(F, inf);
    if mpopt.verbose > 1
        fprintf('\n%3d        %10.3e', i, normF);
    end
    if normF < tol
        converged = 1;
        if mpopt.verbose
            fprintf('\nNewton''s method corrector converged in %d iterations.\n', i);
        end
    end
end

if mpopt.verbose
    if ~converged
        fprintf('\nNewton''s method corrector did not converge in %d iterations.\n', i);
    end
end
