function [V, converged, i] = zgausspf(Ybus, Sbus, V0, ref, pv, pq, mpopt)
%ZGAUSSPF  Solves the power flow using an Implicit Z-bus Gauss method.
%   [V, CONVERGED, I] = ZGAUSSPF(YBUS, SBUS, V0, REF, PV, PQ, MPOPT)
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
%   NOTE: This method does not scale well with the number of generators
%       and seems to have serious problems with some systems with multiple
%       PV buses.
%
%   See also RUNPF.

%   MATPOWER
%   Copyright (c) 1996-2016 by Power System Engineering Research Center (PSERC)
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
tol     = mpopt.pf.tol;
max_it  = mpopt.pf.zg.max_it;
if have_fcn('matlab') && have_fcn('matlab', 'vnum') < 7.3
    lu_vec = 0;     %% lu(..., 'vector') syntax not supported
else
    lu_vec = 1;
end

%% initialize
converged = 0;
i = 0;
j = sqrt(-1);
V = V0;
Vm = abs(V);
npv = length(pv);
npq = length(pq);
nb = length(V);
max_dV = 0;

%% shift voltage angles if necessary
% Varef = angle(V(ref));
% if Varef ~= 0
%     V = V * exp(-j * Varef);
% end
V1 = [real(V(ref)); imag(V(ref))];

%% create sub-matrices
G = real(Ybus);
B = imag(Ybus);
% Y11 = [ G(ref,ref) -B(ref,ref);
%         B(ref,ref)  G(ref,ref) ];
% Y12 = [ G(ref,pv) -B(ref,pv) G(ref,pq) -B(ref,pq);
%         B(ref,pv)  G(ref,pv) B(ref,pq)  G(ref,pq) ];
Y21 = [ G(pv,ref) -B(pv,ref);
        B(pv,ref)  G(pv,ref);
        G(pq,ref) -B(pq,ref);
        B(pq,ref)  G(pq,ref) ];
Y21_V1 = Y21 * V1;
Y22 = [ G(pv,pv) -B(pv,pv) G(pv,pq) -B(pv,pq);
        B(pv,pv)  G(pv,pv) B(pv,pq)  G(pv,pq);
        G(pq,pv) -B(pq,pv) G(pq,pq) -B(pq,pq);
        B(pq,pv)  G(pq,pv) B(pq,pq)  G(pq,pq)  ];
if lu_vec
    [L, U, p, q] = lu(Y22, 'vector');
    [junk, iq] = sort(q);
else
    [L, U, P] = lu(Y22);
end

if npv
    Vmpv0 = Vm(pv);     %% voltage setpoints for PV buses

    Ipv = conj(Sbus(pv) ./ V(pv));
    Ipq = conj(Sbus(pq) ./ V(pq));
    I2 = [real(Ipv); imag(Ipv); real(Ipq); imag(Ipq)];

    %% compute sensitivities for PV buses
    dVdI = zeros(npv, npv);
    ptb = 0.01;
    for k = 1:npv
        if lu_vec
            V2b = U \  (L \ (I2(p) - Y21_V1(p)));
            V2b = V2b(iq);
        else
            V2b = U \  (L \ (P * (I2 - Y21_V1)));
        end
        % V2b = Y22 \ (I2 - Y21_V1);
        I2a = I2;
        I2a(npv+k) = I2a(npv+k) + ptb;
        if lu_vec
            V2a = U \  (L \ (I2a(p) - Y21_V1(p)));
            V2a = V2a(iq);
        else
            V2a = U \  (L \ (P * (I2 - Y21_V1)));
        end
        % V2a = Y22 \ (I2a - Y21_V1);
        dVdI(:,k) = (V2a(1:npv) - V2b(1:npv)) / ptb;
    end
    [LL, UU, PP] = lu(dVdI);    %% not sparse, so don't use 'vector' version
end

%% check tolerance
if mpopt.verbose > 1
    fprintf('\n it      ∆V (p.u.)    max abs(S) mismatch (p.u.)    max V mismatch (PV) (p.u.) ');
    fprintf('\n----    -----------  ----------------------------  ----------------------------');
end

%% do implicit Zbus Gauss iterations
while (~converged && i < max_it)
    %% update iteration counter
    i = i + 1;

    %% save voltage from previous iteration
    Vp = V;

    if npv
% (using full Jacobian does not seem to improve convergence)
%         dVdQ = zeros(npv, npv);
%         for k = 1:npv
%             %% evaluate Jacobian
%             [dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, V);
%     
%             j11 = real(dSbus_dVa([pv; pq], [pv; pq]));
%             j12 = real(dSbus_dVm([pv; pq], [pv; pq]));
%             j21 = imag(dSbus_dVa([pv; pq], [pv; pq]));
%             j22 = imag(dSbus_dVm([pv; pq], [pv; pq]));
%     
%             J = [   j11 j12;
%                     j21 j22;    ];
%             
%             Sb = zeros(2*npv+2*npq, 1);
%             Sb(npv+npq+k) = ptb;
%             dV = J \ Sb;
%             
%             dVdQ(:,k) = dV(npv+npq+(1:npv)) / ptb;
%         end

        %% update Q injections @ PV buses based on V mismatch
        %% compute Q injection at current V (seems to help convergence)
        Qpv = imag( V(pv) .* conj(Ybus(pv, :) * V) );
        Sbus(pv) = Sbus(pv) + j * (Qpv - imag(Sbus(pv)));

        %% compute voltage mismatch at PV buses
        Vmpv = abs(V(pv));
        dV = Vmpv0 - Vmpv;
        [max_dV, k] = max(abs(dV));
%        dV([1:k-1 k+1:end]) = 0;   %% one at a time?

        %% estimate corresponding change in Q injection
        dI = UU \  (LL \ (PP * dV));
%         dI = dVdI \ dV;
%         dQ = dVdQ \ dV;

        %% update Sbus
        Sbus(pv) = Sbus(pv) - j * dI;
%         Sbus(pv) = Sbus(pv) + j * dQ;
    end

    %% set voltage magnitude at PV buses
% (this line does not seem to consistently improve convergence)
%     V(pv) = V(pv) ./ abs(V(pv)) .* abs(V0(pv));

    %% update currents
    Ipv = conj(Sbus(pv) ./ V(pv));
    Ipq = conj(Sbus(pq) ./ V(pq));
    I2 = [real(Ipv); imag(Ipv); real(Ipq); imag(Ipq)];
    
    %% solve for new voltages (except at ref & PV buses)
    if lu_vec
        V2 = U \  (L \ (I2(p) - Y21_V1(p)));
        V2 = V2(iq);
    else
        V2 = U \  (L \ (P * (I2 - Y21_V1)));
    end
    % V2 = Y22 \ (I2 - Y21_V1);
    V(pv) = V2(1:npv) + j * V2(npv+1:2*npv);
    V(pq) = V2(2*npv+1:2*npv+npq) + j * V2(2*npv+npq+1:2*npv+2*npq);

    %% check for convergence
    normV = norm(V-Vp, Inf);
    normS = norm(Sbus([pv; pq]) - V([pv; pq]) .* conj([Ipv; Ipq]), Inf);
    if mpopt.verbose > 1
%        fprintf('\n%3d        %10.3e %10.6g %10.6g %10.6g %10.6g', i, normV, V(pv(1)), V(pv(2)), imag(Sbus(pv(1))), imag(Sbus(pv(2))));
%        fprintf('\n%3d        %10.3e %10.6g %10.6g', i, normV, V(pv(1)), imag(Sbus(pv(1))));
        fprintf('\n%3d     %10.3e             %10.3e                  %10.3e', i, normV, normS, max_dV);
    end
    if normV < tol
        converged = 1;
        if mpopt.verbose
            fprintf('\nImplicit Z-bus Gauss power flow converged in %d iterations.\n', i);
        end
    end
end

%% shift voltage angles back if necessary
% if Varef ~= 0
%     V = V * exp(j * Varef);
% end

if mpopt.verbose
    if ~converged
        fprintf('\nImplicit Z-bus Gauss power flow did not converge in %d iterations.\n', i);
    end
end
