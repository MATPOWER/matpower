function [V, converged, i, Sbus] = zgausspf(Ybus, Sbus, V0, ref, pv, pq, Bpp, mpopt, stage, VPQ, SPQ)
%ZGAUSSPF  Solves the power flow using an Implicit Z-bus Gauss method.
%   [V, CONVERGED, I] = ZGAUSSPF(YBUS, SBUS, V0, REF, PV, PQ, BPP, MPOPT)
%   solves for bus voltages given the full system admittance matrix (for
%   all buses), the complex bus power injection vector (all buses),
%   the initial vector of complex bus voltages, column vectors with the
%   lists of bus indices for the swing bus, PV buses, and PQ buses,
%   respectively, and the fast-decoupled B double-prime matrix (all buses)
%   for Q updates at PV buses. The bus voltage vector contains the set point
%   for generator (including ref bus) buses, and the reference angle of the
%   swing bus, as well as an initial guess for remaining magnitudes and
%   angles. MPOPT is a MATPOWER options struct which can be used to
%   set the termination tolerance, maximum number of iterations, and
%   output options (see MPOPTION for details). Uses default options
%   if this parameter is not given. Returns the final complex voltages,
%   a flag which indicates whether it converged or not, and the number
%   of iterations performed.
%
%   NOTE: This method does not scale well with the number of generators
%       and seems to have serious problems with some systems with many
%       PV buses.
%
%   See also RUNPF.

%   MATPOWER
%   Copyright (c) 1996-2019 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% default arguments
if nargin < 9
    stage = 1;
    if nargin < 8
        mpopt = mpoption;
    end
end

%% options
tol     = mpopt.pf.tol;
max_it  = mpopt.pf.zg.max_it;
complex = 0;    %% use 1 = complex formulation, 0 = use real formulation
switch stage
    case 1
        pv_method = 2;  %% 0 = simple voltage magnitude reset
                        %% 1 = sens real(V) to imag(I)
                        %% 2 = dVm/dQ (fast-decoupled Jac)
                        %% 3 = voltage correction PF by Rajicic, Ackovski, Taleski
                        %% 4 = homotopy method
    case 2
        pv_method = 0;
    case 3
        pv_method = 4;
end
if have_fcn('matlab') && have_fcn('matlab', 'vnum') < 7.3
    lu_vec = 0;     %% lu(..., 'vector') syntax not supported
else
    lu_vec = 1;
end
if mpopt.verbose
    str = sprintf('Implicit Z-bus Gauss power flow (stage %d)', stage);
end

%% initialize
converged_htpy = 0;
i = 0;
npv = length(pv);
npq = length(pq);
max_dV = 0;
lambda = 0;
if stage == 3
    V = VPQ;
    dl = 0.01;
    lambda = -dl;           %-----  TESTING  -----

    %% compute PV current injections from stage 2
    Spq = real(Sbus(pv)) + 1j * imag(SPQ(pv));
    I2PQ = conj(Spq ./ VPQ(pv));
else
    V = V0;
    dl = 0.9;       %% terminate after single outer iteration
end
if stage == 1
    Sbus0 = Sbus;
end

%% shift voltage angles if necessary
% Varef = angle(V(ref));
% if Varef ~= 0
%     V = V * exp(-1j * Varef);
% end
if complex
    V1 = V(ref);
else
    V1 = [real(V(ref)); imag(V(ref))];
end

%% create sub-matrices, separating ref bus from other buses
if complex
%     Y11 = Ybus(ref,ref);
%     Y12 = Ybus(ref, [pv;pq]);
    Y21 = Ybus([pv;pq], ref);
    Y21_V1 = Y21 * V1;
    Y22 = Ybus([pv;pq], [pv;pq]);
else
    G = real(Ybus);
    B = imag(Ybus);
%     Y11 = [ G(ref,ref) -B(ref,ref);
%             B(ref,ref)  G(ref,ref) ];
%     Y12 = [ G(ref,pv) -B(ref,pv) G(ref,pq) -B(ref,pq);
%             B(ref,pv)  G(ref,pv) B(ref,pq)  G(ref,pq) ];
    Y21 = [ G(pv,ref) -B(pv,ref);
            B(pv,ref)  G(pv,ref);
            G(pq,ref) -B(pq,ref);
            B(pq,ref)  G(pq,ref) ];
    Y21_V1 = Y21 * V1;
    Y22 = [ G(pv,pv) -B(pv,pv) G(pv,pq) -B(pv,pq);
            B(pv,pv)  G(pv,pv) B(pv,pq)  G(pv,pq);
            G(pq,pv) -B(pq,pv) G(pq,pq) -B(pq,pq);
            B(pq,pv)  G(pq,pv) B(pq,pq)  G(pq,pq)  ];
end
if lu_vec
    [L, U, p, q] = lu(Y22, 'vector');
    [junk, iq] = sort(q);
else
    [L, U, P] = lu(Y22);
end

%% initialize PV bus handling
if npv
    Vmpv0 = abs(V0(pv));    %% voltage setpoints for PV buses

    switch pv_method
        case 1
            %% Essentially, the following, only more efficiently ...
            % Y22_inv = inv(Y22);
            % dVdI = -imag(Y22_inv(1:npv, 1:npv));      % complex
            % dVdI = Y22_inv(1:npv, npv+(1:npv));       % real
            if complex
                rhs = sparse(1:npv, 1:npv, 1j, npv+npq, npv);
            else
                rhs = sparse(npv+(1:npv), 1:npv, 1, 2*(npv+npq), npv);
            end
            % cols_of_Y22_inv = Y22 \ rhs;
            if lu_vec
                cols_of_Y22_inv = U \ (L \ rhs(p, :));
                cols_of_Y22_inv = cols_of_Y22_inv(iq, :);
            else
                cols_of_Y22_inv = U \ (L \ P * rhs);
            end
            if complex
                dVdI = real(cols_of_Y22_inv(1:npv, :));
            else
                dVdI = cols_of_Y22_inv(1:npv, :);
            end

            [LL, UU, PP] = lu(dVdI);    %% not sparse, so don't use 'vector' version
        case 2
            if lu_vec
                [LBpp, UBpp, pBpp, qBpp] = lu(Bpp(pq, pq), 'vector');
                [junk, iqBpp] = sort(qBpp);
            else
                [LBpp, UBpp, PBpp] = lu(Bpp(pq, pq));
            end
        case 3
            %% compute "breakpoint" impedance matrix Zpv
            npvq = npv + npq;
            if complex
                Z = zeros(npvq, npv);
                II = zeros(npvq, 1);
            else
                Z = zeros(2*npvq, npv);
                II = zeros(2*npvq, 1);
            end
            for k = 1:npv
                II(k) = 1;
                Z(:, k) = Y22 \ II;
                II(k) = 0;
            end
            if complex
                Zpv = Z(1:npv, :);
            else
                Zpv = Z(1:2*npv, :);
            end
    end
end

%% compute initial power mismatch
if npv      %% compute Q injection at PV buses for current V
    Qpv = imag( V(pv) .* conj(Ybus(pv, :) * V) );
    Sbus(pv) = Sbus(pv) + 1j * (Qpv - imag(Sbus(pv)));
end
misS = Sbus - V .* conj(Ybus * V);
normS = norm(misS([pv; pq]), Inf);

%% check tolerance
if mpopt.verbose > 1
    if stage == 3
        fprintf('\n it   it2     ∆V (p.u.)    max abs(S) mis (p.u.)    max V mis (PV) (p.u.)\n');
        fprintf('---- ----    -----------  -----------------------  -----------------------\n');
        fprintf('  0    0                           %10.3e\n', normS);
    else
        fprintf('\n it      ∆V (p.u.)    max abs(S) mismatch (p.u.)    max V mismatch (PV) (p.u.)\n');
        fprintf('----    -----------  ----------------------------  ----------------------------\n');
        fprintf('  0                            %10.3e\n', normS);
    end
end

%% do homotopy iterations
ii = 0;
while (~converged_htpy && lambda + dl < 1)      %% outer loop
    %% update iteration counter & lambda
    ii = ii + 1;
    lambda = lambda + dl;

    i = 0;
    converged = 0;

%% do implicit Zbus Gauss iterations
while (~converged && i < max_it)
    %% update iteration counter
    i = i + 1;

    %% save voltage from previous iteration
    Vp = V;

    if npv && (i > 1 || pv_method ~= 3) %% update Q injections @ PV buses based on V mismatch
        %% for compatibility with voltage correction PF by Rajicic, Ackovski, Taleski
        %% we skip updating Q injections @ PV buses on first iteration

        %% compute voltage mismatch at PV buses
        Vmpv = abs(V(pv));
        dV = Vmpv0 - Vmpv;
        [max_dV, k] = max(abs(dV));

        % Three alternate approaches (that seem to work about equally well):
        %   1 - use precomputed sensitivity of real(V) to imag(I)
        %   2 - use sensitivity of Vm to Q evaluated at current V (fast-decoupled)
        %   3 - use method equivalent to Voltage Correction PF by Rajicic, Ackovski, Taleski
        switch pv_method
            case {1, 2}
                %% compute Q injection at current V
                %% (sometimes improves convergence for pv_method=1,2)
                Qpv = imag( V(pv) .* conj(Ybus(pv, :) * V) );
                Sbus(pv) = Sbus(pv) + 1j * (Qpv - imag(Sbus(pv)));
        end
        switch pv_method
            case 1      %% precomputed sensitivity of real(V) to imag(I)
                %% estimate corresponding change in imag(I) injection
                dQ = -UU \  (LL \ (PP * dV));    %% dQ = -dI = -dVdI \ dV;
            case 2      %% sensitivity of Vm to Q (fast-decoupled Jacobian)
                % dVpq = Bpp(pq, pq) \ (-Bpp(pq, pv) * dV);
                if lu_vec
                    dVpq = UBpp \  (LBpp \ (-Bpp(pq(pBpp), pv) * dV));
                    dVpq = dVpq(iqBpp);
                else
                    dVpq = UBpp \  (LBpp \ (PBpp * (-Bpp(pq, pv) * dV)));
                end
                dQ = Bpp(pv, pq) * dVpq + Bpp(pv, pv) * dV;
            case 3
                dE = (Vmpv0 ./ Vmpv - 1) .* real(V(pv));
                if complex
                    dD = imag(Zpv) \ dE;
                else
                    dD = Zpv(npv+(1:npv), :) \ dE;
                end
                if mpopt.pf.radial.vcorr    %% do voltage correction step?
                    dC = dD .* imag(V(pv)) ./ real(V(pv));
                    if complex
                        dI = -[dC + 1j * dD; zeros(npq, 1)];
                    else
                        dI = -[dC; dD; zeros(2*npq, 1)];
                    end
                    dVV = Y22 \ dI;
                    if complex
                        V([pv; pq]) = V([pv; pq]) + dVV;
                    else
                        V(pv) = V(pv) + (dVV(1:npv) + 1j * dVV(npv+(1:npv)));
                        V(pq) = V(pq) + (dVV(2*npv+(1:npq)) + 1j * dVV(2*npv+npq+(1:npq)));
                    end
                    Vmpv = abs(V(pv));
                end
                dQ = dD .* Vmpv.^2 ./ real(V(pv));
        end
        switch pv_method
            case {1, 2, 3}
                %% update Sbus
                Sbus(pv) = Sbus(pv) + 1j * dQ;

                %% set voltage magnitude at PV buses
                %% doesn't consistently improve convergence for pv_method=1,2
%                V(pv) = V(pv) .* Vmpv0 ./ abs(V(pv));
            case 0
                %% compute Q injection at current V
                %% (updating Q before V converges more consistently)
                Qpv = imag( V(pv) .* conj(Ybus(pv, :) * V) );
                Sbus(pv) = Sbus(pv) + 1j * (Qpv - imag(Sbus(pv)));

                %% set voltage magnitude at PV buses
                V(pv) = V(pv) .* Vmpv0 ./ abs(V(pv));

%                 %% compute Q injection at current V
%                 Qpv = imag( V(pv) .* conj(Ybus(pv, :) * V) );
%                 Sbus(pv) = Sbus(pv) + 1j * (Qpv - imag(Sbus(pv)));
            case 4
                %% compute Q injection at current V
                Qpv = imag( V(pv) .* conj(Ybus(pv, :) * V) );
                Sbus(pv) = Sbus(pv) + 1j * (Qpv - imag(Sbus(pv)));

                %% set voltage magnitude at PV buses
%                 Vs = V;
%                 Vs(pv) = Vs(pv) .* Vmpv0 ./ abs(Vs(pv));
        end
    end

    %% update currents
    if stage == 3
        Vspec = V(pv) .* Vmpv0 ./ abs(V(pv));
        Ispec = conj(Sbus(pv) ./ Vspec);
%         Ispec = conj(Sbus(pv) ./ Vs(pv));
%         I2PQ = conj(SPQ(pv) ./ V(pv));
        Ipv = lambda * Ispec + (1-lambda) * I2PQ;
    else
        Ipv = conj(Sbus(pv) ./ V(pv));
    end
    Ipq = conj(Sbus(pq) ./ V(pq));
    if complex
        I2 = [Ipv; Ipq];
    else
        I2 = [real(Ipv); imag(Ipv); real(Ipq); imag(Ipq)];
    end

    %% solve for new voltages (except at ref & PV buses)
    if lu_vec
        V2 = U \  (L \ (I2(p) - Y21_V1(p)));
        V2 = V2(iq);
    else
        V2 = U \  (L \ (P * (I2 - Y21_V1)));
    end
    % V2 = Y22 \ (I2 - Y21_V1);
    if complex
        V(pv) = V2(1:npv);
        V(pq) = V2(npv+1:npv+npq);
    else
        V(pv) = V2(1:npv) + 1j * V2(npv+1:2*npv);
        V(pq) = V2(2*npv+1:2*npv+npq) + 1j * V2(2*npv+npq+1:2*npv+2*npq);
    end

    %% check for convergence
    normV = norm(V-Vp, Inf);
    misS = Sbus - V .* conj(Ybus * V);
    normS = norm(misS([pv; pq]), Inf);
% if stage ~= 2
%     P1 = real(Sbus);
%     P2 = real(V .* conj(Ybus * V));
%     Q1 = imag(Sbus);
%     Q2 = imag(V .* conj(Ybus * V));
% %     normS2 = norm(Sbus([pv; pq]) - SPQ([pv; pq]), Inf)
% %     P = [P1(pv) P2(pv)]
% %     Q = [Q1(pv) Q2(pv)]
% %     Sbus([pv; pq]) - SPQ([pv; pq])
% %     [imag(Sbus(pv)) imag(SPQ(pv))]
% end

    %% print progress
    if mpopt.verbose > 1
        if stage == 3
            fprintf('%3d %4d     %10.3e            %10.3e             %10.3e\n', ii, i, normV, normS, max_dV);
        else
            fprintf('%3d     %10.3e             %10.3e                  %10.3e\n', i, normV, normS, max_dV);
        end
    end

    %% check Z-Gauss convergence
    ttol = tol; if stage == 3, ttol = tol / 10; end %-----  TESTING  -----
    if normV < ttol                                 %-----  TESTING  -----
%     if normV < tol
        converged = 1;
    elseif normV > 10   %% diverging, time to bail out
        converged = -1;
    elseif stage == 3
        %% homotopy voltage update
        Vspec = V(pv) .* Vmpv0 ./ abs(V(pv));
        V(pv) = lambda * Vspec + (1-lambda) * V(pv);
%         V(pv) = lambda * Vs(pv) + (1-lambda) * V(pv);
    end
end     %% end Z-Gauss iterations

    %% check for homotopy convergence
    if stage == 3
        normV = norm(abs(V(pv)) - Vmpv0);
        if normV < tol
            converged_htpy = 1;
        elseif normV > 10   %% diverging, time to bail out
            converged_htpy = -1;
        end
    end
end

%% shift voltage angles back if necessary
% if Varef ~= 0
%     V = V * exp(1j * Varef);
% end

if mpopt.verbose
    switch converged
        case 1
            fprintf('%s converged in %d iterations.\n', str, i);
        case 0
            fprintf('%s did not converge in %d iterations.\n', str, i);
        case -1
            fprintf('%s diverged in %d iterations.\n', str, i);
    end
end

%% use homotopy method
if stage == 1 && converged ~= 1
    %% call stage 2
    [VPQ, converged, iterations, SPQ] = zgausspf(Ybus, Sbus0, V0, ref, [], sort([pq; pv]), Bpp, mpopt, 2);

    %% call homotopy
    [V, converged, i] = zgausspf(Ybus, Sbus0, V0, ref, pv, pq, Bpp, mpopt, 3, VPQ, SPQ);
end

if converged == -1
    converged = 0;
end
