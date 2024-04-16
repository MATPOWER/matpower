classdef math_model_pf_acps < mp.math_model_pf_ac & mp.mm_shared_pfcpf_acps
% mp.math_model_pf_acps - Power flow (PF) **math model** for AC-polar-power formulation.
%
% Implements formulation-specific node balance constraints and inherits
% from formulation-specific class for shared PF/CPF code.
%
% Also includes implementations of methods specific to fast-decoupled
% power flow.

%   MATPOWER
%   Copyright (c) 2021-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        function tag = form_tag(obj)
            %

            tag = 'acps';
        end

        function name = form_name(obj)
            %

            name = 'AC-polar-power';
        end

        function obj = add_node_balance_constraints(obj, nm, dm, mpopt)
            %

            alg = mpopt.pf.alg;
            ad = obj.aux_data;

            %% power balance constraints
            switch alg
                case  {'FDXB', 'FDBX'}
                    fcn = @(x)node_balance_equations(obj, x, nm, 1);
                otherwise
                    fcn = @(x)node_balance_equations(obj, x, nm);
            end
            obj.add_nln_constraint({'Pmis', 'Qmis'}, [ad.npv+ad.npq;ad.npq], 1, fcn, []);
        end

        function x = gs_x_update(obj, x, f, nm, dm, mpopt)
            %

            alg = mpopt.pf.alg;
            ad = obj.aux_data;

            %% update network model state ([v_; z_]) from math model state (x)
            [v_, z_] = obj.convert_x_m2n(x, nm, 1);

            [pv, pq, npv, npq, Y] = deal(ad.pv, ad.pq, ad.npv, ad.npq, ad.Y);

            %% total nodal complex bus power extractions
            SS = zeros(size(v_));
            SS(pv) = f(1:npv);
            SS(pq) = f(npv+1:npv+npq) + 1j * f(npv+npq+1:npv+2*npq);
%            SS = C * nm.port_inj_power([v_; z_], 1);

            %% complex net nodal injection (from all but constant Z elements)
            S0 = v_ .* conj(Y * v_) - SS;

            %% update voltage
            %% at PQ buses
            for k = pq'
                v_(k) = v_(k) + (conj(S0(k)/v_(k)) - Y(k,:) * v_) / Y(k, k);
            end

            %% at PV buses
            if npv
                for k = pv'
                    S0(k) = real(S0(k)) + 1j * imag( v_(k) * conj(Y(k,:) * v_) );
                    v_(k) = v_(k) + (conj(S0(k)/v_(k)) - Y(k,:) * v_) / Y(k, k);
                end
            end

            x = [angle(v_([pv; pq])); abs(v_(pq))];
        end

        function x = zg_x_update(obj, x, f, nm, dm, mpopt)
            %

            alg = mpopt.pf.alg;
            ad = obj.aux_data;

            %% update network model state ([v_; z_]) from math model state (x)
            [v_, z_] = obj.convert_x_m2n(x, nm, 1);

            [pv, pq, ref, npv, npq] = deal(ad.pv, ad.pq, ad.ref, ad.npv, ad.npq);
            pvq = [pv; pq];

            %% build and cache S0 and Y21_v1
            if isfield(ad, 'S0')
                [Y21_v1, S0] = deal(ad.Y21_v1, ad.S0);
            else
                Y21_v1 = ad.Y21 * v_(ref);
                SS = zeros(size(v_));
                SS(pv) = f(1:npv);
                SS(pq) = f(npv+1:npv+npq) + 1j * f(npv+npq+1:npv+2*npq);
    %            SS = C * nm.port_inj_power([v_; z_], 1);

                %% complex net nodal injection (from all but constant Z elements)
                S0 = v_ .* conj(ad.Y * v_) - SS;
                S0(ref) = 0;

                %% cache 'em
                [ad.Y21_v1, ad.S0] = deal(Y21_v1, S0);
                obj.aux_data = ad;
            end

            if npv  %% update Q injections at PV buses based on vm mismatch
                %% compute and cache initial vm at PV buses and
                %% factored fast-decoupled Bpp matrix
                if isfield(ad, 'vmpv0')
                    [vmpv, vmpv0, Bpp, LBpp, UBpp, pBpp, iqBpp] = ...
                        deal(ad.vmpv, ad.vmpv0, ad.Bpp, ad.LBpp, ad.UBpp, ad.pBpp, ad.iqBpp);
                else
                    vmpv0 = abs(v_(pv));
                    vmpv = vmpv0;

                    %% modify data model to form Bpp (B double prime)
                    dm2 = obj.fdpf_B_matrix_models(dm, 'FDBX');

                    %% build network models and get admittance matrices
                    nm2 = feval(class(nm));
                    nm2.build(dm2);
                    [Y2, L, M] = nm2.get_params([], {'Y', 'L', 'M'});
                    if any(any(L)) || any(any(M))
                        error('mp.math_model_pf_acps.zg_x_update: B matrix for Z-bus Gauss w/PV buses not implemented for models with non-zero L and/or M matrices.')
                    end
                    Bpp = -nm2.C * imag(Y2) * nm2.C';

                    [LBpp, UBpp, pBpp, qBpp] = lu(Bpp(pq, pq), 'vector');
                    [junk, iqBpp] = sort(qBpp);

                    %% cache 'em
                    [ad.vmpv, ad.vmpv0, ad.Bpp, ad.LBpp, ad.UBpp, ad.pBpp, ad.iqBpp] = ...
                        deal(vmpv, vmpv0, Bpp, LBpp, UBpp, pBpp, iqBpp);
                    obj.aux_data = ad;
                end

                %% compute voltage mismatches at PV buses
                v_(pv) = vmpv .* v_(pv) ./ abs(v_(pv));
                dV = vmpv0 - vmpv;
%                 [max_dV, k] = max(abs(dV));
%                 fprintf('       %10.3e', max_dV)

                %% compute Q injection at current V
                %% (sometimes improves convergence)
                Qpv = imag( v_(pv) .* conj(ad.Y(pv, :) * v_) );
                S0(pv) = S0(pv) + 1j * (Qpv - imag(S0(pv)));

                % dVpq = Bpp(pq, pq) \ (-Bpp(pq, pv) * dV);
                dVpq = UBpp \  (LBpp \ (-Bpp(pq(pBpp), pv) * dV));
                dVpq = dVpq(iqBpp);
                dQ = Bpp(pv, pq) * dVpq + Bpp(pv, pv) * dV;

                %% update S0
                S0(pv) = S0(pv) + 1j * dQ;
            end

            %% complex current injections
            I2 = conj(S0(pvq) ./ v_(pvq));

            V2 = ad.U \  (ad.L \ (I2(ad.p) - Y21_v1(ad.p)));
            V2 = V2(ad.iq);

            v_(pv) = V2(1:npv);
            v_(pq) = V2(npv+1:npv+npq);
            obj.aux_data.vmpv = abs(v_(pv));

            x(1:ad.npv+ad.npq) = angle(v_(pvq));
            x(ad.npv+ad.npq+1:ad.npv+2*ad.npq) = abs(v_(pq));
        end

        function JJ = fd_jac_approx(obj, nm, dm, mpopt)
            %

            alg = mpopt.pf.alg;

            %% create copies of data model for building B prime, B double prime
            [dm1, dm2] = obj.fdpf_B_matrix_models(dm, alg);

            %% build network models and get admittance matrices
            nm1 = feval(class(nm));
            nm2 = feval(class(nm));
            nm1.build(dm1);
            nm2.build(dm2);
            [Y1, L, M] = nm1.get_params([], {'Y', 'L', 'M'});
            Y2 = nm2.get_params();
            if any(any(L)) || any(any(M))
                error('mp.math_model_pf_acps.fd_jac_approx: fast-decoupled Jacobian approximation not implemented for models with non-zero L and/or M matrices.')
            end

            %% form reduced Bp and Bpp matrices
            ad = obj.aux_data;
            Cp  = nm1.C([ad.pv; ad.pq], :);
            Cpp = nm2.C(ad.pq, :);
            Bp  = -imag( Cp  * Y1 * Cp' );
            Bpp = -imag( Cpp * Y2 * Cpp' );
            JJ = {Bp, Bpp};
        end

        function [dm1, dm2] = fdpf_B_matrix_models(obj, dm, alg)
            %

            %% [dmp, dmpp] = obj.fdpf_B_matrix_models(dm, alg)
            %% dmpp = obj.fdpf_B_matrix_models(dm, alg)
            %% returns copies of dm used for building B prime, B double prime
            %% for fast-decoupled power flow

            %% modify data model to form Bp (B prime)
            if nargout > 1      %% for both Bp and Bpp
                dm1 = dm.copy();
                if dm1.elements.has_name('shunt')
                    dm1.elements.shunt.tab.bs(:) = 0;   %% zero out shunts at buses
                end
                dm2 = dm1.copy();
                dm1.elements.branch.tab.b_fr(:) = 0;    %% zero out line charging shunts
                dm1.elements.branch.tab.b_to(:) = 0;
                dm1.elements.branch.tab.tm(:) = 1;      %% cancel out taps
                if strcmp(alg, 'FDXB')                  %% if XB method
                    dm1.elements.branch.tab.r(:) = 0;   %% zero out line resistance
                end
                dm1 = dm1.build_params();
            else
                dm2 = dm.copy();
            end

            %% modify data model to form Bpp (B double prime)
            dm2.elements.branch.tab.ta(:) = 0;      %% zero out phase shifters
            if strcmp(alg, 'FDBX')                  %% if BX method
                dm2.elements.branch.tab.r(:) = 0;   %% zero out line resistance
            end

            if nargout > 1      %% for both Bp and Bpp
                dm2 = dm2.build_params();
            else                %% for just Bpp
                dm1 = dm2.build_params();
            end
        end
    end     %% methods
end         %% classdef
