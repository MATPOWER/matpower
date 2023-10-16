classdef form_acp < mp.form_ac
% mp.form_acp - Base class for |MATPOWER| AC polar **formulations**.
%
% Used as a mix-in class for all **network model element** classes
% with an AC network model formulation with a **polar** repesentation for
% voltages. That is, each concrete network model element class with an AC
% polar formulation must inherit, at least indirectly, from
% both mp.nm_element and mp.form_acp.
%
% Provides implementation of evaluation of voltage-related Jacobian and
% Hessian terms needed by some mp.form_ac methods.
%
% mp.form_dc Methods:
%   * form_name - get char array w/name of formulation (``'AC-polar'``)
%   * form_tag - get char array w/short label of formulation (``'acp'``)
%   * model_vvars - get cell array of names of voltage state variables (``{'va', 'vm'}``)
%   * port_inj_current_jac - compute voltage-related terms of current injection Jacobian
%   * port_inj_current_hess_v - compute voltage-related terms of current injection Hessian
%   * port_inj_current_hess_vz - compute voltage/non-voltage-related terms of current injection Hessian
%   * port_inj_power_jac - compute voltage-related terms of power injection Jacobian
%   * port_inj_power_hess_v - compute voltage-related terms of power injection Hessian
%   * port_inj_power_hess_vz - compute voltage/non-voltage-related terms of power injection Hessian
%   * aux_data_va_vm - return voltage angles/magnitudes from auxiliary data
%
% For more details, see the :ref:`sec_nm_formulations_ac` section in the
% |MATPOWER-Dev-Manual| and the derivations in |TN5|.
%
% See also mp.form, mp.form_ac, mp.form_acc, mp.nm_element.

%   MATPOWER
%   Copyright (c) 2019-2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        function name = form_name(obj)
            % Get user-readable name of formulation, i.e. ``'AC-polar'``.
            %
            % See :meth:`mp.form.form_name`.

            name = 'AC-polar';
        end

        function tag = form_tag(obj)
            % Get short label of formulation, i.e. ``'acp'``.
            %
            % See :meth:`mp.form.form_tag`.

            tag = 'acp';
        end

        function vtypes = model_vvars(obj)
            % Get cell array of names of voltage state variables, i.e. ``{'va', 'vm'}``.
            %
            % See :meth:`mp.form.model_vvars`.

            vtypes = {'va', 'vm'};
        end

        function [Iva, Ivm] = port_inj_current_jac(obj, n, v_, Y, M, invdiagvic, diagSlincJ)
            % Compute voltage-related terms of current injection Jacobian.
            % ::
            %
            %   [Iva, Ivm] = nme.port_inj_current_jac(n, v_, Y, M, invdiagvic, diagSlincJ)
            %
            % Called by mp.form_ac.port_inj_current to compute
            % voltage-related Jacobian terms.

            %% intermediate terms
            diagv = sparse(1:n, 1:n, v_, n, n);
            C = invdiagvic * (diagSlincJ - conj(M * diagv));
            D = sparse(1:n, 1:n, 1 ./ abs(v_), n, n);

%             %% linear current term
%             Ivm = Y * diagv;
%             Iva = 1j * Ivm;
%             Ivm = Ivm * D;
%
%             %% + current from linear power term
%             Iva = Iva + 1j * C;
%             Ivm = Ivm - C * D;

            A = Y * diagv;
            Iva = 1j * (A + C);
            Ivm = (A - C) * D;
        end

        function [Ivava, Ivavm, Ivmvm] = port_inj_current_hess_v(obj, x_, lam, v_, z_, diaginvic, Y, M, diagSlincJ, dlamJ)
            % Compute voltage-related terms of current injection Hessian.
            % ::
            %
            %   [Ivava, Ivavm, Ivmvm] = nme.port_inj_current_hess_v(x_, lam)
            %   [Ivava, Ivavm, Ivmvm] = nme.port_inj_current_hess_v(x_, lam, sysx)
            %   [Ivava, Ivavm, Ivmvm] = nme.port_inj_current_hess_v(x_, lam, sysx, idx)
            %   [...] = nme.port_inj_current_hess_vz(x_, lam, v_, z_, diaginvic, Y, M, diagSlincJ, dlamJ)
            %
            % Called by mp.form_ac.port_inj_current_hess to compute
            % voltage-related Hessian terms.

            if nargin < 10
                sysx = v_;
                if nargin < 5
                    idx = []
                else
                    idx = z_;
                end
                [Y, M, N, s] = obj.get_params(idx, {'Y', 'M', 'N', 's'});
                [v_, z_, vi_] = obj.x2vz(x_, sysx, idx);

                %% compute linear power injections
                if isempty(z_)
                    Slin = M*v_ + s;
                else
                    Slin = M*v_ + N*z_ + s;
                end

                n  = length(v_);    %% number of all port voltages
                ni = length(vi_);   %% number of selected port voltages
                diaginvic = sparse(1:ni, 1:ni, 1 ./ conj(vi_), ni, ni);
                if isempty(idx)     %% all ports
                    diagSlincJ = sparse(1:n, 1:n, conj(Slin), n, n);
                    dlamJ      = sparse(1:n, 1:n, lam, n, n);
                else                %% selected ports
                    diagSlincJ = sparse(1:ni, idx, conj(Slin), ni, n);
                    dlamJ      = sparse(1:ni, idx, lam, ni, n);
                end
            else
                n  = length(v_);    %% number of all port voltages
                ni = length(lam);
            end

            %% intermediate terms
            diagv  = sparse(1:n, 1:n, v_, n, n);
            A = diaginvic * diagSlincJ;
            B = diaginvic * conj(M);
            % B2 = B * conj(diagv);
            % C = A - B2;
            D = sparse(1:n, 1:n, 1 ./ abs(v_), n, n);
            % E = diaginvic * conj(N);
            F = dlamJ.';
            G = F * B * conj(diagv);    %% F * B2
            dBtlam = sparse(1:n, 1:n, B.' * lam, n, n);
            H = dBtlam * conj(diagv);
            K = (F * A).';
            GG = G + G.';
            LL = GG - H - K;
            MM = D * (2*K - GG) * D;

            %% linear current term
            dYtlam = sparse(1:n, 1:n, Y.' * lam, n, n);
            Ivava = -dYtlam * diagv;
            Ivavm = -j * Ivava * D;

            %% + current from linear power term
            Ivava = Ivava + LL;
            Ivavm = Ivavm + 1j * LL.' * D;
            Ivmvm = MM;
        end

        function [Ivazr, Ivazi, Ivmzr, Ivmzi] = port_inj_current_hess_vz(obj, x_, lam, v_, z_, diaginvic, N, dlamJ)
            % Compute voltage/non-voltage-related terms of current injection Hessian.
            % ::
            %
            %   [Ivazr, Ivazi, Ivmzr, Ivmzi] = nme.port_inj_current_hess_vz(x_, lam)
            %   [...] = nme.port_inj_current_hess_vz(x_, lam, sysx)
            %   [...] = nme.port_inj_current_hess_vz(x_, lam, sysx, idx)
            %   [...] = nme.port_inj_current_hess_vz(x_, lam, v_, z_, diaginvic, N, dlamJ)
            %
            % Called by mp.form_ac.port_inj_current_hess to compute
            % voltage/non-voltage-related Hessian terms.

            if nargin < 8
                sysx = v_;
                if nargin < 5
                    idx = []
                else
                    idx = z_;
                end
                N = obj.get_params(idx, 'N');
                [v_, z_, vi_] = obj.x2vz(x_, sysx, idx);

                n  = length(v_);    %% number of all port voltages
                ni = length(vi_);   %% number of selected port voltages
                diaginvic = sparse(1:ni, 1:ni, 1 ./ conj(vi_), ni, ni);
                if isempty(idx)     %% all ports
                    dlamJ = sparse(1:n, 1:n, lam, n, n);
                else                %% selected ports
                    dlamJ = sparse(1:ni, idx, lam, ni, n);
                end
            else
                n  = length(v_);    %% number of all port voltages
                ni = length(lam);
            end

            %% intermediate terms
            D = sparse(1:n, 1:n, 1 ./ abs(v_), n, n);
            E = diaginvic * conj(N);
            NN = dlamJ.' * E;

            %% current from linear power term
            Ivazr = 1j * NN;
            Ivazi = NN;
            Ivmzr = -D * NN;
            Ivmzi = -1j * Ivmzr;
        end

        function [Sva, Svm] = port_inj_power_jac(obj, n, v_, Y, M, diagv, diagvi, diagIlincJ)
            % Compute voltage-related terms of power injection Jacobian.
            % ::
            %
            %   [Sva, Svm] = nme.port_inj_power_jac(...)
            %
            % Called by mp.form_ac.port_inj_power to compute
            % voltage-related Jacobian terms.

            %% intermediate terms
            A = diagvi * diagIlincJ;
%             B = diagvi * conj(Y);
%             C = B * conj(diagv);
            C = diagvi * conj(Y * diagv);
            D = sparse(1:n, 1:n, 1 ./ abs(v_), n, n);

%             %% linear power term
%             Svm = M * diagv;
%             Sva = 1j * Svm;
%             Svm = Svm * D;
% 
%             %% + power from linear current term
%             Sva = Sva + 1j * (A - C);
%             Svm = Svm + (A + C) * D;

            Svm = M * diagv + A;
            Sva = 1j * (Svm - C);
            Svm = (Svm + C) * D;
        end

        function [Svava, Svavm, Svmvm] = port_inj_power_hess_v(obj, x_, lam, v_, z_, diagvi, Y, M, diagIlincJ, dlamJ)
            % Compute voltage-related terms of power injection Hessian.
            % ::
            %
            %   [Svava, Svavm, Svmvm] = nme.port_inj_power_hess_v(x_, lam)
            %   [Svava, Svavm, Svmvm] = nme.port_inj_power_hess_v(x_, lam, sysx)
            %   [Svava, Svavm, Svmvm] = nme.port_inj_power_hess_v(x_, lam, sysx, idx)
            %   [...] = nme.port_inj_power_hess_v(x_, lam, v_, z_, diagvi, Y, M, diagIlincJ, dlamJ)
            %
            % Called by mp.form_ac.port_inj_power_hess to compute
            % voltage-related Hessian terms.

            if nargin < 10
                sysx = v_;
                if nargin < 5
                    idx = []
                else
                    idx = z_;
                end
                [Y, L, M, i] = obj.get_params(idx, {'Y', 'L', 'M', 'i'});
                [v_, z_, vi_] = obj.x2vz(x_, sysx, idx);

                %% compute linear current injections
                if isempty(z_)
                    Ilin = Y*v_ + i;
                else
                    Ilin = Y*v_ + L*z_ + i;
                end

                n  = length(v_);    %% number of all port voltages
                ni = length(vi_);   %% number of selected port voltages
                diagvi = sparse(1:ni, 1:ni, vi_, ni, ni);
                if isempty(idx)     %% all ports
                    diagIlincJ = sparse(1:n, 1:n, conj(Ilin), n, n);
                    dlamJ      = sparse(1:n, 1:n, lam, n, n);
                else                %% selected ports
                    diagIlincJ = sparse(1:ni, idx, conj(Ilin), ni, n);
                    dlamJ      = sparse(1:ni, idx, lam, ni, n);
                end
            else
                n  = length(v_);    %% number of all port voltages
                ni = length(lam);
            end

            %% intermediate terms
            diagv  = sparse(1:n, 1:n, v_, n, n);

            A = diagvi * diagIlincJ;
            B = diagvi * conj(Y);
            C = B * conj(diagv);
            D = sparse(1:n, 1:n, 1 ./ abs(v_), n, n);
            % E = diagvi * conj(L);
            F = dlamJ.';
            G = F * C;
            dBtlam = sparse(1:n, 1:n, B.' * lam, n, n);
            H = conj(diagv) * ((F*B).' - dBtlam);
            K = F * (C - A);

            %% linear power term
            dMtlam = sparse(1:n, 1:n, M.' * lam, n, n);
            Svava = -dMtlam * diagv;
            Svavm = -j * Svava * D;

            %% + power from linear current term
            Svava = Svava + H + K;
            Svavm = Svavm + 1j * (H - K).' * D;
            Svmvm = D * (G + G.') * D;
        end

        function [Svazr, Svazi, Svmzr, Svmzi] = port_inj_power_hess_vz(obj, x_, lam, v_, z_, diagvi, L, dlamJ)
            % Compute voltage/non-voltage-related terms of power injection Hessian.
            % ::
            %
            %   [Svazr, Svazi, Svmzr, Svmzi] = nme.port_inj_power_hess_vz(x_, lam)
            %   [...] = nme.port_inj_power_hess_vz(x_, lam, sysx)
            %   [...] = nme.port_inj_power_hess_vz(x_, lam, sysx, idx)
            %   [...] = nme.port_inj_power_hess_vz(x_, lam, v_, z_, diagvi, L, dlamJ)
            %
            % Called by mp.form_ac.port_inj_power_hess to compute
            % voltage/non-voltage-related Hessian terms.

            if nargin < 8
                sysx = v_;
                if nargin < 5
                    idx = []
                else
                    idx = z_;
                end
                L = obj.get_params(idx, 'L');
                [v_, z_, vi_] = obj.x2vz(x_, sysx, idx);

                n  = length(v_);    %% number of all port voltages
                ni = length(vi_);   %% number of selected port voltages
                diagvi = sparse(1:ni, 1:ni, vi_, ni, ni);
                if isempty(idx)     %% all ports
                    dlamJ = sparse(1:n, 1:n, lam, n, n);
                else                %% selected ports
                    dlamJ = sparse(1:ni, idx, lam, ni, n);
                end
            else
                n  = length(v_);    %% number of all port voltages
                ni = length(lam);
            end

            %% intermediate terms
            D = sparse(1:n, 1:n, 1 ./ abs(v_), n, n);
            E = diagvi * conj(L);
            LL = dlamJ.' * E;
            M = D * LL;

            %% power from linear current term
            Svazr = 1j * LL;
            Svazi = LL;
            Svmzr = M;
            Svmzi = -1j * M;
        end

        function [va, vm] = aux_data_va_vm(obj, ad)
            % Return voltage angles/magnitudes from auxiliary data.
            % ::
            %
            %   [va, vm] = nme.aux_data_va_vm(ad)
            %
            % Simply returns voltage data stored in ``ad.va`` and
            % ``ad.vm``.
            %
            % Input:
            %   ad (struct) : struct of auxiliary data
            %
            % Outputs:
            %   va (double) : vector of voltage angles corresponding to
            %       voltage information stored in auxiliary data
            %   vm (double) : vector of voltage magnitudes corresponding to
            %       voltage information stored in auxiliary data

            va = ad.va;
            vm = ad.vm;
        end
    end     %% methods
end         %% classdef
