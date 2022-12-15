classdef form_acp < mp.form_ac
%MP.FORM_ACP  MATPOWER Formulation class for AC polar voltage formulations
%   Each concrete Network Model Element class must inherit, at least
%   indirectly, from both MP.NM_ELEMENT and MP.FORM.
%
%   Subclass of MP.FORM_AC.
%   MP.FORM provides properties and methods related to the specific
%   formulation (e.g. DC version, AC polar power version, etc.)
%
%   Properties
%       (model parameters inherited from MP.FORM_AC)
%
%   Methods
%       form_name() - returns string w/name of formulation ('AC-polar formulation')
%       form_tag() - returns string w/short label for formulation ('acp')

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        function name = form_name(obj)
            name = 'AC-polar';
        end

        function tag = form_tag(obj)
            tag = 'acp';
        end

        function vtypes = model_vvars(obj)
            vtypes = {'va', 'vm'};
        end

        function [Iva, Ivm] = port_inj_current_jac(obj, ...
                n, v_, Y, M, invdiagvic, diagSlincJ)
            % [Iva, Ivm] = obj.port_inj_current_jac(...)

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
            % [Ivava, Ivavm, Ivmvm] = obj.port_inj_current_hess_v(x_, lam)
            % [Ivava, Ivavm, Ivmvm] = obj.port_inj_current_hess_v(x_, lam, sysx)
            % [Ivava, Ivavm, Ivmvm] = obj.port_inj_current_hess_v(x_, lam, sysx, idx)
            % [...] = obj.port_inj_current_hess_vz(x_, lam, v_, z_, diaginvic, Y, M, diagSlincJ, dlamJ)

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
            % [Ivazr, Ivazi, Ivmzr, Ivmzi] = obj.port_inj_current_hess_vz(x_, lam)
            % [...] = obj.port_inj_current_hess_vz(x_, lam, sysx)
            % [...] = obj.port_inj_current_hess_vz(x_, lam, sysx, idx)
            % [...] = obj.port_inj_current_hess_vz(x_, lam, v_, z_, diaginvic, N, dlamJ)

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

        function [Sva, Svm] = port_inj_power_jac(obj, ...
                n, v_, Y, M, diagv, diagvi, diagIlincJ)
            % [Sva, Svm] = obj.port_inj_power_jac(...)

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
            % [Svava, Svavm, Svmvm] = obj.port_inj_power_hess_v(x_, lam)
            % [Svava, Svavm, Svmvm] = obj.port_inj_power_hess_v(x_, lam, sysx)
            % [Svava, Svavm, Svmvm] = obj.port_inj_power_hess_v(x_, lam, sysx, idx)
            % [...] = obj.port_inj_power_hess_v(x_, lam, v_, z_, diagvi, Y, M, diagIlincJ, dlamJ)

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
            % [Svazr, Svazi, Svmzr, Svmzi] = obj.port_inj_power_hess_vz(x_, lam)
            % [...] = obj.port_inj_power_hess_vz(x_, lam, sysx)
            % [...] = obj.port_inj_power_hess_vz(x_, lam, sysx, idx)
            % [...] = obj.port_inj_power_hess_vz(x_, lam, v_, z_, diagvi, L, dlamJ)

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
            va = ad.va;
            vm = ad.vm;
        end
    end     %% methods
end         %% classdef
