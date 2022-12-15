classdef form_acc < mp.form_ac
%MP.FORM_ACC  MATPOWER Formulation class for AC cartesian voltage formulations
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
%       form_name() - returns string w/name of formulation ('AC-cartesian formulation')
%       form_tag() - returns string w/short label for formulation ('acc')

%   MATPOWER
%   Copyright (c) 2019-2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        function name = form_name(obj)
            name = 'AC-cartesian';
        end

        function tag = form_tag(obj)
            tag = 'acc';
        end

        function vtypes = model_vvars(obj)
            vtypes = {'vr', 'vi'};
        end

        function [Iu, Iw] = port_inj_current_jac(obj, ...
                n, v_, Y, M, invdiagvic, diagSlincJ)
            % [Iu, Iw] = obj.port_inj_current_jac(...)

            %% intermediate terms
            E = invdiagvic * (conj(M) - invdiagvic * diagSlincJ);

%             %% linear current term
%             Iu = Y;
%             Iw = 1j * Y;
% 
%             %% + current from linear power term
%             Iu = Iu + E;
%             Iw = Iw - 1j * E;

            Iu = Y + E;
            Iw = 1j * (Y - E);
        end

        function [Iuu, Iuw, Iww] = port_inj_current_hess_v(obj, x_, lam, v_, z_, diaginvic, Y, M, diagSlincJ, dlamJ)
            % [Iuu, Iuw, Iww] = obj.port_inj_current_hess_v(x_, lam)
            % [Iuu, Iuw, Iww] = obj.port_inj_current_hess_v(x_, lam, sysx)
            % [Iuu, Iuw, Iww] = obj.port_inj_current_hess_v(x_, lam, sysx, idx)
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
            % A = diaginvic * conj(M);
            % B = diaginvic * conj(N);
            % C = diaginvic * diaginvic * diagSlincJ;
            D = (diaginvic * dlamJ).';
            E = diaginvic * (conj(M) - diaginvic * diagSlincJ);
            % E = A - C;
            F = D * E;
            G = -(F.' + F);
            % H = -D * B;

            %% linear current term
            %% second derivatives all zero

            %% current from linear power term
            Iuu = G;
            Iuw = -1j * G;
            Iww = -G;
        end

        function [Iuzr, Iuzi, Iwzr, Iwzi] = port_inj_current_hess_vz(obj, x_, lam, v_, z_, diaginvic, N, dlamJ)
            % [Iuzr, Iuzi, Iwzr, Iwzi] = obj.port_inj_current_hess_vz(x_, lam)
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
            % A = diaginvic * conj(M);
            B = diaginvic * conj(N);
            % C = diaginvic * diaginvic * diagSlincJ;
            D = (diaginvic * dlamJ).';
            % E = diaginvic * (conj(M) - diaginvic * diagSlincJ);
            % E = A - C;
            % F = D * E;
            % G = -(F.' + F);
            H = -D * B;

            %% current from linear power term
            Iuzr = H;
            Iuzi = -1j * H;
            Iwzr = Iuzi;
            Iwzi = -H;
        end

        function [Su, Sw] = port_inj_power_jac(obj, ...
                n, v_, Y, M, diagv, diagvi, diagIlincJ)
            % [Su, Sw] = obj.port_inj_power_jac(...)

            %% intermediate terms
%             A = diagIlincJ;
            B = diagvi * conj(Y);

%             %% linear power term
%             Su = M;
%             Sw = 1j * M;
% 
%             %% + power from linear current term
%             Su = Su + A + B;
%             Sw = Sw + 1j * (A - B);

            A = M + diagIlincJ;
            Su = A + B;
            Sw = 1j * (A - B);
        end

        function [Suu, Suw, Sww] = port_inj_power_hess_v(obj, x_, lam, v_, z_, diagvi, Y, M, diagIlincJ, dlamJ)
            % [Suu, Suw, Sww] = obj.port_inj_power_hess_v(x_, lam)
            % [Suu, Suw, Sww] = obj.port_inj_power_hess_v(x_, lam, sysx)
            % [Suu, Suw, Sww] = obj.port_inj_power_hess_v(x_, lam, sysx, idx)
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
            % D = dlamJ.';
            E = dlamJ.' * conj(Y);  %% D * conj(Y);
            % F = E + E.';
            % G = 1j * (E - E.');

            %% linear power term
            %% second derivatives all zero

            %% power from linear current term
            Suu = E + E.';          %% F
            Suw = 1j * (E.' - E);   %% G.'
            Sww = Suu;
        end

        function [Suzr, Suzi, Swzr, Swzi] = port_inj_power_hess_vz(obj, x_, lam, v_, z_, diagvi, L, dlamJ)
            % [Suzr, Suzi, Swzr, Swzi] = obj.port_inj_power_hess_vz(x_, lam)
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
            D = dlamJ.';
            H = D * conj(L);

            %% power from linear current term
            Suzr = H;
            Suzi = -1j * H;
            Swzr = 1j * H;
            Swzi = H;
        end

        function [va, vm] = aux_data_va_vm(obj, ad)
            v_ = ad.vr + 1j * ad.vi;
            va = angle(v_);
            vm = abs(v_);
        end

        function [g, dg] = va_fcn(obj, xx, idx, lim)
            %% lim can be a vector value for equality constraint or
            %% upper bound on va, or a cell array with {vamin, vamax}
            %% for double bounds

            %% unpack data
            [vr, vi] = deal(xx{:});

            %% compute voltage angle mismatch
            if isempty(idx)
                va = angle(vr + 1j * vi);
            else
                va = angle(vr(idx) + 1j * vi(idx));
            end
            if iscell(lim)
                g = [ lim{1} - va;
                      va - lim{2}   ];
            else
                g = va - lim;
            end

            if nargout > 1
                %% compute partials of voltage angle w.r.t vr and vi
                nn = length(vr);
                if isempty(idx)
                    idx = 1:nn;
                end
                n = length(idx);
                vm2 = vr(idx).^2 + vi(idx).^2;
                dva_dvr = sparse(1:n, idx, -vi(idx) ./ vm2, n, nn);
                dva_dvi = sparse(1:n, idx,  vr(idx) ./ vm2, n, nn);
                if iscell(lim)
                    dg = [ -dva_dvr -dva_dvi;       %% va w.r.t vr, vi
                            dva_dvr  dva_dvi  ];    %% va w.r.t vr, vi
                else
                    dg = [dva_dvr dva_dvi];         %% va w.r.t vr, vi
                end
            end
        end

        function d2G = va_hess(obj, xx, lam, idx)
            %% unpack data
            [vr, vi] = deal(xx{:});
            nn = length(vr);

            %% evaluate Hessian of voltage angle function
            if isempty(idx)
                vvr = vr;
                vvi = vi;
                idx = 1:nn;
            else
                vvr = vr(idx);
                vvi = vi(idx);
            end
            vvr2 = vvr.^2;
            vvi2 = vvi.^2;
            vvm4 = (vvr2 + vvi2).^2;
            n = length(idx);
            if length(lam) == n     %% upper bound or equality
                lamvm4 = lam ./ vvm4;
            else                    %% doubly bounded (use lam_ub-lam_lb)
                lamvm4 = (lam(n+1:2*n) - lam(1:n)) ./ vvm4;
            end
            d2vref_rr = sparse(idx, idx, 2 * lamvm4 .*  vvr .* vvi,   nn, nn);
            d2vref_ri = sparse(idx, idx,     lamvm4 .* (vvi2 - vvr2), nn, nn);

            %% construct Hessian
            d2G = [ d2vref_rr  d2vref_ri;
                    d2vref_ri -d2vref_rr ];
        end

        function [g, dg] = vm2_fcn(obj, xx, idx, lim)
            %% lim can be a scalar value for equality constraint or
            %% upper bound on vm^2, or a cell array with {vm2min, vm2max}
            %% for double bounds

            %% unpack data
            [vr, vi] = deal(xx{:});

            %% compute voltage magnitude^2 mismatch
            if isempty(idx)
                vm2 = vr.^2 + vi.^2;
            else
                vm2 = vr(idx).^2 + vi(idx).^2;
            end
            if iscell(lim)
                g = [ lim{1} - vm2;
                      vm2 - lim{2}   ];
            else
                g = vm2 - lim;
            end

            if nargout > 1
                %% compute partials of voltage magnitude^2 w.r.t vr and vi
                nn = length(vr);
                if isempty(idx)
                    idx = 1:nn;
                end
                n = length(idx);
                dvm_dvr = sparse(1:n, idx, 2 * vr(idx), n, nn);
                dvm_dvi = sparse(1:n, idx, 2 * vi(idx), n, nn);
                if iscell(lim)
                    dg = [ -dvm_dvr -dvm_dvi;
                            dvm_dvr  dvm_dvi  ];    %% vm2 w.r.t vr, vi
                else
                    dg = [ dvm_dvr dvm_dvi ];       %% vm2 w.r.t vr, vi
                end
            end
        end

        function d2G = vm2_hess(obj, xx, lam, idx)
            %% unpack data
            [vr, vi] = deal(xx{:});
            nn = length(vr);

            %% evaluate Hessian of voltage magnitude^2 function
            if isempty(idx)
                idx = 1:nn;
            end
            n = length(idx);
            if length(lam) == n     %% upper bound or equality
                dlam = sparse(idx, idx, 2*lam, nn, nn);
            else                    %% doubly bounded (use lam_ub-lam_lb)
                dlam = sparse(idx, idx, 2*(lam(n+1:2*n) - lam(1:n)), nn, nn);
            end

            %% construct Hessian
            zz = sparse(nn, nn);
            d2G = [dlam zz; zz dlam];
        end
    end     %% methods
end         %% classdef
